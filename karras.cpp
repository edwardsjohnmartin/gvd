#include "./karras.h"

#include <iostream>
#include <algorithm>
#include <memory>

#include "./bb.h"
#include "./gpu.h"

using std::cout;
using std::endl;
using std::vector;
using std::shared_ptr;

// kWidth is the quantized width in one dimension.
// Number of bits per dimension is kNumBits = log(kWidth).
// Total number of bits for morton code is kNumBits * DIM.
// static const int kWidth = 8;
// static const int kNumBits = 3;

namespace Karras {

int xyz2z(intn p, const Resln& resln) {
  int ret = 0;
  for (int i = 0; i < resln.bits; ++i) {
    for (int j = 0; j < DIM; ++j) {
      if (p.s[j] & (1<<i))
        ret |= 1 << (i*DIM+j);
    }
  }
  return ret;
}

intn z2xyz(const int z, const Resln& resln) {
  intn p = make_intn(0);
  for (int i = 0; i < resln.bits; ++i) {
    for (int j = 0; j < DIM; ++j) {
      if (z & (1 << (i*DIM+j)))
        p.s[j] |= (1<<i);
    }
  }
  return p;
}

int sign(const int i) {
  return (i<0) ? -1 : ((i>0) ? +1 : 0);
}

// Longest common prefix
//
// Suppose mbits = 6, then morton code is
//   ______
// 00011010
//
// Suppose length = 3, then lcp (masked) is
//   ___
// 00011000
//
// Now shift, and lcp is
//      ___
// 00000011
int compute_lcp(const int value, const int length, const Resln& resln) {
  int mask = 0;
  for (int i = 0; i < length; ++i) {
    mask |= (1 << (resln.mbits - 1 - i));
  }
  const int lcp = (value & mask) >> (resln.mbits - length);
  return lcp;
}

// Longest common prefix, denoted \delta in karras2014
int compute_lcp_length(const int a, const int b, const Resln& resln) {
  for (int i = resln.mbits-1; i >= 0; --i) {
    const int mask = 1 << i;
    if ((a & mask) != (b & mask)) {
      return resln.mbits - i - 1;
    }
  }
  return resln.mbits;
}

class LcpLength {
 public:
  LcpLength(const vector<int>& mpoints, const Resln& resln)
      : _mpoints(mpoints), _resln(resln) {}

  int operator()(const int i, const int j) const {
    if (i < 0 || i >= _mpoints.size() || 
        j < 0 || j >= _mpoints.size()) {
      throw logic_error("Illegal indices into lcp_length");
    }
    return compute_lcp_length(_mpoints[i], _mpoints[j], _resln);
  }

 private:
  const vector<int>& _mpoints;
  const Resln& _resln;
};

struct BrtNode {
  // If lcp is 10011 and DIM == 2 then the last bit is dropped
  // and return is [0b10, 0b01] = [2, 1].
  vector<int> oct_nodes() const {
    static const int mask = (DIM == 2) ? 3 : 7;
    if (DIM > 3)
      throw logic_error("BrtNode::oct_nodes not yet supported for D>3");
    const int bias = lcp_length % DIM;
    const int n = lcp_length / DIM;
    vector<int> ret(n);
    for (int i = 0; i < n; ++i) {
      const int offset = DIM * (n-i-1) + bias;
      ret[i] = (lcp >> offset) & mask;
    }
    return ret;
  }

  // left child (right child = left+1)
  int left;
  // Whether the left (resp. right) child is a leaf or not
  bool left_leaf, right_leaf;
  // The longest common prefix
  int lcp;
  // Number of bits in the longest common prefix
  int lcp_length;

  // Secondary - computed in a second pass
  int parent;
};

vector<intn> Quantize(const vector<floatn>& points, const Resln& resln) {
  if (points.empty())
    return vector<intn>();

  BoundingBox<floatn> bb;
  for (const floatn& p : points) {
    bb(p);
  }
  const float dwidth = bb.max_size();
  if (dwidth == 0) {
    vector<intn> ret;
    ret.push_back(make_intn(0));
    return ret;
  }

  // Quantize points to integers
  vector<intn> qpoints(points.size());
  for (int i = 0; i < points.size(); ++i) {
    const floatn& p = points[i];
    intn q = make_intn(0);
    for (int k = 0; k < DIM; ++k) {
      const double d =
          (resln.width-1) * ((p.s[k] - bb.min().s[k]) / dwidth);
      const int v = static_cast<int>(d+0.5);
      if (v < 0) {
        cerr << "Coordinate in dimension " << k << " is less than zero.  d = "
             << d << " v = " << v << endl;
        cerr << "  p[k] = " << p.s[k]
             << " bb.min()[k] = " << bb.min().s[k] << endl;
        cerr << "  dwidth = " << dwidth << " kwidth = " << resln.width << endl;
        throw logic_error("bad coordinate");
      }
      q.s[k] = v;
    }
    qpoints[i] = q;
  }
  
  return qpoints;
}

vector<OctNode> BuildOctree(
    const vector<intn>& points, const Resln& resln, const bool verbose) {
  if (points.empty())
    throw logic_error("Zero points not supported");
  
  vector<int> mpoints(points.size());
  for (int i = 0; i < points.size(); ++i) {
    mpoints[i] = xyz2z(points[i], resln);
  }

  sort(mpoints.begin(), mpoints.end());

  // Make sure points are unique
  std::vector<int>::iterator it;
  it = std::unique (mpoints.begin(), mpoints.end());
  mpoints.resize(std::distance(mpoints.begin(),it));

  // // Send mpoints to gpu
  // if (o.gpu) {
  //   oct::Gpu gpu;
  //   gpu.CreateMPoints(mpoints.size());
  // }

  const int n = mpoints.size();
  const LcpLength lcp_length(mpoints, resln);
  vector<BrtNode> I(n-1);
  vector<BrtNode> L(n);
  for (int i = 0; i < n-1; ++i) {
    // Determine direction of the range (+1 or -1)
    const int d = (i==0) ? 1 : sign(lcp_length(i, i+1) - lcp_length(i, i-1));
    // Compute upper bound for the length of the range
    int l;
    if (i == 0) {
      l = n-1;
    } else {
      const int lcp_min = lcp_length(i, i-d);
      int l_max = 2;
      while (i+l_max*d >= 0 && i+l_max*d <= n-1 &&
             lcp_length(i, i+l_max*d) > lcp_min) {
        l_max = l_max << 1;
      }
      // Find the other end using binary search.
      // In some cases, the search can go right off the end of the array.
      // l_max likely is beyond the end of the array, but we need it to be
      // since it's a power of 2. So define a max length that we call l_cutoff.
      const int l_cutoff = (d==-1) ? i : n - i - 1;
      l = 0;
      for (int t = l_max / 2; t >= 1; t /= 2) {
        if (l + t <= l_cutoff) {
          if (lcp_length(i, i+(l+t)*d) > lcp_min) {
            l = l + t;
          }
        }
      }
    }

    // j is the index of the other end of the range. In other words,
    // range = [i, j] or range = [j, i].
    const int j = i + l * d;
    // Find the split position using binary search
    const int lcp_node = lcp_length(i, j);

    const int s_cutoff = (d==-1) ? i - 1 : n - i - 2;
    int s = 0;
    for (int den = 2; den < 2*l; den *= 2) {
      const int t = static_cast<int>(ceil(l/(float)den));
      if (s + t <= s_cutoff) {
        if (lcp_length(i, i+(s+t)*d) > lcp_node) {
          s = s + t;
        }
      }
    }
    const int split = i + s * d + min(d, 0);
    // Output child pointers
    I[i].left = split;
    I[i].left_leaf = (min(i, j) == split);
    I[i].right_leaf = (max(i, j) == split+1);
    I[i].lcp = compute_lcp(mpoints[i], lcp_node, resln);
    I[i].lcp_length = lcp_node;
  }

  // Set parents
  if (n > 0)
    I[0].parent = -1;
  for (int i = 0; i < n-1; ++i) {
    const int left = I[i].left;
    const int right = left+1;
    if (!I[i].left_leaf) {
      I[left].parent = i;
    }
    if (!I[i].right_leaf) {
      I[right].parent = i;
    }
  }

  // Debug output
  if (verbose) {
    cout << endl;
    for (int i = 0; i < n-1; ++i) {
      cout << i << ": left = " << I[i].left << (I[i].left_leaf ? "L" : "I")
           << " right = " << I[i].left+1 << (I[i].right_leaf ? "L" : "I")
           << " lcp = " << I[i].lcp
           << " oct nodes: (";
      vector<int> onodes = I[i].oct_nodes();
      for (const int onode : onodes) {
        cout << onode << " ";
      }
      cout << ") lcp_length = " << I[i].lcp_length
           << " parent = " << I[i].parent
           << endl;
    }
  }

  // Determine number of octree nodes necessary
  // First pass - initialize temporary array
  vector<int> local_splits(n-1, 0);
  if (n > 0)
    local_splits[0] = 1;
  for (int i = 0; i < n-1; ++i) {
    const int local = I[i].lcp_length / DIM;
    const int left = I[i].left;
    const int right = left+1;
    if (!I[i].left_leaf) {
      local_splits[left] = I[left].lcp_length/DIM - local;
    }
    if (!I[i].right_leaf) {
      local_splits[right] = I[right].lcp_length/DIM - local;
    }
  }
  // Second pass - calculate prefix sums
  vector<int> prefix_sums(local_splits.size()+1);
  std::copy(local_splits.begin(), local_splits.end(), prefix_sums.begin()+1);
  prefix_sums[0] = 0;
  for (int i = 1; i < prefix_sums.size(); ++i) {
    prefix_sums[i] = prefix_sums[i-1] + prefix_sums[i];
  }

  if (verbose) {
    for (int i = 0; i < prefix_sums.size(); ++i) {
      cout << prefix_sums[i] << " ";
    }
    cout << endl;
  }

  const int splits = prefix_sums.back();
  if (verbose) {
    cout << "# octree splits = " << splits << endl;
  }

  // Set parent for each octree node
  vector<OctNode> octree(splits);
  for (int brt_i = 1; brt_i < n-1; ++brt_i) {
    if (local_splits[brt_i] > 0) {
      // m = number of local splits
      const int m = local_splits[brt_i];
      // onodes = vector of octree indices \in [0, 2^DIM]
      const vector<int> onodes = I[brt_i].oct_nodes();
      // current octree node index
      int oct_i = prefix_sums[brt_i];
      for (int j = 0; j < m-1; ++j) {
        const int oct_parent = oct_i+1;
        const int onode = onodes[onodes.size() - j - 1];
        octree[oct_parent].set_child(onode, oct_i);
        oct_i = oct_parent;
      }
      int brt_parent = I[brt_i].parent;
      while (local_splits[brt_parent] == 0) {
        brt_parent = I[brt_parent].parent;
      }
      const int oct_parent = prefix_sums[brt_parent];
      if (brt_parent < 0 || brt_parent >= prefix_sums.size()) {
        throw logic_error("error 0");
      }
      if (oct_parent < 0 || oct_parent >= octree.size()) {
        throw logic_error("error 1");
      }
      if (int(onodes.size()) - m < 0) {
        throw logic_error("error 2");
      }
      if (onodes[onodes.size() - m] < 0 ||
          onodes[onodes.size() - m] >= 4) {
        throw logic_error("error 3");
      }
      octree[oct_parent].set_child(onodes[onodes.size() - m], oct_i);
    }
  }

  if (verbose) {
    for (int i = 0; i < octree.size(); ++i) {
      cout << i << ": ";
      for (int j = 0; j < 4; ++j) {
        cout << octree[i][j] << " ";
      }
      cout << endl;
    }
    OutputOctree(octree);
  }

  return octree;
}

// Debug output
void OutputOctreeNode(
    const int node, const std::vector<OctNode>& octree, vector<int> path) {
  for (int i = 0; i < 4; ++i) {
    vector<int> p(path);
    p.push_back(i);

    for (int i : p) {
      cout << i;
    }
    cout << endl;

    if (!octree[node].is_leaf(i))
      OutputOctreeNode(octree[node][i], octree, p);
  }
}

void OutputOctree(const std::vector<OctNode>& octree) {
  if (!octree.empty()) {
    cout << endl;
    vector<int> p;
    p.push_back(0);
    OutputOctreeNode(0, octree, p);
  }
}

} // namespace
