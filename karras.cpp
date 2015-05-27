#include "./karras.h"

#include <iostream>
#include <algorithm>
#include <memory>

#include "./bb.h"
#include "./opencl/defs.h"

using std::cout;
using std::endl;
using std::vector;
using std::shared_ptr;

// static const int kNumBits = oct::kMaxLevel;
// static const int kNumBits = 8;

// kWidth is the quantized width in one dimension.
// Number of bits per dimension is kNumBits = log(kWidth).
// Total number of bits for morton code is kNumBits * DIM.
static const int kWidth = 8;
static const int kNumBits = 3;

int xyz2z(intn p) {
  int ret = 0;
  for (int i = 0; i < kNumBits; ++i) {
    for (int j = 0; j < DIM; ++j) {
      if (p.s[j] & (1<<i))
        ret |= 1 << (i*DIM+j);
    }
  }
  return ret;
}

intn z2xyz(const int z) {//, intn& p) {
  intn p = make_intn(0);
  for (int i = 0; i < kNumBits; ++i) {
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

// Longest common prefix, denoted \delta in karras2014
int compute_lcp(const int a, const int b) {
  for (int i = kNumBits*DIM-1; i >= 0; --i) {
    const int mask = 1 << i;
    if ((a & mask) != (b & mask)) {
      return kNumBits*DIM - i - 1;
    }
  }
  return kNumBits*DIM;
}

class LCP {
 public:
  LCP(const vector<int>& mpoints) : _mpoints(mpoints) {}

  int operator()(const int i, const int j) const {
    return compute_lcp(_mpoints[i], _mpoints[j]);
  }

 private:
  const vector<int>& _mpoints;
};

struct BrtNode {
  int left; // right = left+1
  bool left_leaf, right_leaf;
  int lcp;
};

struct OctNode {
  OctNode() : children(0) {}
  ~OctNode() {
    if (children)
      delete children;
  }
  bool is_leaf() const { return !children; }
  int& operator[](const int i) {
    if (!children)
      children = new int[1<<DIM];
    return children[i];
  }
  const int& operator[](const int i) const {
    if (!children) throw logic_error("leaf can't get children");
    return children[i];
  }
  int* children;
};

void ktest(const vector<floatn>& points) {
  BoundingBox<floatn> bb;
  for (const floatn& p : points) {
    bb(p);
  }
  const float dwidth = bb.max_size();
  if (dwidth == 0) return;

  // Quantize points to integers and morton codes
  vector<intn> qpoints(points.size());
  // vector<int> mpoints(points.size());
  for (int i = 0; i < points.size(); ++i) {
    const floatn& p = points[i];
    intn q = make_intn(0);
    for (int k = 0; k < DIM; ++k) {
      const double d =
          (kWidth-1) * ((p.s[k] - bb.min().s[k]) / dwidth);
      const int v = static_cast<int>(d+0.5);
      if (v < 0) {
        cerr << "Coordinate in dimension " << k << " is less than zero.  d = "
             << d << " v = " << v << endl;
        cerr << "  p[k] = " << p.s[k]
             << " bb.min()[k] = " << bb.min().s[k] << endl;
        cerr << "  dwidth = " << dwidth << " kwidth = " << kWidth << endl;
        throw logic_error("bad coordinate");
      }
      q.s[k] = v;
    }
    qpoints[i] = q;
    // mpoints[i] = xyz2z(q);
  }
  
  ktest(qpoints);
}

void ktest(const vector<intn>& points) {
  vector<int> mpoints(points.size());
  for (int i = 0; i < points.size(); ++i) {
    mpoints[i] = xyz2z(points[i]);
  }

  sort(mpoints.begin(), mpoints.end());

  // cout << endl;
  // for (const int m : mpoints) {
  //   cout << m << endl;
  // }

  const int n = points.size();
  const LCP lcp(mpoints);
  vector<BrtNode> I(n-1);
  vector<BrtNode> L(n);
  for (int i = 0; i < n-1; ++i) {
    // Determine direction of the range (+1 or -1)
    const int d = (i==0) ? 1 : sign(lcp(i, i+1) - lcp(i, i-1));
    // Compute upper bound for the length of the range
    int l;
    if (i == 0) {
      l = n-1;
    } else {
      const int lcp_min = lcp(i, i-d);
      int l_max = 2;
      while (i+l_max*d > 0 && i+l_max*d <= n-1 && lcp(i, i+l_max*d) > lcp_min) {
        l_max = l_max << 1;
      }
      // Find the other end using binary search
      l = 0;
      for (int t = l_max / 2; t >= 1; t /= 2) {
        if (lcp(i, i+(l+t)*d) > lcp_min) {
          l = l + t;
        }
      }
    }
    const int j = i + l * d;
    // Find the split position using binary search
    const int lcp_node = lcp(i, j);
    int s = 0;
    for (int den = 2; den < 2*l; den *= 2) {
      const int t = static_cast<int>(ceil(l/(float)den));
      if (lcp(i, i+(s+t)*d) > lcp_node) {
        s = s + t;
      }
    }
    const int split = i + s * d + min(d, 0);
    // Output child pointers
    I[i].left = split;
    I[i].left_leaf = (min(i, j) == split);
    I[i].right_leaf = (max(i, j) == split+1);
    I[i].lcp = lcp_node;
  }

  // Debug output
  for (int i = 0; i < n-1; ++i) {
    cout << i << ": left = " << I[i].left << (I[i].left_leaf ? "L" : "I")
         << " right = " << I[i].left+1 << (I[i].right_leaf ? "L" : "I")
         << " lcp = " << I[i].lcp
         << endl;
  }

  // Determine number of octree nodes necessary
  int sum = 0;
  for (int i = 0; i < n-1; ++i) {
    const int local = I[i].lcp / DIM;
    const int left = I[i].left;
    const int right = left+1;
    const int left_lcp = I[i].left_leaf ? 0 : I[left].lcp/DIM;
    const int right_lcp = I[i].right_leaf ? 0 : I[right].lcp/DIM;
    if (left_lcp > 0 || right_lcp > 0) {
      sum += max(left_lcp, right_lcp);
    }
  }
  cout << "# octree nodes = " << sum << endl;

  shared_ptr<OctNode> octree(new OctNode[sum]);
  
}
