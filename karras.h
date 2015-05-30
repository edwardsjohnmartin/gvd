#ifndef __KARRAS_H__
#define __KARRAS_H__

#include <vector>
#include <stdexcept>

#include "./opencl/defs.h"
#include "./opencl/vec.h"

namespace Karras {

struct OctNode {
  OctNode() : children(0) {}
  ~OctNode() {
    if (children)
      delete children;
  }
  bool is_leaf() const { return !children; }
  bool is_leaf(const int i) const {
    EnsureChildren();
    return children[i] == -1;
  }
  int& operator[](const int i) {
    EnsureChildren();
    return children[i];
  }
  const int& operator[](const int i) const {
    if (!children) throw std::logic_error("leaf can't get children");
    return children[i];
  }

 private:
  void EnsureChildren() const {
    if (!children) {
      children = new int[1<<DIM];
      std::fill(children, children + (1<<DIM), -1);
    }
  }

 public:
  mutable int* children;
};

// Stores resolution and octree height values
struct Resln {
  Resln()
      : width(8), volume(64), bits(3), mbits(3*DIM) {}
  Resln(const int width_)
      : width(width_) {
    if (width == 0) {
      throw std::logic_error("No support for width of 0");
    }
    volume = width;
    for (int i = 1; i < DIM; ++i) {
      volume *= width;
    }
    bits = 0;
    int w = width;
    while (!(w & 1)) {
      ++bits;
      w = w >> 1;
    }
    mbits = bits * DIM;
  }

  // width is the quantized width in one dimension.
  int width;
  int volume;
  // Number of bits per dimension is bits = log(width).
  int bits;
  // Total number of bits for morton code is bits * DIM.
  int mbits;
};

int xyz2z(intn p, const Resln& r);
intn z2xyz(const int z, const Resln& r);

std::vector<intn> Quantize(const std::vector<floatn>& points, const Resln& r);

std::vector<OctNode> BuildOctree(
    const std::vector<intn>& opoints, const Resln& r, const bool verbose=false);

void OutputOctree(const std::vector<OctNode>& octree);

} // namespace

inline std::ostream& operator<<(std::ostream& out, const Karras::Resln& resln) {
  out << "width=" << resln.width << ", volume=" << resln.volume
      << ", bits=" << resln.bits
      << ", mbits=" << resln.mbits;
  return out;
}

#endif
