#ifndef __RESLN_H__
#define __RESLN_H__

#include "./bigint/BigUnsigned.hh"
#include "./bigint/BigIntegerUtils.hh"

namespace Karras {

// typedef unsigned int Morton;
// typedef unsigned long Morton;
typedef BigUnsigned Morton;

// Stores resolution and octree height values
struct Resln {
 public:
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

  friend std::ostream& operator<<(std::ostream& out, const Resln& resln);
  friend std::istream& operator>>(std::istream& in, Resln& resln);

  // width is the quantized width in one dimension.
  int width;
  int volume;
  // Number of bits per dimension is bits = log(width).
  int bits;
  // Total number of bits for morton code is mbits = bits * DIM.
  int mbits;
};

inline std::ostream& operator<<(std::ostream& out, const Resln& resln) {
  out << resln.width << " " << resln.volume << " "
      << resln.bits << " " << resln.mbits;
  return out;
}

inline std::istream& operator>>(std::istream& in, Resln& resln) {
  in >> resln.width >> resln.volume >> resln.bits >> resln.mbits;
  return in;
}

} // namespace

#endif
