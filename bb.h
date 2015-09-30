/*******************************************************
 ** Generalized Voronoi Diagram Project               **
 ** Copyright (c) 2015 John Martin Edwards            **
 ** Scientific Computing and Imaging Institute        **
 ** 72 S Central Campus Drive, Room 3750              **
 ** Salt Lake City, UT 84112                          **
 **                                                   **
 ** For information about this project contact        **
 ** John Edwards at                                   **
 **    edwardsjohnmartin@gmail.com                    **
 ** or visit                                          **
 **    sci.utah.edu/~jedwards/research/gvd/index.html **
 *******************************************************/

#ifndef _123_BOUNDING_BOX_H_
#define _123_BOUNDING_BOX_H_

#include <algorithm>

#include "./opencl/vec.h"

using namespace std;

template <typename Vec>
class VecTypes {
};
template<> class VecTypes<int2> {
 public:
  typedef int NumType;
  static const int D = 2;
};
template<> class VecTypes<float2> {
 public:
  typedef float NumType;
  static const int D = 2;
};
template<> class VecTypes<double2> {
 public:
  typedef double NumType;
  static const int D = 2;
};
template<> class VecTypes<int3> {
 public:
  typedef int NumType;
  static const int D = 3;
};
template<> class VecTypes<float3> {
 public:
  typedef float NumType;
  static const int D = 3;
};
template<> class VecTypes<double3> {
 public:
  typedef double NumType;
  static const int D = 3;
};

// template <typename NumType, int D>
template <typename Point>
class BoundingBox {
 public:
  // typedef MyVec<NumType, D> Point;
  typedef typename VecTypes<Point>::NumType NumType;
  static const int D = VecTypes<Point>::D;

 public:
  BoundingBox() : _initialized(false) {}
  BoundingBox(const Point& a, const Point& b) : _initialized(false) {
    (*this)(a);
    (*this)(b);
  }

  Point center() const {
    return Point((_min+_max)/2);
  }

  Point size() const { return _max - _min; }
  Point min() const { return _min; }
  Point max() const { return _max; }
  NumType max_size() const {
    const Point s = size();
    NumType m = s.s[0];
    for (int i = 1; i < D; ++i) {
      m = (s.s[i]>m)?s.s[i]:m;
    }
    return m;
  }

  bool in_open(const Point& p, const NumType epsilon = 0) const {
    for (int i = 0; i < D; ++i) {
      if (p.s[i] <= _min.s[i]-epsilon ||
          p.s[i] >= _max.s[i]+epsilon) return false;
    }
    return true;
  }

  bool in_half_open(const Point& p, const NumType epsilon = 0) const {
    for (int i = 0; i < D; ++i) {
      if (p.s[i] < _min.s[i]-epsilon ||
          p.s[i] >= _max.s[i]+epsilon) return false;
    }
    return true;
  }

  bool in_closed(const Point& p, const NumType epsilon = 0) const {
    for (int i = 0; i < D; ++i) {
      if (p.s[i] < _min.s[i]-epsilon ||
          p.s[i] > _max.s[i]+epsilon) return false;
    }
    return true;
  }

  // Returns the smallest square bounding box that contains
  // this and has identical origin.
  BoundingBox<Point> Square() const {
    const Point pwidth = max() - min();
    NumType dwidth = pwidth.s[0];
    for (int i = 1; i < D; ++i) {
      dwidth = std::max(dwidth, pwidth.s[i]);
    }

    if (dwidth == 0) return *this;

    return BoundingBox<Point>(min(), min()+dwidth);
  }

  // Returns the smallest square bounding box that contains
  // this and is centered in the same place.
  BoundingBox<Point> CenteredSquare() const {
    const Point pwidth = max() - min();
    NumType dwidth = pwidth.s[0];
    for (int i = 1; i < D; ++i) {
      dwidth = std::max(dwidth, pwidth.s[i]);
    }

    if (dwidth == 0) return *this;

    Point origin = make_floatn(0);
    for (int i = 0; i < D; ++i) {
      origin.s[i] = min().s[i] - (dwidth-size().s[i])/2;
    }
    return BoundingBox<Point>(origin, origin+dwidth);
  }

  BoundingBox<Point> Scale(const double& f) const {
    return BoundingBox<Point>(min(), min()+size()*f);
  }

  BoundingBox<Point> ScaleCentered(const double& f) const {
    const Point s = size()*f;
    const Point c = center();
    const Point minp = c - s/2;
    const Point maxp = c + s/2;
    return BoundingBox<Point>(minp, maxp);
  }

  bool IsSquare() const {
    const Point s = size();
    const NumType a = s.s[0];
    for (int i = 1; i < D; ++i) {
      if (s.s[i] != a) return false;
    }
    return true;
  }

  void operator()(const Point& v) {
    if (_initialized) {
      _min = vec_min(_min, v);
      _max = vec_max(_max, v);
    } else {
      _min = v;
      _max = v;
      _initialized = true;
    }
  }

  void operator+=(const BoundingBox& rhs) {
    if (_initialized) {
      (*this)(rhs._min);
      (*this)(rhs._max);
    } else {
      *this = rhs;
    }
  }

  // Resize by f in each dimension
  BoundingBox operator*(const NumType f) const {
    if (_initialized) {
      const Point s = size() * f;
      const Point c = center();
      Point mi = c - s/2;
      Point ma = c + s/2;
      BoundingBox ret;
      ret(mi);
      ret(ma);
      return ret;
    } else {
      return *this;
    }
  }

 private:
  bool _initialized;
  Point _min;
  Point _max;
};

template <typename Point>
inline ostream & operator<<(ostream & out, const BoundingBox<Point>& b) {
  out << "[box min=" << b.min()
      << " max=" << b.max()
      << " center=" << b.center() << "]";
  return out;
}

typedef BoundingBox<float3> BoundingBox3f;

#endif
