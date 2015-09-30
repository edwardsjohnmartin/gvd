#ifndef __SEGMENT_H__
#define __SEGMENT_H__

#include "./vec.h"

class float_seg {
 public:
  float_seg() {}
  float_seg(const floatn& a_, const floatn& b_)
      : _a(a_), _b(b_) {}

  const floatn& operator[](int i) const { return (i==0) ? _a : _b; }
  floatn& operator[](int i) { return (i==0) ? _a : _b; }
  const floatn& a() const { return _a; }
  const floatn& b() const { return _b; }
  floatn& a() { return _a; }
  floatn& b() { return _b; }
  floatn vector() const { return _b - _a; }
  const floatn& p() const { return _a; }
  floatn unit() const { return vector() / length(); }
  float length() const { return ::length(vector()); }
  float length2() const { return ::length2(vector()); }
  void reverse() { std::swap(_a, _b); }
  bool is_degenerate() const { return _a == _b; }

  friend std::ostream& operator<<(std::ostream& out, const float_seg& s);
  friend std::istream& operator>>(std::istream& in, float_seg& s);

 private:
  floatn _a;
  floatn _b;
};

inline std::ostream& operator<<(std::ostream& out, const float_seg& s) {
  out << s.a() << " " << s.b();
  return out;
}

inline std::istream& operator>>(std::istream& in, float_seg& s) {
  in >> s._a >> s._b;
  return in;
}

#endif
