#ifndef __TRIANGLE_CPP_H__
#define __TRIANGLE_CPP_H__

#include <vector>
#include <map>

#include "./opencl/triangle.h"

inline Triangle make_triangle(int v[]) {
  return *(Triangle*)(v);
}

//----------------------------------------
// Operators
//----------------------------------------

inline bool operator==(const Triangle& lhs, const Triangle& rhs) {
  for (int i = 0; i < 3; ++i) {
    int li = -1;
    for (int j = 0; li == -1 && j < 3; ++j) {
      if (lhs.s[i] == rhs.s[j]) li = j;
    }
    if (li == -1) return false;
  }
  return true;
}

inline bool operator!=(const Triangle& lhs, const Triangle& rhs) {
  return !(lhs == rhs);
}

inline bool operator<(const Triangle& lhs, const Triangle& rhs) { 
  int sorted[3] = { lhs.s[0], lhs.s[1], lhs.s[2] };
  int rsorted[3] = { rhs.s[0], rhs.s[1], rhs.s[2] };
  std::sort(sorted, sorted+3);
  std::sort(rsorted, rsorted+3);
  for (int i = 0; i < 3; ++i) {
    if (sorted[i] != rsorted[i]) {
      return sorted[i] < rsorted[i];
    }
  }
  return false;
}

inline std::ostream& operator<<(std::ostream& out, const Triangle& t) {
  out << t.s[0] << " " << t.s[1] << " " << t.s[2];
  return out;
}

inline std::istream& operator>>(std::istream& in, Triangle& t) {
  int i, j, k;
  in >> i >> j >> k;
  t = make_triangle(i, j, k);
  return in;
}

template <typename Out_iter>
void get_edges(const Triangle& t, Out_iter out) {
  for (int i = 0; i < 3; ++i) {
    *out++ = make_edge(t.s[i], t.s[(i+1)%3]);
  }
}

template <typename Iter>
std::vector<std::vector<Triangle> > build_v2t(Iter begin, Iter end)
{
  int m = -1;
  for (Iter it = begin; it != end; ++it) {
    const Triangle& t = *it;
    for (int i = 0; i < 3; ++i) {
      m = std::max(m, t.s[i]);
    }
  }
  std::vector<std::vector<Triangle> > ret(m+1);
  for (Iter it = begin; it != end; ++it) {
    const Triangle& t = *it;
    for (int i = 0; i < 3; ++i) {
      ret[t.s[i]].push_back(t);
    }
  }
  return ret;
}

template <typename Iter>
std::map<Edge, std::vector<Triangle> > build_e2t(Iter begin, Iter end)
{
  std::map<Edge, std::vector<Triangle> > ret;
  for (Iter it = begin; it != end; ++it) {
    const Triangle& t = *it;
    for (int i = 0; i < 3; ++i) {
      ret[triangle_opposite(t.s[i], t)].push_back(t);
    }
  }
  return ret;
}

#endif
