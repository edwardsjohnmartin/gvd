#ifndef __EDGE_CPP_H__
#define __EDGE_CPP_H__

#include "./opencl/triangle.h"

inline std::ostream& operator<<(std::ostream& out, const Edge& e) {
  out << e.s[0] << " " << e.s[1];
  return out;
}

inline std::istream& operator>>(std::istream& in, Edge& e) {
  int a, b;
  in >> a >> b;
  e = make_edge(a, b);
  return in;
}

inline bool operator==(const Edge& lhs, const Edge& rhs) {
  for (int i = 0; i < 2; ++i) {
    int li = -1;
    for (int j = 0; li == -1 && j < 2; ++j) {
      if (lhs.s[i] == rhs.s[j]) li = j;
    }
    if (li == -1) return false;
  }
  return true;
}

inline bool operator!=(const Edge& lhs, const Edge& rhs) {
  return !(lhs == rhs);
}

inline bool operator<(const Edge& lhs, const Edge& rhs) {
  int sorted[2] = { rhs.s[0], rhs.s[1] };
  int rsorted[2] = { lhs.s[0], lhs.s[1] };
  if (sorted[0] > sorted[1]) {
    std::swap(sorted[0], sorted[1]);
  }
  if (rsorted[0] > rsorted[1]) {
    std::swap(rsorted[0], rsorted[1]);
  }
  for (int i = 0; i < 2; ++i) {
    if (sorted[i] != rsorted[i]) {
      return sorted[i] < rsorted[i];
    }
  }
  return false;
}

#endif
