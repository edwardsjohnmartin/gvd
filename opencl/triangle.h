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

#ifndef __TRIANGLE_H__
#define __TRIANGLE_H__

#include "./edge.h"

typedef struct {
  int s[3];
} Triangle;

//------------------------------------------------------------
// Functions
//------------------------------------------------------------

inline Triangle make_triangle(int v0, int v1, int v2) {
  const int v[3] = { v0, v1, v2 };
  return *(Triangle*)(v);
}

inline Triangle triangle_inverted(const Triangle t) {
  return make_triangle(t.s[2], t.s[1], t.s[0]);
}

// Given a local vertex index, return the
// global vertex index.
inline int triangle_l2g(int i, const Triangle t) {
  return t.s[i];
}

// Given a global vertex index, return the
// local vertex index.  Returns -1 if vi isn't
// a vertex.
inline int triangle_g2l(int vi, const Triangle t) {
  for (int i = 0; i < 3; ++i) {
    if (t.s[i] == vi) {
      return i;
    }
  }
  return -1;
}

inline int triangle_next(int vi, const Triangle t) {
  return triangle_l2g((triangle_g2l(vi, t) + 1) % 3, t);
}

inline int triangle_prev(int vi, const Triangle t) {
  return triangle_l2g((triangle_g2l(vi, t) + 2) % 3, t);
}

// Given a global vertex index, return the
// remaining two global vertex indices in
// their original order.
inline Edge triangle_opposite(int vi, const Triangle t) {
  return make_edge(triangle_prev(vi, t), triangle_next(vi, t));
}

inline bool triangle_has_edge(int vi, int vj, const Triangle t) {
  int vii = -1;
  for (int i = 0; vii == -1 && i < 3; ++i) {
    if (t.s[i] == vi) vii = i;
  }
  if (vii == -1) return false;

  int vji = -1;
  for (int i = 0; vji == -1 && i < 3; ++i) {
    if (t.s[i] == vj) vji = i;
  }
  return (vji != -1);
}

#endif
