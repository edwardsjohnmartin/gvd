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

#ifndef __VERTEX_H__
#define __VERTEX_H__

#include "./defs.h"
#include "./bit.h"
#include "./vec.h"

NAMESPACE_OCT_BEGIN

//----------------------------------------
// struct Vertex
//----------------------------------------
typedef struct {
  uchar base;
  intn position;
  // level is valid only if base is true
  level_t level;
  int _neighbors[2*DIM];
  // levels (distances) to the neighbors
  level_t neighbor_level[2*DIM];
  int _closest_point;
  int corners[1<<DIM];
} Vertex;

//------------------------------------------------------------
// Function definitions
//------------------------------------------------------------
inline Vertex make_vertex(intn position) {
  Vertex v;
  v.base = false;
  v.position = position;
  v.level = 0;
  for (int i = 0; i < 2*DIM; ++i) {
    v._neighbors[i] = -1;
    v.neighbor_level[i] = 0;
  }
  v._closest_point = -1;
  for (int i = 0; i < (1<<DIM); ++i) {
    v.corners[i] = -1;
  }
  return v;
}

#ifndef OPEN_CL
static bool operator==(Vertex a, Vertex b) {
  if (a.base != b.base) return false;
  if (a.position != b.position) return false;
  if (a.level != b.level) return false;
  for (int i = 0; i < 2*DIM; ++i) {
    if (a._neighbors[i] != b._neighbors[i]) return false;
    if (a.neighbor_level[i] != b.neighbor_level[i]) return false;
  }
  if (a._closest_point != b._closest_point) return false;
  for (int i = 0; i < (1<<DIM); ++i) {
    if (a.corners[i] != b.corners[i]) return false;
  }
  return true;
}
static bool operator!=(Vertex a, Vertex b) {
  return !(a == b);
}
static std::ostream& operator<<(std::ostream& out, const Vertex& v) {
  out << "base = " << v.base << ", position = " << v.position
      << " level = " << (int)v.level << " neighbors = { ";
  for (int i = 0; i < 2*DIM; ++i) {
    out << v._neighbors[i] << "(" << (int)v.neighbor_level[i] << ")";
    if (i < 2*DIM-1)
      out << ", ";
  }
  out << " }, closest point = " << v._closest_point << " corners = { ";
  for (int i = 0; i < (1<<DIM); ++i) {
    out << v.corners[i];
    if (i < (1<<DIM)-1)
      out << ", ";
  }
  out << " }";
  return out;
}
#endif

inline int v_neighbor_index(const Direction d, const Vertex v) {
  return v._neighbors[IndexCardinal(d.pos, d.neg)];
}

inline char v_neighbor_level(const Direction d, const Vertex v) {
  // return v._neighbor_level[IndexCardinal(d.pos, d.neg)] & (~level_mask);
  return v.neighbor_level[IndexCardinal(d.pos, d.neg)];
}

inline void v_set_neighbor(
    const int index, const Direction d, const level_t level,
    __GLOBAL__ Vertex* v) {
  const int nidx = IndexCardinal(d.pos, d.neg);
  v->_neighbors[nidx] = index;
  // v->_neighbor_level[nidx] = (v->_neighbor_level[nidx] & level_mask) + level;
  v->neighbor_level[nidx] = level;
}

inline int v_closest_point(const Vertex v) {
  return v._closest_point;
}

inline void v_set_closest_point(int i, __GLOBAL__ Vertex* v) {
  v->_closest_point = i;
}

inline uchar v_is_base(const Vertex v) {
  // return v._neighbor_level[0] & base_mask;
  return v.base;
}

inline void v_set_is_base(const uchar leaf, __GLOBAL__ Vertex* v) {
  // if (leaf)
  //   v->_neighbor_level[0] |= base_mask;
  // else
  //   v->_neighbor_level[0] &= ~base_mask;
  v->base =  leaf;
}

inline level_t v_cell_level(const Vertex v) {
  assert(v_is_base(v));
  // return
  //     ((v._neighbor_level[1] & level_mask) >> 2) +
  //     ((v._neighbor_level[2] & level_mask) >> 4) +
  //     ((v._neighbor_level[3] & level_mask) >> 6);
  return v.level;
}

inline void v_set_cell_level(level_t level, __GLOBAL__ Vertex* v) {
  assert(v_is_base(*v));
  // v->_neighbor_level[1] = (v->_neighbor_level[1] & (~level_mask)) +
  //     (((level & (3 << 4)) << 2));
  // v->_neighbor_level[2] = (v->_neighbor_level[2] & (~level_mask)) +
  //     (((level & (3 << 2)) << 4));
  // v->_neighbor_level[3] = (v->_neighbor_level[3] & (~level_mask)) +
  //     ((level & 3) << 6);
  v->level = level;
}

NAMESPACE_OCT_END

#endif
