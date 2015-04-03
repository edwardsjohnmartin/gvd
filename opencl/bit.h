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

#ifndef __BIT_H__
#define __BIT_H__

#include "./defs.h"

NAMESPACE_OCT_BEGIN

// Usage:
//   kLviStrides[axis]
// extern __CONST__ int kLviStrides[3];
// static __CONST__ int kLviStrides[3] = { 1, 2, 4 };

// Usage:
//   kSlviStrides[axis]
// extern __CONST__ int kSlviStrides[3];
static __CONST__ int kSlviStrides[3] = { 1, 3, 9 };

// Stores the local index of the center vertex in a subdivided cell
// Usage:
//   kSlviCenter[dim]
// extern __CONST__ int kSlviCenter[4];
// static __CONST__ int kSlviCenter[4] = { 0, 1, 4, 13 };

// Local vertex index (lvi)
//
//        *-------------------*
//        |2                  |3
//        |                   |
//        |                   |
//        |                   |
//        |                   |
//        |                   |
//        |                   |
//        *-------------------*
//         0                   1
//       
//             
//             *------------------*
//           / |6                /|7
//          /                   / |
//         /   |               /  |
//        *-------------------*   |
//        |2   |              |3  |
//        |                   |   |
//        |    |              |   |
//        |    *- - - - - - - |- -*
//        |  /  4             |  / 5
//        | /                 | /
//        |/                  |/
//        *-------------------*
//         0                   1
//  y ^   
//    | _/ z
//    |/
//    --------> x

// Subdivided local vertex index (slvi)
//
//        *---------*---------*
//        |6        |7        |8
//        |         |         |
//        |         |         |
//        *---------*---------* 
//        |3        |4        |5
//        |         |         |
//        |         |         |
//        *---------*---------*
//         0         1         2
//       
//             
//             *---------*--------*
//           /  24     /  25     /|26
//          *---------*---------* |
//         / 15      / 16      /| |
//        *---------*---------* | *
//        |6        |7        |8|/|23
//        |         |         | * |
//        |         |         |/| |
//        *---------*---------* | *
//        |3        |4        |5|/ 20
//        |         |         | *
//        |         |         |/ 11
//        *---------*---------*
//         0         1         2
//  y ^   
//    | _/ z
//    |/
//    --------> x

// Edge index
//
//          __________7___________
//         /|                    /|
//        /                     / |
//      1/  |                 3/  |
//      /                     /   |
//     /  10|                /  11|
//    /___________6_________/     |
//   |      |               |     |
//   |                      |     |
//   |      |               |     |
//   |       _  _  _ 5_  _  |_  _ |
//   |     /                |     /
//  8|   0                 9|   2/
//   |   /                  |   /
//   |                      |  /
//   | /                    | /
//   |__________4___________|/
//

// extern __CONST__ int lvi2slvi[8];
static __CONST__ int lvi2slvi[8] = { 0, 2, 6, 8, 18, 20, 24, 26 };
// static __CONST__ int slvi2lvi[8] = {
//    0, -1,  1, -1, -1, -1,  2, -1,  3,
//   -1, -1, -1, -1, -1, -1, -1, -1, -1,
//    4, -1,  5, -1, -1, -1,  6, -1,  7 };

// Returns the local index of the neighbor vertex in the positive direction
// in the given axis.  Returns -1 if no such vertex exists (such as is the
// case when moving in the positive direction along the axis moves outside
// the cell boundary).
// lvi = [0, 1, ..., 2^D]
// axis = [0, 1, ..., D-1]
// template <int D>
inline int NbrVertexPositive(const int lvi, const int axis) {
  if (lvi & (1 << axis)) return -1;
  const int index = lvi + (1 << axis);
  return index;
}

// Finds the index in a cardinal direction.  In other words, no diagonals
// are allowed.  The number of such indices is 2 * D.
//
//         _ negative (0) or positive (1) bit
//        | |
//   0 0 1 1
//  |_____|
//  dimension bits
//
// In the example above, 3 = neighbor in positive first dimension.
// Other examples:
//     0 0 0 0 - negative 0th dimension
//     0 1 0 1 - positive 2nd dimension
// template <int D>
inline int IndexCardinal(const int pos, const int neg) {
  // no index exists for no direction
  assert(pos || neg);

  int index = -1;
  for (int i = 0; i < DIM; ++i) {
    const int test = 0x1 << i;
    if (pos & test) {
      assert(index == -1);
      index = (i << 1) | 1;
    } else if (neg & test) {
      assert(index == -1);
      index = (i << 1) | 0;
    }
  }
  return index;
}

// Finds the index of an edge.  Let d1 be primary axis and d2 be secondary
// axis.  d2 = (d1+1)%D.
//
//        |_negative (0) or positive (1) in primary dimension
//        | | _negative (0) or positive (1) in secondary dimension
//        | | |
//   0 0 1 1 1
//  |_____|
//  primary dimension bits
//
// In the example above, 3 = neighbor in positive first dimension.
//  0: 0 0|0|0
//  1: 0 0|0|1
//  2: 0 0|1|0
//  3: 0 0|1|1
//  4: 0 1|0|0
//  5: 0 1|0|1
//  6: 0 1|1|0
//  7: 0 1|1|1
//  8: 1 0|0|0
//  9: 1 0|0|1
// 10: 1 0|1|0
// 11: 1 0|1|1
// Tests:
// assert(IndexEdge(0, 0, 0) == 0);
// assert(IndexEdge(0, 0, 1) == 1);
// assert(IndexEdge(0, 1, 0) == 2);
// assert(IndexEdge(0, 1, 1) == 3);
// assert(IndexEdge(1, 0, 0) == 4);
// assert(IndexEdge(1, 0, 1) == 5);
// assert(IndexEdge(1, 1, 0) == 6);
// assert(IndexEdge(1, 1, 1) == 7);
// assert(IndexEdge(2, 0, 0) == 8);
// assert(IndexEdge(2, 0, 1) == 9);
// assert(IndexEdge(2, 1, 0) == 10);
// assert(IndexEdge(2, 1, 1) == 11);
inline int IndexEdge(const int primary_axis,
              const bool primary_pos, const bool secondary_pos) {
  return (primary_axis<<2) + (primary_pos?2:0) + (secondary_pos?1:0);
}

// To understand the direction of the neighbor, consider the node in
// quadrant 0.  The neighbor in the "1 positive direction" will be the node
// in quadrant 1.  The neighbor in the "3 positive direction" will be the
// node in quadrant 3.  Now consider the node in quadrant 3.  The neighbor
// in the "1 negative direction" will be the node in quadrant 2.  And so on.
// The following table gives further illustration.
//
//   pos_dir | neg_dir | neighbor
//  ------------------------------
//      0    |    0    | itself
//      1    |    0    | right
//      0    |    1    | left
//      2    |    0    | above
//      0    |    2    | below
//      3    |    0    | right/above
//      0    |    3    | left/below
//      2    |    1    | left/above
//      1    |    2    | right/below
typedef struct {
  int pos, neg;
} Direction;

inline Direction DirectionFromPosNeg(int pos, int neg) {
  const Direction d = { pos, neg };
  return d;
}

inline Direction DirectionFromPosAxis(const int axis) {
  const Direction d = { (1 << axis), 0 };
  return d;
}

inline Direction DirectionFromAxis(const int axis, const bool pos) {
  const int i = (1 << axis);
  if (pos) return DirectionFromPosNeg(i, 0);
  return DirectionFromPosNeg(0, i);
}

inline Direction DirectionFromLvi(const int lvi) {
  const int pos = lvi;
  const int neg = lvi ^ ((1<<DIM)-1);
  return DirectionFromPosNeg(pos, neg);
}

inline Direction DirectionFromSlvi(int slvi) {
  int pos = 0, neg = 0;
  for (int i = 0; i < DIM; ++i) {
    const int mask = (1<<i);
    if (slvi % 3 == 0) {
      neg += mask;
    } else if (slvi % 3 == 1) {
      // do nothing
    } else {
      pos += mask;
    }
    slvi /= 3;
  }
  return DirectionFromPosNeg(pos, neg);
}

inline Direction Reversed(const Direction* d) {
  return DirectionFromPosNeg(d->neg, d->pos);
}

// <0, 0>.Add(1, 0) = <1, 0>
// <1, 0>.Add(0, 1) = <0, 0>
// Tests:
// assert(oct::Direction<3>(0, 0).Add(0, 0) == oct::Direction<3>(0, 0));
// assert(oct::Direction<3>(0, 0).Add(1, 0) == oct::Direction<3>(1, 0));
// assert(oct::Direction<3>(0, 0).Add(2, 1) == oct::Direction<3>(2, 1));
// assert(oct::Direction<3>(0, 0).Add(4, 2) == oct::Direction<3>(4, 2));
// assert(oct::Direction<3>(0, 7).Add(0, 0) == oct::Direction<3>(0, 7));
// assert(oct::Direction<3>(0, 7).Add(1, 0) == oct::Direction<3>(0, 6));
// assert(oct::Direction<3>(0, 7).Add(2, 0) == oct::Direction<3>(0, 5));
// assert(oct::Direction<3>(0, 7).Add(3, 0) == oct::Direction<3>(0, 4));
// assert(oct::Direction<3>(4, 2).Add(0, 0) == oct::Direction<3>(4, 2));
// assert(oct::Direction<3>(4, 2).Add(1, 0) == oct::Direction<3>(5, 2));
// assert(oct::Direction<3>(4, 2).Add(0, 1) == oct::Direction<3>(4, 3));
// assert(oct::Direction<3>(4, 2).Add(0, 5) == oct::Direction<3>(0, 3));
inline Direction Add(const Direction* d, int pos, int neg) {
  int new_pos = d->pos;
  int new_neg = d->neg;
  for (int i = 0; i < DIM; ++i) {
    const int mask = (1<<i);
    if (pos & mask) {
      if (new_neg & mask)
        new_neg = new_neg - mask;
      else
        new_pos += mask;
    }
    if (neg & mask) {
      if (new_pos & mask)
        new_pos = new_pos - mask;
      else
        new_neg += mask;
    }
  }
  return DirectionFromPosNeg(new_pos, new_neg);
}

// Finds the index in a cardinal direction.  In other words, no diagonals
// are allowed.  The number of such indices is 2 * D.
//
//         _ negative (0) or positive (1) bit
//        | |
//   0 0 1 1
//  |_____|
//  dimension bits
//
// In the example above, 3 = neighbor in positive first dimension.
// Other examples:
//     0 0 0 0 - negative 0th dimension
//     0 1 0 1 - positive 2nd dimension
// int IndexCardinal(const Direction* d) {
//   return oct::IndexCardinal<DIM>(d->pos, d->neg);
// }

// Finds the index in any direction, cardinal or diagonal.  There are
// 3^D directions.
//
// Tests:
// assert(oct::Direction<2>(0, 3).Index() == 0);
// assert(oct::Direction<2>(0, 2).Index() == 1);
// assert(oct::Direction<2>(1, 2).Index() == 2);
// assert(oct::Direction<2>(0, 1).Index() == 3);
// assert(oct::Direction<2>(0, 0).Index() == 4);
// assert(oct::Direction<2>(1, 0).Index() == 5);
// assert(oct::Direction<2>(2, 1).Index() == 6);
// assert(oct::Direction<2>(2, 0).Index() == 7);
// assert(oct::Direction<2>(3, 0).Index() == 8);
// assert(oct::Direction<3>(0, 7).Index() == 0);
// assert(oct::Direction<3>(0, 6).Index() == 1);
// assert(oct::Direction<3>(1, 6).Index() == 2);
// assert(oct::Direction<3>(0, 5).Index() == 3);
// assert(oct::Direction<3>(2, 5).Index() == 6);
// assert(oct::Direction<3>(0, 3).Index() == 9);
// assert(oct::Direction<3>(0, 0).Index() == 13);
// assert(oct::Direction<3>(4, 3).Index() == 18);
// assert(oct::Direction<3>(7, 0).Index() == 26);
inline int Slvi(const Direction* d) {
  int index = 0;
  int f = 1;
  for (int i = 0; i < DIM; ++i) {
    const int mask = 0x1 << i;
    if (d->neg & mask) {
      // add nothing
    } else if (d->pos & mask) {
      index += (f<<1);
    } else {
      index += f;
    }
    f *= 3;
  }
  return index;
}

NAMESPACE_OCT_END

#endif
