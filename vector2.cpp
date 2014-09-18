#include "./vector2.h"

#include <assert.h>

#include <list>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <vector>
#include <utility>

using namespace std;

namespace oct {

const int LabeledGeometry2::D = 2;

typedef int2 int2;
typedef float2 float2;

struct Distance2f {
  float d;
  float2 p;
};

struct Distance2i {
  int d;
  int2 p;
};

Distance2f make_dist2f(float d, const float2 p) {
  Distance2f dist = { d, p };
  return dist;
}

Distance2i make_dist2i(int d, const int2 p) {
  Distance2i dist = { d, p };
  return dist;
}

Distance2i min_pair(const Distance2i& a, const Distance2i& b) {
  if (a.d < b.d) return a;
  return b;
}

int convert_int_rte(const float f) {
  return (int)(f + 0.5f);
}

int2 convert_int2_rte(const float2 v) {
  return make_int2(convert_int_rte(v.s[0]),
              convert_int_rte(v.s[1]));
}

// assert(oct::distance(
//     float2(1, 0), float2(0, 0), float2(2, 0)).first == 0);
// assert(oct::distance(
//     float2(1, 0), float2(0, 0), float2(2, 0)).second == float2(1, 0));
// assert(oct::distance(
//     float2(1, 1), float2(0, 0), float2(2, 0)).first == 1);
// assert(oct::distance(
//     float2(1, 1), float2(0, 0), float2(2, 0)).second == float2(1, 0));
// assert(oct::distance(
//     float2(1, -1), float2(0, 0), float2(2, 0)).first == 1);
// assert(oct::distance(
//     float2(1, -1), float2(0, 0), float2(2, 0)).second == float2(1, 0));
// assert(oct::distance(
//     float2(-1, 0), float2(0, 0), float2(2, 0)).first == 1);
// assert(oct::distance(
//     float2(-1, 0), float2(0, 0), float2(2, 0)).second == float2(0, 0));
// assert(oct::distance(
//     float2(3, 0), float2(0, 0), float2(2, 0)).first == 1);
// assert(oct::distance(
//     float2(3, 0), float2(0, 0), float2(2, 0)).second == float2(2, 0));
// assert(oct::distance(
//     float2(-1, 1), float2(0, 0), float2(2, 0)).first == sqrt(2));
// assert(oct::distance(
//     float2(-1, 1), float2(0, 0), float2(2, 0)).second == float2(0, 0));
// pair<double, float2> distance_line2f(
Distance2f distance_line2f(
    const float2& p, const float2& a, const float2& b) {
  // const double dot = (p-a)*((b-a).unit());
  const double dotf = dot((p-a), (normalize(b-a)));
  if (dotf < 0) {
    // return make_dist2f((p-a).norm(), a);
    return make_dist2f(length(p-a), a);
  // } else if (dot > (b-a).norm()) {
  } else if (dotf > length(b-a)) {
    // return make_dist2f((p-b).norm(), b);
    return make_dist2f(length(p-b), b);
  } else {
    double3 pa = make_double3(convert_double2(p-a));
    double3 ba = normalize(make_double3(convert_double2(b-a)));
    const double dist = length(cross(pa,ba));
    const float2 pi = a + normalize(b-a)*dotf;
    return make_dist2f(dist, pi);
  }
}

Distance2i distance_line2i(const int2& p, const int2& a, const int2& b) {
  Distance2f d = distance_line2f(make_float2(p.s[0], p.s[1]),
                                 make_float2(a.s[0], a.s[1]),
                                 make_float2(b.s[0], b.s[1]));
  return make_dist2i(convert_int_rte(d.d),
                     convert_int2_rte(d.p));
}

// Finds the distance from p to q where q is the closest point
// on polyline.  Returns <|q-p|, q>.
Distance2i distance_geom(const int2& p, const LabeledGeometry2& geometry) {
  Distance2i best = { INT_MAX, int2() };
  const int2* vertices = geometry.GetVertices();
  for (int i = 0; i < geometry.size(); ++i) {
    const Edge& e = geometry[i];
    best = min_pair(best, distance_line2i(
        p, vertices[e.s[0]], vertices[e.s[1]]));
  }
  return best;
}

Distance2i distance_geom(const int2& p, const int* geometry, int2* verts) {
  Distance2i best = { INT_MAX, int2() };
  // const int2* vertices = geometry.GetVertices();
  const int2* vertices = verts;//geometry.GetVertices();
  const int num_edges = geometry[1];
  const Edge* edges = (Edge*)(geometry+2);
  // for (int i = 0; i < geometry.size(); ++i) {
  for (int i = 0; i < num_edges; ++i) {
    // const Edge& e = geometry[i];
    const Edge& e = edges[i];
    best = min_pair(best, distance_line2i(
        p, vertices[e.s[0]], vertices[e.s[1]]));
  }
  return best;
}

PointAndLabel distance_geoms(
    const int2& p, const vector<LabeledGeometry2>& geometries,
    int2* verts, int* verts_offsets) {
  if (geometries.empty()) {
    throw std::logic_error("zero geometries not supported");
  }

  Distance2i min_dist = { INT_MAX, int2() };
  int idx = -1;
  for (int i = 0; i < geometries.size(); ++i) {
    const Distance2i dist = distance_geom(p, geometries[i]);
    if (dist.d < min_dist.d) {
      min_dist = dist;
      idx = i;
    }
  }
  if (min_dist.d == INT_MAX) {
    cerr << "No closest point" << endl;
    throw logic_error("No closest point");
  }
  return PointAndLabel(min_dist.p, geometries[idx].GetLabel());
}

PointAndLabel distance_geoms(
    const int2& p, const int* geometries,
    int2* verts, int* verts_offsets) {
  // if (geometries.empty()) {
  if (geometries[0] == 0) {
    throw std::logic_error("zero geometries not supported");
  }

  Distance2i min_dist = { INT_MAX, int2() };
  int idx = -1;
  // for (int i = 0; i < geometries.size(); ++i) {
  for (int i = 0; i < geometries[0]; ++i) {
    int offset = geometries[i];
    // const Distance2i dist = distance_geom(p, geometries[i]);
    const Distance2i dist = distance_geom(
        p, geometries+offset, verts+verts_offsets[geometries[offset]]);
    if (dist.d < min_dist.d) {
      min_dist = dist;
      idx = i;
    }
  }
  int idx_offset = geometries[idx];
  return PointAndLabel(min_dist.p, geometries[idx_offset]);
}

// dim is the dimension.
//   dim = 0 for split into left/right
//   dim = 1 for split into bottom/top
void Split(vector<int2>& vertices,
           const vector<Edge>& edges,
           const int2& mid,
           vector<Edge>& lower, vector<Edge>& upper,
           const int dim) {
  for (int i = 0; i < edges.size(); ++i) {
    const Edge& e = edges[i];
    int pi = e.s[0];
    int qi = e.s[1];
    if (vertices[qi].s[dim] < vertices[pi].s[dim]) swap(pi, qi);
    int2 p = vertices[pi];
    int2 q = vertices[qi];
    if (q.s[dim] < mid.s[dim]) {
      lower.push_back(e);
    } else if (p.s[dim] > mid.s[dim]) {
      upper.push_back(e);
    } else {
      // split
      const double m =
          (q.s[1-dim]-p.s[1-dim]) / static_cast<double>((q.s[dim]-p.s[dim]));
      int2 r = make_int2(0);
      r.s[dim] = mid.s[dim];
      r.s[1-dim] = p.s[1-dim] + static_cast<int>((mid.s[dim] - p.s[dim]) * m + 0.5);

      const int ri = vertices.size();
      vertices.push_back(r);
      lower.push_back(make_edge(pi, ri));
      upper.push_back(make_edge(ri, qi));
    }
  }
}

void ClipGeometry(
    const LabeledGeometry2& geometry, const int2& base_point,
    const level_t level,
    vector<vector<LabeledGeometry2> >& clipped,
    int2* vertices_) {
  const index_t width = CellWidth(level);
  const index_t width2 = width >> 1;  // divides by two
  const int2 base = base_point;
  const int2 mid = base_point + width2;

  // shared_ptr<vector<int2> > vertices = geometry.GetVertices();
  int2* vertices = geometry.GetVertices();
  const vector<Edge>& edges = geometry.GetEdges();
  const vector<SATData2>& axes = geometry.GetAxes();

  // The current geometry split into 4
  LabeledGeometry2 split[4];
  for (int i = 0; i < 4; ++i) {
    split[i] = LabeledGeometry2(vertices, geometry.GetLabel());
    // split[i] = LabeledGeometry2(geometry.GetLabel());
    for (int j = 0; j < edges.size(); ++j) {
      const int2 v = base + make_int2(width2*(i&1), width2*((i&2)>>1));
      if (Intersects(v, width2, axes[j])) {
        split[i].Add(edges[j], axes[j]);
      }
    }
  }
  for (int i = 0; i < 4; ++i) {
    if (!split[i].empty()) {
      clipped[i].push_back(split[i]);
    }
  }
}

// 1010 |  1110 |  1011 |  1111
//      |       |       |
// -----------------------------
// 0010 |  0110 |  0011 |  0111
//      |       |       |
// -----------------------------
// 1000 |  1100 |  1001 |  1101
//      |       |       |
// -----------------------------
// 0000 |  0100 |  0001 |  0101
//      |       |       |
void outcode16(const int2&a, const int2& base, const int width) {
  // int code = 0;
  // const int w2 = width >> 1;
  // const int w4 = width >> 2;
  // if (a[0] > base[0] + w2) {
  //   code += 1;
  //   if (a[0] > base[0] + w2 + w4) {
  //     code += 4;
  //   }
  // } else {
  //   if (a[0] > base[0] + w4) {
  //     code += 4;
  //   }
  // }
  // if (a[1] > base[1] + w2) {
  //   code += 2;
  //   if (a[1] > base[1] + w2 + w4) {
  //     code += 8;
  //   }
  // } else {
  //   if (a[1] > base[1] + w4) {
  //     code += 8;
  //   }
  // }
  // return code;
}

// 1010 |  1110 |  1011 |  1111
//      |       |       |
// -----------------------------
// 0010 |  0110 |  0011 |  0111
//      |       |       |
// -----------------------------
// 1000 |  1100 |  1001 |  1101
//      |       |       |
// -----------------------------
// 0000 |  0100 |  0001 |  0101
//      |       |       |
void Intersects(int2 a, int2 b,
                const int2& base, const int width,
                bool intersects[]) {
  // if (a[0] > b[0]) swap(a, b);
  // const double m = (b[1]-a[1]) / static_cast<double>(b[0]-a[0]);
  // const int acode = outcode16(a, base, width);
  // const int bcode = outcode16(b, base, width);
  // if (acode & bcode & 1) {
  //   // right-hand side
  //   acode = EraseBit(acode, 0);
  //   bcode = EraseBit(bcode, 0);
  //   if (acode & bcode & 1) {
  //     // top
  //     intersects[3] = true;
  //   } else if (~acode & ~bcode & 1) {
  //     // bottom
  //     intersects[1] = true;
  //   } else {
  //     // compute intersection with y = base[1] + w2
  //     const int x = (base[1]-a[1]) * m;
  //     if (x
  //   }
  // }
}

// bool Test(int base_x, int base_y, int width,
//           const int2& v0, const int2& v1, const int2& v2) {
//   int2 base(base_x, base_y);
//   int2 verts[] = { v0, v1, v2 };
//   vector<int2> vertices(verts, verts+3);
//   return oct::Intersects(base, width, Triangle(0, 1, 2), vertices);
// }

// bool Test(int base_x, int base_y, int width,
//           const int2& v0, const int2& v1) {
//   int2 base(base_x, base_y);
//   int2 verts[] = { v0, v1 };
//   vector<int2> vertices(verts, verts+2);
//   return oct::Intersects(base, width, Edge(0, 1), vertices);
// }

// int main(int argc, char** argv) {
//   assert(Test(0, 0, 10, int2(0, 0), int2(10, 0), int2(10, 10)) == true);
//   assert(Test(0, 0, 10, int2(5, 5), int2(10, 5), int2(10, 10)) == true);
//   assert(Test(0, 0, 5, int2(5, 5), int2(10, 5), int2(10, 10)) == true);
//   assert(Test(0, 0, 2, int2(5, 5), int2(10, 5), int2(10, 10)) == false);
//   assert(Test(0, 0, 10, int2(11, 0), int2(15, 0), int2(11, 10)) == false);
//   assert(Test(0, 0, 10, int2(11, 0), int2(15, 0), int2(-1, 15)) == true);
//   assert(Test(0, 0, 10, int2(11, 0), int2(12, 0), int2(0, 15)) == true);
//   assert(Test(0, 0, 10, int2(21, -1), int2(25, -1), int2(-1, 25)) ==
//          false);
//   assert(Test(0, 0, 10, int2(-1, 21), int2(-1, 25), int2(25, -1)) ==
//          false);
//   assert(Test(0, 0, 10, int2(-1, 25), int2(-1, 21), int2(25, -1)) ==
//          false);
//   assert(Test(0, 0, 10, int2(5, 5), int2(5, 6), int2(6, 6)) == true);
//   assert(Test(0, 0, 10, int2(5, 5), int2(15, 5), int2(5, 6)) == true);
//   assert(Test(0, 0, 10, int2(-1, -1), int2(100, -1), int2(-1, 100)) ==
//          true);

//   assert(Test(0, 0, 10, int2(-1, -1), int2(10, -1)) == false);
//   assert(Test(0, 0, 10, int2(-1, -1), int2(10, 5)) == true);
//   return 0;


//                            a6
//     ------------  a3       |\
//    |            |          | \_
//    |            |          |   \_
//    |            |          |     \_
//    |            |          |       \_
//     ------------  a2       |_________\
//   a0            a1         a4         a5
//
// Uses the Separating Axis Theorem (SAT)
// bool Intersects(const int2& base, const int width,
//                 const Triangle& tri, const std::vector<int2>& vertices) {
//   bool left[7];
//   bool right[7];
//   std::fill(left, left+7, false);
//   std::fill(right, right+7, false);
//   vector<int2> t;
//   for (int i = 0; i < 3; ++i) {
//     t.push_back(vertices[tri[i]]);
//   }
//   for (int i = 0; i < 3; ++i) {
//     const int2& v = t[i];
//     for (int axis = 0; axis < 2; ++axis) {
//       // a0/a2
//       if (v[axis] < base[axis])
//         left[axis<<1] = true;
//       else
//         right[axis<<1] = true;
//       // a1/a3
//       if (v[axis] > base[axis] + width)
//         right[(axis<<1)+1] = true;
//       else
//         left[(axis<<1)+1] = true;
//     }
//   }

//   // Get the 3 axis vectors for the triangle
//   float2 axes[3];
//   float norms[3];
//   Range ranges[3];
//   float2 axes_n[3];
//   float dots[3][4];  // dot products
//   for (int i = 0; i < 3; ++i) {
//     axes[i] = Ortho(ToFloat(t[(i+1)%3] - t[i]));
//     norms[i] = axes[i].norm();
//     axes_n[i] = axes[i] / norms[i];
//     ranges[i] = Range();
//     for (int j = 0; j < 3; ++j) {
//       ranges[i](ToFloat(t[j]) * axes_n[i]);
//     }
//     for (int j = 0; j < 4; ++j) {
//       const int2 v = base + int2(width*(j&1), width*((j&2)>>1));
//       dots[i][j] = ToFloat(v) * axes_n[i];
//     }
//   }
  
//   if (!right[0] || !left[1] || !right[2] || !left[3])
//     return false;

//   // Now for each of the 3 triangle axes
//   for (int i = 0; i < 3; ++i) {
//     Range r;
//     for (int j = 0; j < 4; ++j) {
//       r(dots[i][j]);
//     }
//     if (!ranges[i].Intersects(r))
//       return false;
//   }
//   return true;
// }

//                            a5
//     ------------  a3        \
//    |            |            \_
//    |            |              \_
//    |            |                \_
//    |            |                  \_
//     ------------  a2                 \
//   a0            a1                    a4
//
// Uses the Separating Axis Theorem (SAT)
bool Intersects(const int2& base, const int width,
                const SATData2& d) {
  // This strange way of doing it is do serve as a placeholder
  // in case we want to add optimizations when computing
  // intersections with sub-cells.
  bool left[7];
  bool right[7];
  std::fill(left, left+7, false);
  std::fill(right, right+7, false);
  for (int i = 0; i < 2; ++i) {
    const int2& v = d[i];
    for (int axis = 0; axis < 2; ++axis) {
      // a0/a2
      if (v.s[axis] < base.s[axis])
        left[axis<<1] = true;
      else
        right[axis<<1] = true;
      // a1/a3
      if (v.s[axis] > base.s[axis] + width)
        right[(axis<<1)+1] = true;
      else
        left[(axis<<1)+1] = true;
    }
  }

  if (!right[0] || !left[1] || !right[2] || !left[3])
    return false;

  // dot products with cube
  Range r;
  for (int j = 0; j < 4; ++j) {
    const int2 v = base + make_int2(width*(j&1), width*((j&2)>>1));
    r(dot(SATData2::ToFloat(v), d.plane));
  }

  return d.range.Intersects(r);
}
}
