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

#include "./vector3.h"
#include "./opencl/distance3.h"

#include <assert.h>
#include <limits.h>

#include <list>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <vector>
#include <set>
#include <utility>

#include "./opencl/tribox.h"

#include "./timer.h"
#include "./bb.h"

using namespace std;

namespace oct {

// const int LabeledGeometry3::D = 3;

// Finds the closest distance to an individual geometry
Distance3i distance_geom(const int3& p, const LabeledGeometry3& geometry,
                    int3* verts) {
  Distance3i best = { INT_MAX, int3() };
  const vector<Triangle>& tris = geometry.GetTriangles();
  for (int j = 0; j < tris.size(); ++j) {
    vector<int3> poly_verts;
    set<int3> set_verts;
    const Triangle& t = tris[j];
    for (int k = 0; k < 3; ++k) {
      poly_verts.push_back(verts[t.s[k]]);
      set_verts.insert(verts[t.s[k]]);
    }
    if (poly_verts.size() == set_verts.size()) {
      int3 pverts[] = { poly_verts[0], poly_verts[1], poly_verts[2] };
      best = min_pair3i(best, distance_trianglei(p, pverts));
    } else {
      // Degenerate triangle
      if (set_verts.size() == 1) {
        int3 closest = verts[tris[j].s[0]];
        float3 closest_d =
            make_float3(closest.s[0], closest.s[1], closest.s[2]);
        float3 p_d = make_float3(p.s[0], p.s[1], p.s[2]);
        int dist = static_cast<int>(length(p_d-closest_d)+0.5);
        best = min_pair3i(best, make_dist3(dist, closest));
      } else {
        int3 a = *set_verts.begin();
        int3 b = *set_verts.rbegin();
        Distance3f dist = distance_line3(make_float3(p.s[0], p.s[1], p.s[2]),
                                   make_float3(a.s[0], a.s[1], a.s[2]),
                                   make_float3(b.s[0], b.s[1], b.s[2]));
        Distance3i disti = make_dist3(
            (int)(dist.d+0.5),
            make_int3((int)(dist.p.s[0]),
                     (int)(dist.p.s[1]),
                      (int)(dist.p.s[2])));
        best = min_pair3i(best, disti);
      }
    }
  }
  return best;
}

// Finds the closest distance to a set of geometries
PointAndLabel distance_geoms(
    const int3& p, const vector<LabeledGeometry3>& geometries,
    int3* verts, int* verts_offsets) {
  if (geometries.empty()) {
    throw std::logic_error("zero geometries not supported");
  }
  // int3 t[] = { int3(), int3(), int3() };
  // Distance3i min_dist = make_dist3(INT_MAX, make_int3(0, 0, 0), t);
  Distance3i min_dist = make_dist3(INT_MAX, make_int3(0, 0, 0));
  int idx = -1;
  for (int i = 0; i < geometries.size(); ++i) {
    const Distance3i dist = distance_geom(
        p, geometries[i], verts+verts_offsets[geometries[i].GetLabel()]);

    if (dist.d < min_dist.d) {
      min_dist = dist;
      idx = i;
    }
  }
  PointAndLabel pl(min_dist.p, geometries[idx].GetLabel());
  return pl;
}

// Finds the closest distance to an individual geometry
Distance3i distance_geom(const int3& p, const int* geometry, int3* verts) {
  Distance3i best = { INT_MAX, int3() };
  const int num_tris = geometry[1];
  const Triangle* tris = (Triangle*)(geometry+2);
  for (int j = 0; j < num_tris; ++j) {
    const Triangle& t = tris[j];
    const int3 tri_verts[3] = { verts[t.s[0]], verts[t.s[1]], verts[t.s[2]] };
    int3 set_verts[3];
    const int num_unique = find_unique(tri_verts, set_verts);

    if (num_unique == 3) {
      best = min_pair3i(best, distance_trianglei(p, tri_verts));
    } else {
      // Degenerate triangle
      if (num_unique == 1) {
        // Degenerate to a point
        int3 closest = verts[tris[j].s[0]];
        float3 closest_d =
            make_float3(closest.s[0], closest.s[1], closest.s[2]);
        float3 p_d = make_float3(p.s[0], p.s[1], p.s[2]);
        int dist = (int)(length(p_d - closest_d)+0.5f);
        best = min_pair3i(best, make_dist3(dist, closest));
      } else {
        // Degenerate to a line
        int3 a = set_verts[0];
        int3 b = set_verts[1];
        Distance3f dist = distance_line3(make_float3(p.s[0], p.s[1], p.s[2]),
                                   make_float3(a.s[0], a.s[1], a.s[2]),
                                   make_float3(b.s[0], b.s[1], b.s[2]));
        Distance3i disti = make_dist3(
            (int)(dist.d+0.5),
            make_int3((int)(dist.p.s[0]),
                     (int)(dist.p.s[1]),
                      (int)(dist.p.s[2])));
        best = min_pair3i(best, disti);
      }
    }
  }
  return best;
}

// Finds the closest distance to a set of geometries
PointAndLabel distance_geoms(
    const int3& p, const int* geometries,
    int3* verts, int* verts_offsets) {
  if (geometries[0] == 0) {
    throw std::logic_error("zero geometries not supported");
  }

  // int3 t[] = { int3(), int3(), int3() };
  // Distance3i min_dist = make_dist3(INT_MAX, make_int3(0, 0, 0), t);
  Distance3i min_dist = make_dist3(INT_MAX, make_int3(0, 0, 0));
  int idx = -1;
  // Loop through each geometry
  for (int i = 0; i < geometries[0]; ++i) {
    int offset = geometries[i];
    const Distance3i dist = distance_geom(
        p, geometries+offset, verts+verts_offsets[geometries[offset]]);

    if (dist.d < min_dist.d) {
      min_dist = dist;
      idx = i;
    }
  }
  int idx_offset = geometries[idx];
  PointAndLabel pl(min_dist.p, geometries[idx_offset]);
  return pl;
}

// * are induced points between points a and b
//  _____________________
// |          |          |
// |          |    _a    |
// |          |  _/      |
// |          *_/        |
// |________*/|__________|
// |      _/  |          |0
// |    _/    |          |
// |   b      |          |
// |          |          |
// |__________|__________|
void ClipGeometry(
    const LabeledGeometry3& geometry, const int3& base_point,
    const level_t level,
    vector<vector<LabeledGeometry3> >& clipped,
    int3* vertices) {
  const int3 offsets = base_point;
  const index_t width = Level2CellWidth(level);
  const index_t width2 = width >> 1;  // divides by two
  const index_t width4 = width2 >> 1;  // half the size of child cell
  const vector<Triangle>& triangles = geometry.GetTriangles();

  const size_t t_size = triangles.size();
  vector<Triangle> new_tris[8];

  for (int ti = 0; ti < t_size; ++ti) {
    for (int k = 0; k < 2; ++k) {
      for (int j = 0; j < 2; ++j) {
        for (int i = 0; i < 2; ++i) {
          const int3 center(offsets + make_int3(i, j, k)*width2 + width4);
          if (TriBoxOverlap(center, width4, vertices, triangles[ti])) {
            const int child_idx = 4*k + 2*j + i;
            new_tris[child_idx].push_back(triangles[ti]);
          }
        }
      }
    }
  }

  for (int i = 0; i < 8; ++i) {
    if (!new_tris[i].empty()) {
      clipped[i].push_back(LabeledGeometry3(
          new_tris[i], geometry.GetLabel()));
    }
  }
}
}
