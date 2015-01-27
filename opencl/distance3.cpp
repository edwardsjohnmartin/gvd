#include <limits.h>

#include "./distance3.h"

#include "./geometry.h"

#ifndef OPEN_CL
#include "./vec.h"
#endif

NAMESPACE_OCT_BEGIN

/***************************************************
 * Utility functions
 ***************************************************/
Distance3i make_dist3(int d, const int3 p) {
  Distance3i dist = { d, p };
  return dist;
}

Distance3f min_pair3f(const Distance3f a, const Distance3f b) {
  if (a.d < b.d) return a;
  return b;
}

Distance3i min_pair3i(const Distance3i a, const Distance3i b) {
  if (a.d < b.d) return a;
  return b;
}

/***************************************************
 * Distance functions
 ***************************************************/

Distance3f distance_line3(const float3 p, const float3 a, const float3 b) {
  const float3 bma = b-a;
  const float3 bma_u = fast_normalize(bma);
  const float ab = fast_length(bma);
  const float3 pma = p - a;
  const float dotf = dot(pma, bma_u);
  if (dotf < 0) {
    Distance3f d = { fast_length(pma), a };
    return d;
  } else if (dotf > ab) {
    Distance3f d = { fast_length(pma), b };
    return d;
  } else {
    const float dist = fast_length(cross(pma, bma_u));
    const float3 pi = a + (bma_u * dotf);
    Distance3f d = { dist, pi };
    return d;
  }
}

#ifdef OPEN_CL
// Squared length
inline float length2(const float3 v) {
  return dot(v, v);
}
#endif

Distance3f distance_trianglef(const float3 p, const float3* poly) {
  const float3 normal = normalize(cross(poly[1]-poly[0], poly[2]-poly[0]));

  // Degenerate triangle
  if (normal.x != normal.x) {
    // Hack to safely fail if triangle is degenerate
    float max_d = length2(poly[1]-poly[0]);
    float3 a = poly[1];
    float3 b = poly[0];
    if (length2(poly[2]-poly[1]) > max_d) {
      max_d = length2(poly[2]-poly[1]);
      a = poly[2];
      b = poly[1];
    }
    if (length2(poly[2]-poly[0]) > max_d) {
      // No need to reset max_d
      a = poly[2];
      b = poly[0];
    }
    return distance_line3(p, a, b);
  }

  const float proj_dist = dot((p - poly[0]), normal);

  const float3 proj_pnt = p - (normal * proj_dist);
  float dist = fabs(proj_dist);

  int drop_dim = 0;
  float drop_dim_val = fabs(normal.x);
  if (fabs(normal.y) > drop_dim_val) {
    drop_dim = 1;
    drop_dim_val = fabs(normal.y);
  }
  if (fabs(normal.z) > drop_dim_val) {
    drop_dim = 2;
    drop_dim_val = fabs(normal.z);
  }

  float2 poly_proj[3];
  float2 proj_proj_pnt;

#ifdef OPEN_CL
  for (int i = 0; i < 3; ++i) {
    if (drop_dim == 0) {
      poly_proj[i] = poly[i].yz;
    } else if (drop_dim == 1) {
      poly_proj[i] = poly[i].xz;
    } else {
      poly_proj[i] = poly[i].xy;
    }
  }
  if (drop_dim == 0) {
    proj_proj_pnt = proj_pnt.yz;
  } else if (drop_dim == 1) {
    proj_proj_pnt = proj_pnt.xz;
  } else {
    proj_proj_pnt = proj_pnt.xy;
  }
#else
  for (int i = 0; i < 3; ++i) {
    if (drop_dim == 0) {
      poly_proj[i] = make_float2(poly[i].y, poly[i].z);
    } else if (drop_dim == 1) {
      poly_proj[i] = make_float2(poly[i].x, poly[i].z);
    } else {
      poly_proj[i] = make_float2(poly[i].x, poly[i].y);
    }
  }
  if (drop_dim == 0) {
    proj_proj_pnt = make_float2(proj_pnt.y, proj_pnt.z);
  } else if (drop_dim == 1) {
    proj_proj_pnt = make_float2(proj_pnt.x, proj_pnt.z);
  } else {
    proj_proj_pnt = make_float2(proj_pnt.x, proj_pnt.y);
  }
#endif

  float2 poly_shift[3];
  for (int i = 0; i < 3; ++i) {
    poly_shift[i] = poly_proj[i] - proj_proj_pnt;
  }

  bool test_val;
  bool first_time = true;
  bool inside = true;
  for (int i = 0; i < 3; ++i) {
    float2 v = poly_shift[i];
    float2 vn = poly_shift[(i+1) % 3];
    float area = vn.x*v.y - v.x*vn.y;  // actually is twice area
    if (first_time) {
      test_val = area > 0;
      first_time = false;
    } else {
      if (test_val != area > 0) {
        inside = false;
        break;
      }
    }
  }

  if (inside) {
    // nan!
    /* assert(dist == dist); */

    // dist is set at start of function to be proj distance
    Distance3f d = { dist, proj_pnt };
    return d;
  } else {
    bool unset = true;
    Distance3f best;
    for (int i = 0; i < 3; ++i) {
      if (unset) {
        best = distance_line3(p, poly[i], poly[(i+1)%3]);
        unset = false;
      } else {
        Distance3f dist = distance_line3(p, poly[i], poly[(i+1)%3]);
        best = min_pair3f(best, dist);
      }
    }
    // nan!
    /* assert(best.d == best.d) { */
    return best;
  }
}

// function used to go from int distances to float and then back
// when finding distance to an individual polygon
Distance3i distance_trianglei(const int3 p, const int3* poly) {
  float3 poly_float3[3];
  for (int i = 0; i < 3; ++i) {
    /* poly_float3[i] = (float3)(poly[i].x, poly[i].y, poly[i].z); */
    poly_float3[i] = convert_float3(poly[i]);
  }
  Distance3f d =
      /* distance_trianglef((float3)(p.x, p.y, p.z), poly_float3); */
      distance_trianglef(convert_float3(p), poly_float3);
  return make_dist3(
      // convert_int_rte(d.d),
      // convert_int3_rte(d.p));
      convert_int(d.d),
      convert_int3(d.p));
}

int find_unique(const int3* verts, int3* unique_verts) {
  int num_unique = 0;
  for (int k = 0; k < 3; ++k) {
    int is_unique = 1;
    for (int l = 0; is_unique && l < num_unique; ++l) {
      is_unique = any(verts[k] != unique_verts[l]);
    }
    if (is_unique) {
      unique_verts[num_unique++] = verts[k];
    }
  }
  return num_unique;
}

Distance3i distance_geom3(
    const int3 p, global const int* geometry, global const int3* verts) {
  // Distance3i best = { INT_MAX, (int3)(0, 0, 0) };
  Distance3i best = { INT_MAX, make_int3(0, 0, 0) };
  const int num_tris = geometry[1];
  const Triangle* tris = (Triangle*)(geometry+2);
  for (int j = 0; j < num_tris; ++j) {
    Triangle t = tris[j];
    const int3 tri_verts[3] =
        { (verts[t.s[0]]),
          (verts[t.s[1]]),
          (verts[t.s[2]]) };
    int3 set_verts[3];
    const int num_unique = find_unique(tri_verts, set_verts);

    if (num_unique == 3) {
      best = min_pair3i(best, distance_trianglei(p, tri_verts));
    } else {
      // Degenerate triangle
      if (num_unique == 1) {
        // Degenerate to a point
        int3 closest = tri_verts[0];
        float3 closest_d = convert_float3(closest);
        float3 p_d = convert_float3(p);
        int dist = (int)(fast_length(p_d-closest_d)+0.5f);
        best = min_pair3i(best, make_dist3(dist, closest));
      } else {
        // Degenerate to a line
        int3 a = set_verts[0];
        int3 b = set_verts[1];
        Distance3f dist = distance_line3(convert_float3(p),
                                        convert_float3(a),
                                        convert_float3(b));
        Distance3i disti = make_dist3(
            // convert_int_rte(dist.d),
            // convert_int3_rte(dist.p));
            convert_int(dist.d),
            convert_int3(dist.p));
        best = min_pair3i(best, disti);
      }
    }
  }
  return best;
}

/* Finds the closest distance to a set of geometries */
PointAndLabel3 distance_geoms3(
    const int3 p_, global const int* geometries,
    global const int3* verts, global const int* verts_offsets) {
  const int3 p = p_;
  Distance3i min_dist = make_dist3(INT_MAX, make_int3(0, 0, 0));
  int idx = -1;
  for (int i = 0; i < geometries[0]; ++i) {
    int offset = geometries[i];
    const Distance3i dist = distance_geom3(
        p, geometries+offset, verts+verts_offsets[geometries[offset]]);

    if (dist.d < min_dist.d) {
      min_dist = dist;
      idx = i;
    }
  }
  // idx, idx_offset are correct.  Must be the actual point from
  // a higher-level dist function
  int offset = geometries[idx];
  const int label = geometries[offset];
  PointAndLabel3 pl = { min_dist.p, label };
  return pl;
}

// Computes the distance of a point from a geometry
Distance3i distance_geometry(
    const int3 p, Geometry geometry, __GLOBAL__ const int3* verts) {
  Distance3i best = { INT_MAX, make_int3(0, 0, 0) };
  const int num_tris = g_m(geometry);
  __GLOBAL__ const Triangle* tris = geometry.faces;
  for (int j = 0; j < num_tris; ++j) {
    Triangle t = tris[j];
    const int3 tri_verts[3] =
        { (verts[t.s[0]]),
          (verts[t.s[1]]),
          (verts[t.s[2]]) };
    int3 set_verts[3];
    const int num_unique = find_unique(tri_verts, set_verts);

    if (num_unique == 3) {
      best = min_pair3i(best, distance_trianglei(p, tri_verts));
    } else {
      // Degenerate triangle
      if (num_unique == 1) {
        // Degenerate to a point
        int3 closest = tri_verts[0];
        float3 closest_d = convert_float3(closest);
        float3 p_d = convert_float3(p);
        int dist = (int)(fast_length(p_d-closest_d)+0.5f);
        best = min_pair3i(best, make_dist3(dist, closest));
      } else {
        // Degenerate to a line
        int3 a = set_verts[0];
        int3 b = set_verts[1];
        Distance3f dist = distance_line3(convert_float3(p),
                                        convert_float3(a),
                                        convert_float3(b));
        Distance3i disti = make_dist3(
            // convert_int_rte(dist.d),
            // convert_int3_rte(dist.p));
            convert_int(dist.d),
            convert_int3(dist.p));
        best = min_pair3i(best, disti);
      }
    }
  }
  return best;
}

NAMESPACE_OCT_END
