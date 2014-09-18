#ifndef __VECTOR2_CL__
#define __VECTOR2_CL__

#include "./edge.h"

typedef struct {
  float d;
  float2 p;
} Distance2f;

typedef struct {
  int d;
  int2 p;
} Distance2i;

typedef struct {
  int2 p;
  int l;
} PointAndLabel2;


/***************************************************
 * Declarations
 ***************************************************/
PointAndLabel2 distance_geoms2(
    /* const Vec2i p, global const int* geometries, */
    const int2 p, global const int* geometries,
    /* global const Vec2i* verts, global const int* verts_offsets); */
    global const int2* verts, global const int* verts_offsets);
Distance2i distance_geom2(
    /* const int2 p, global const int* geometry, global const Vec2i* verts); */
    const int2 p, global const int* geometry, global const int2* verts);
Distance2i distance_line2i(const int2 p, const int2 a, const int2 b);
Distance2f distance_line2f(const float2 p, const float2 a, const float2 b);
Distance2i min_pair2(const Distance2i a, const Distance2i b);
Distance2i make_dist2i(int d, const int2 p);
Distance2f make_dist2f(float d, const float2 p);

/***************************************************
 * Definitions
 ***************************************************/

/* Vec2i to_vec2i(const int2 p) { */
/*   return make_vec2i(p.x, p.y); */
/* } */

/* int2 to_int2(const Vec2i p) { */
/*   return (int2)(p.x[0], p.x[1]); */
/* } */

Distance2f make_dist2f(float d, const float2 p) {
  Distance2f dist = { d, p };
  return dist;
}

Distance2i make_dist2i(int d, const int2 p) {
  Distance2i dist = { d, p };
  return dist;
}

Distance2i min_pair2(const Distance2i a, const Distance2i b) {
  if (a.d < b.d) return a;
  return b;
}

Distance2f distance_line2f(const float2 p, const float2 a, const float2 b) {
  const float dotf = dot(p-a, normalize(b-a));
  if (dotf < 0) {
    return make_dist2f(fast_length(p-a), a);
  } else if (dotf > fast_length(b-a)) {
    return make_dist2f(fast_length(p-b), b);
  } else {
    const float dist = fast_length(
        cross((float3)(p-a, 0), normalize((float3)(b-a, 0))));
    const float2 pi = a + normalize(b-a)*dotf;
    return make_dist2f(dist, pi);
  }
}

Distance2i distance_line2i(const int2 p, const int2 a, const int2 b) {
  Distance2f d = distance_line2f(convert_float2(p),
                                 convert_float2(a),
                                 convert_float2(b));
  return make_dist2i(convert_int_rte(d.d),
                     convert_int2_rte(d.p));
}

Distance2i distance_geom2(
    const int2 p, global const int* geometry, global const int2* verts) {
  Distance2i best = { INT_MAX, int2() };
  const int num_edges = geometry[1];
  const Edge* edges = (Edge*)(geometry+2);
  for (int i = 0; i < num_edges; ++i) {
    const Edge e = edges[i];
    best = min_pair2(best, distance_line2i(
        p, verts[e.s[0]], verts[e.s[1]]));
  }
  return best;
}

PointAndLabel2 distance_geoms2(
    const int2 p_, global const int* geometries,
    global const int2* verts, global const int* verts_offsets) {

  const int2 p = p_;
  Distance2i min_dist = { INT_MAX, int2() };
  int idx = -1;
  for (int i = 0; i < geometries[0]; ++i) {
    int offset = geometries[i];
    const Distance2i dist = distance_geom2(
        p, geometries+offset, verts+verts_offsets[geometries[offset]]);
    if (dist.d < min_dist.d) {
      min_dist = dist;
      idx = i;
    }
  }
  int offset = geometries[idx];
  const int label = geometries[offset];
  PointAndLabel2 pl = { min_dist.p, label };
  /* PointAndLabel2 pl = { (int2)(geometries[3], geometries[1]), label }; */
  return pl;
}


#endif
