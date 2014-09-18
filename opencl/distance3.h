#ifndef __DISTANCE3_H__
#define __DISTANCE3_H__

#include "./defs.h"
#include "./triangle.h"
#include "./geometry.h"
#include "./vec.h"

#ifdef OPEN_CL
// OpenCL
#else
// C++
// #define global
// inline int convert_int(const float f) {
//   return (int)(f+0.5);
// }
// inline int convert_int(const double d) {
//   return (int)(d+0.5);
// }
// inline bool any(const bool b) {
//   return b;
// }
#endif

NAMESPACE_OCT_BEGIN

typedef struct {
  float d;
  float3 p;
} Distance3f;

typedef struct {
  int d;
  int3 p;
} Distance3i;

typedef struct {
  int3 p;
  int l;
} PointAndLabel3;

/***************************************************
 * Declarations
 ***************************************************/
int find_unique(const int3* verts, int3* unique_verts);
Distance3f distance_line3(const float3 p, const float3 a, const float3 b);
Distance3i make_dist3(int d, const int3 p);
Distance3f min_pair3f(const Distance3f a, const Distance3f b);
Distance3i min_pair3i(const Distance3i a, const Distance3i b);
Distance3f distance_trianglef(const float3 p, const float3* poly);
Distance3i distance_trianglei(const int3 p, const int3* poly);
Distance3i distance_geom3(
    const int3 p, global const int* geometry, global const int3* verts);
PointAndLabel3 distance_geoms3(
    const int3 p, global const int* geometries,
    global const int3* verts, global const int* verts_offsets);
Distance3i distance_geometry(
    const int3 p, Geometry geometry, __GLOBAL__ const int3* verts);

NAMESPACE_OCT_END

#endif
