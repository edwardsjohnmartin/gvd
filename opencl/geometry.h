#ifndef __LGEOMETRY_H__
#define __LGEOMETRY_H__

#include "./defs.h"
#include "./vec.h"

NAMESPACE_OCT_BEGIN

//----------------------------------------------------------------------
// Notes
// * N = number of octree cells (base vertices)
// * n = number of geometries for a given cell
// * m = number of faces in a geometry
//----------------------------------------------------------------------

//------------------------------------------------------------
// GeomPoint
// A point on a geometry.  This is just a location with a
// label.  A geometry point can be anywhere on a geometry, not
// just at vertices.
//------------------------------------------------------------
typedef struct {
  // The location on a geometry
  intn p;
  // The label of the geometry
  int l;
} GeomPoint;

inline GeomPoint make_gpoint(intn p, int l) {
  GeomPoint gp = { p, l };
  return gp;
}

//------------------------------------------------------------
// GeomVertices
// Vertices of a geometry mesh.
//
//  0          num labels/offset0
//  1          offset1
//  ...
//  offset0    num vertices
//  offset0+1  vertex0
//  offset0+2  vertex1
//  ...
//------------------------------------------------------------
typedef struct {
  int num_vertices;
  __GLOBAL__ intn* vertices;
  int num_offsets;
  __GLOBAL__ int* offsets;
} GeomVertices;

inline GeomVertices make_geom_vertices(
    int num_vertices, __GLOBAL__ intn* vertices,
    int num_offsets, __GLOBAL__ int* offsets) {
  GeomVertices gv = { num_vertices, vertices, num_offsets, offsets };
  return gv;
}

inline __GLOBAL__ intn* get_geom_vertices(const int label, GeomVertices gv) {
  return gv.vertices + gv.offsets[label];
}

// inline int get_num_geom_vertices(GeomVertices gv) {
//   return gv._num_vertices;
// }

//------------------------------------------------------------
// Geometry struct
//------------------------------------------------------------
// The geometry is stored in memory as a single array as follows:
// 
//  offset    item
//  ------    -----
//    0       label
//    1       count
//    2       face0
//    4       face1
//    6 -     face2-facen
//------------------------------------------------------------
typedef struct {
  __GLOBAL__ Face* faces;
  __GLOBAL__ int* array;
} Geometry;

inline const int g_label(const Geometry g) {
  return g.array[0];
}

inline void g_set_label(int label, const Geometry g) {
  g.array[0] = label;
}

inline const int g_m(const Geometry g) {
  return g.array[1];
}

inline void g_set_m(int m, const Geometry g) {
  g.array[1] = m;
}

inline Geometry make_geometry(__GLOBAL__ int* array) {
  Geometry g = { (__GLOBAL__ Face*)(array+2), array };
  return g;
}

// Returns the number of integers in the array representation
inline int geometry_size(const Geometry geom) {
  // 1 for label
  // 1 for count
  // n * facesize for faces
  assert(g_m(geom) > 0);
  return 2 + g_m(geom) * DIM;
}

//------------------------------------------------------------
// geometries struct
//------------------------------------------------------------
// Each octree cell has a set of labeled geometries.
// The geometries are stored in memory as a single array as follows:
// 
//  offset    item
//  ------    -----
//    0       count/offset0
//    1       offset1
//    2       offset2
//    3 -     offset3-offsetn
//    count   geometry0
//    ...
//------------------------------------------------------------
typedef struct {
  __GLOBAL__ int* offsets;
  __GLOBAL__ int* array;
} Geometries;

inline int g_n(const Geometries g) {
  if (!g.array) return 0;
  return g.array[0];
}

inline void g_set_n(int n, const Geometries g) {
  assert(g.array);
  g.array[0] = n;
}

inline Geometries make_geometries(__GLOBAL__ int* array) {
  Geometries g = { array, array };
  return g;
}

inline Geometry get_geometry(const int i, const Geometries geoms) {
  return make_geometry(geoms.array+geoms.offsets[i]);
}

// Returns the number of integers in the array representation
static int geometries_size(const Geometries geoms) {
  const int n = g_n(geoms);
  int size = n; // offsets
  for (int i = 0; i < n; ++i) {
    size += geometry_size(get_geometry(i, geoms));
  }
  if (size == 0) return 1;
  return size;
}

inline bool cell_is_empty(const Geometries geoms) {
  return g_n(geoms) == 0;
}

inline bool cell_has_multiple_labels(const Geometries geoms) {
  if (g_n(geoms) < 2)
    return false;
  const int label = g_label(get_geometry(0, geoms));
  for (int i = 1; i < g_n(geoms); ++i) {
    if (g_label(get_geometry(i, geoms)) != label)
      return true;
  }
  return false;
}

//------------------------------------------------------------
// Vi2Geometries struct
//------------------------------------------------------------
// Stores the geometries for every cell in the octree
// The geometries are stored in memory as a single array as follows:
// 
//  offset    item
//  ------    -----
//    0       array size
//    1       free_offset
//    2       count/offset0
//    3       offset1
//    4       offset2
//    5 -     offset3-offsetn
//    offset0 cell_geometry0
//    ...
//
// If an offset is -1, then the vertex is either not a base vertex
// or the cell is empty.  Cell geometries are not necessarily in
// order.  That is, offseti may be greater than offsetj for i < j.
//------------------------------------------------------------
typedef struct {
  __GLOBAL__ int* offsets;
  __GLOBAL__ int* array;
} Vi2Geometries;

static __CONST__ int vi2g_offset = 2;

inline int g_array_size(const Vi2Geometries g) {
  assert(g.array);
  return g.array[0];
}

inline int g_free_offset(const Vi2Geometries g) {
  assert(g.array);
  return g.array[1];
}

inline int g_set_free_offset(const int free_offset, const Vi2Geometries g) {
  assert(g.array);
  // g.array[1] = free_offset;
#ifdef OPEN_CL
  const int ret = atomic_xchg(&(g.array[1]), free_offset);
#else
  const int ret = g.array[1];
  g.array[1] = free_offset;
#endif
  return ret;
}

inline int g_add_free_offset(const int val, const Vi2Geometries g) {
  assert(g.array);
#ifdef OPEN_CL
  const int ret = atomic_add(&(g.array[1]), val);
#else
  const int ret = g.array[1];
  g.array[1] += val;
#endif
  return ret;
}

// N is the number of vertices
static const int g_N(const Vi2Geometries g) {
  assert(g.array);
  // for (int vi = 0; vi < g.array_size - vi2g_offset; ++vi) {
  for (int vi = 0; vi < g_array_size(g) - vi2g_offset; ++vi) {
    if (g.offsets[vi] != -1)
      return g.offsets[vi] - vi2g_offset;
  }
  // No cell is non-empty
  assert(0);
  return 0;
}

inline bool gcell_is_empty(
    const int vi, const Vi2Geometries geoms) {
  return geoms.offsets[vi] == -1;
}

inline Geometries get_geometries(
    const int vi, const Vi2Geometries vi2g) {
  if (gcell_is_empty(vi, vi2g)) {
    assert(0);
    // Geometries g = { 0, 0 };
    // return g;
  }
  return make_geometries(vi2g.array+vi2g.offsets[vi]);
}

// Precondition: the array is in proper format
inline Vi2Geometries make_vi2geometries(__GLOBAL__ int* array) {
  Vi2Geometries g = { array+vi2g_offset, array };
  return g;
}

// Returns the number of integers in the array representation
static int vi2geometries_size(const Vi2Geometries geoms) {
  const int n = g_N(geoms);
  int size = n; // offsets
  for (int vi = 0; vi < n; ++vi) {
    if (!gcell_is_empty(vi, geoms))
      size += geometries_size(get_geometries(vi, geoms));
  }
  return size;
}

// Suppose m cells are candidates for subdivision.  Then each of those
// cells must have 2^D * geometries_size(cell[i]) space to write
// new geometries to.
//
// candidates are the cells that are candidates for subdivision.
int compute_sparse_vi2geometries_size(
    const Vi2Geometries dense_geoms,
    __GLOBAL__ uchar* candidates);

// Suppose m cells are candidates for subdivision.  Then the number of
// new cells will be M = m * (2^D - 1).  Each of these new cells is allocated
// space at the end of the array, with enough space to have the same
// geometries as its parent cell.
//
// candidates is of size dense_geoms.n (the number of vertices) and
// candidates[vi] is true if vi is a candidate for subdivision.
Vi2Geometries make_sparse_vi2geometries(
    const Vi2Geometries dense_geoms,
    __GLOBAL__ uchar* candidates,
    __GLOBAL__ int* array,
    int array_size);

// Given a sparse array of octree geometries, defragment (condense)
// so that all free space is at the end.  This is done by writing to array.
// N is the current number of vertices.  Truncate the number of
// vertices in geoms to N.
Vi2Geometries condense_vi2geometries(
    const Vi2Geometries sparse, const int N, __GLOBAL__ int* array,
    int array_size);

void allocate(const Geometries source, const int vi,
              Vi2Geometries vi2geometries);

//------------------------------------------------------------
// ClipGeometries
//------------------------------------------------------------
void ClipGeometries(
    const Geometries source, Geometries target, const int3 base_point,
    const index_t cell_width, const GeomVertices geom_vertices);

// inline void ClipGeometries1(
//     int* source_array, int* target_array, const int3 base_point,
//     const index_t cell_width, 
//     int gv_num_vertices, intn* gv_vertices,
//     int gv_num_offsets, int* gv_offsets) {
//   Geometries source = make_geometries(source_array);
//   Geometries target = make_geometries(target_array);
//   GeomVertices gv = make_geom_vertices(
//       gv_num_vertices, gv_vertices, gv_num_offsets, gv_offsets);
//   ClipGeometries(source, target,
//                  base_point, cell_width, gv);
// }

NAMESPACE_OCT_END

#endif
