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

#include "./opencl/vertex.h"
#include "./opencl/vector2.cl"
#include "./opencl/distance3.h"
#include "./opencl/cl_ambiguous.h"
#include "./opencl/cl_octree.h"

inline int3 Position3(const int lvi, const int3 base_point, const int width);
inline int2 Position2(const int lvi, const int2 base_point, const int width);

inline int3 Position3(const int lvi, const int3 base_point, const int width) {
  int3 p = base_point;
  if (lvi & 1)
    p.x += width;
  if (lvi & 2)
    p.y += width;
  if (lvi & 4)
    p.z += width;
  return p;
}

inline int2 Position2(const int lvi, const int2 base_point, const int width) {
  int2 p = base_point;
  if (lvi & 1)
    p.x += width;
  if (lvi & 2)
    p.y += width;
  return p;
}

/* Compute the distances for the vertices of non-empty leafs */
kernel void ComputeLeafVertexDistances2(
    global const int2* base_points,
    global const level_t* all_cell_levels,
    global const int* all_geom_array,
    global const int* all_corner_vertices,
    global int2* points,
    global int* labels,
    private const level_t max_level,
    global const int2* all_geom_vertices,
    global const int* all_geom_vertex_offsets) {
  size_t i = get_global_id(0);
  const level_t level = all_cell_levels[i];
  const index_t width = Level2CellWidth(level);
  const int2 base_point = ((int2*)base_points)[i];
  global const int* corner_vertices = all_corner_vertices+i*4;
  const int base_vi = corner_vertices[0];
  int offset = all_geom_array[base_vi];
  global const int* geom_array = all_geom_array+offset;

  for (int lvi = 0; lvi < 4; ++lvi) {
    const int2 v_point = Position2(lvi, base_point, width);
    int2 point;
    int label;
    const PointAndLabel2 dpair =
        distance_geoms2(v_point, geom_array,
                       all_geom_vertices, all_geom_vertex_offsets);
    point = dpair.p;
    label = dpair.l;
    points[i*4+lvi] = point;
    labels[i*4+lvi] = label;
  }
}

/* Compute the distances for the vertices of non-empty leafs */
kernel void ComputeLeafVertexDistances3(
    global const int3* base_points,
    global const level_t* all_cell_levels,
    global const int* all_geom_array,
    global const int* all_corner_vertices,
    global int3* points,
    global int* labels,
    private const level_t max_level,
    global const int3* all_geom_vertices,
    global const int* all_geom_vertex_offsets) {
  size_t i = get_global_id(0);
  const level_t level = all_cell_levels[i];
  const index_t width = Level2CellWidth(level);
  const int3 base_point = base_points[i];
  global const int* corner_vertices = all_corner_vertices+i*8;
  const int base_vi = corner_vertices[0];
  int offset = all_geom_array[base_vi];
  global const int* geom_array = all_geom_array+offset;

  for (int lvi = 0; lvi < 8; ++lvi) {
    const int3 v_point = Position3(lvi, base_point, width);
    int3 point;
    int label;
    const PointAndLabel3 dpair =
        distance_geoms3(v_point, geom_array,
                       all_geom_vertices, all_geom_vertex_offsets);
    point = dpair.p;
    label = dpair.l;
    points[i*8+lvi] = point;
    labels[i*8+lvi] = label;
  }
}

kernel void UpdateVertices3(
    global Vertex* vertices,
    global const int* point_labels,
    global const int* vi_changed,
    global const Vertex* v_changed) {
  const size_t i = get_global_id(0);
  const int vi = vi_changed[i];
  vertices[vi] = v_changed[i];
  /* if (i == 0) */
  /*   printf("Updating vi = %d\n", vi); */
}

kernel void CountAmbiguous(
    global int* ambiguous,
    private const int array_size,
    global int* count) {
  *count = 0;
  int beg = 0;
  int end = array_size-1;
  while (beg < end) {
    // For each ambiguous vertex, replace the value with the index.  If
    // an unambiguous vertex is reached, swap it with the last ambiguous
    // vertex.
    while (beg < end && ambiguous[beg] != 0) {
      ambiguous[beg] = beg;
      ++(*count);
      ++beg;
    }
    // Find last ambiguous vertex
    while (beg < end && ambiguous[end] == 0) {
      ambiguous[end] = -1;
      --end;
    }
    if (beg < end) {
      // swap
      ambiguous[beg] = end;
      ambiguous[end] = -1;
      ++(*count);
      ++beg; --end;
      /* const int temp = ambiguous[beg]; */
      /* ambiguous[beg] = ambiguous[end]; */
      /* ambiguous[end] = temp; */
    }
  }
  /* for (int i = 0; i < array_size; ++i) { */
  /*   if (ambiguous[i] != 0) { */
  /*     ++(*count); */
  /*   } */
  /* } */
}

kernel void k_FindToSubdivide1(
    global int* vn_header, global Vertex* vertices,
    global int* vi2geometries_array,
    int max_level,
    global uchar* to_subdivide) {
  const size_t vi = get_global_id(0);
  FindToSubdivide1(vi, vn_header, vertices,
                   vi2geometries_array, max_level, to_subdivide);
}

kernel void k_FindNbrs(
    global int* vn_header, global Vertex* vertices,
    __GLOBAL__ int* vi2nbrs) {
  const size_t vi = get_global_id(0);
  UVertexNetwork vn = make_uvertex_network(vn_header, vertices, 0);
  if (!IsBase(vi, vn))
    return;
  const int num_nbrs =
      GetNeighborCells(vi, vi2nbrs+vi*kNumIncidentCells, vn);
  if (num_nbrs < kNumIncidentCells)
    vi2nbrs[vi*kNumIncidentCells+num_nbrs] = -1;
}

kernel void k_FindToSubdivide2(
    global int* vn_header, global Vertex* vertices,
    __GLOBAL__ int* vi2nbrs,
    __GLOBAL__ int* vi2geometries_array,
    int max_level,
    __GLOBAL__ uchar* to_subdivide,
    __GLOBAL__ uchar* changed) {
  const size_t vi = get_global_id(0);
  FindToSubdivide2(vi, vn_header, vertices, vi2nbrs,
                   vi2geometries_array, max_level, to_subdivide, changed);
}

kernel void k_CountToSubdivide(
    int num_vertices, __GLOBAL__ uchar* to_subdivide, __GLOBAL__ int* count) {
  const int c = CountToSubdivide(num_vertices, to_subdivide);
  *count = c;
}

kernel void k_CopyVertexArray(
    __GLOBAL__ Vertex* target, __GLOBAL__ Vertex* source, int num_vertices) {
  const size_t vi = get_global_id(0);
  target[vi] = source[vi];
}

kernel void k_ComputeSparseVi2GeometriesSize(
    __GLOBAL__ int* vi2geometries_array,
    __GLOBAL__ uchar* to_subdivide,
    __GLOBAL__ int* size) {
  Vi2Geometries vi2geometries = make_vi2geometries(vi2geometries_array);
  *size = compute_sparse_vi2geometries_size(vi2geometries, to_subdivide);
}

kernel void k_MakeSparseVi2Geometries(
    __GLOBAL__ int* dense_array,
    __GLOBAL__ uchar* candidates,
    __GLOBAL__ int* array,
    int array_size) {
  Vi2Geometries dense = make_vi2geometries(dense_array);
  make_sparse_vi2geometries(dense, candidates, array, array_size);
}

kernel void k_ComputeDenseVi2GeometriesSize(
    __GLOBAL__ int* vi2geometries_array,
    __GLOBAL__ int* size) {
  Vi2Geometries vi2geometries = make_vi2geometries(vi2geometries_array);
  *size = g_array_size(vi2geometries);
}

kernel void k_CondenseVi2Geometries(
    __GLOBAL__ int* sparse_array,
    __GLOBAL__ int* vn_header,
    __GLOBAL__ int* array,
    int array_size) {
  const int N = NumVerticesFromHeader(vn_header);
  Vi2Geometries sparse = make_vi2geometries(sparse_array);
  condense_vi2geometries(sparse, N, array, array_size);
}

kernel void k_SubdivideCell(
    int filter,
    __GLOBAL__ int* vn_header,
    __GLOBAL__ Vertex* vertex_array,
    const __GLOBAL__ uchar* to_subdivide) {
  // const int vi = get_global_id(0);
  // SubdivideCell(
  //     vi, filter, vn_header, vertex_array, to_subdivide);

  // if (vi == 0 && filter == 0)
  //   printf("** Finished subdivide\n");
}

kernel void k_SubdivideCell_A(
    int filter,
    __GLOBAL__ int* vn_header,
    __GLOBAL__ Vertex* vertex_array,
    const __GLOBAL__ uchar* to_subdivide,
    __GLOBAL__ int* size) {
  const int vi = get_global_id(0);
  SubdivideCell_A(
      vi, filter, vn_header, vertex_array, to_subdivide);
  // if (filter == 0) {
  //   atomic_inc(size);
  // }
}

kernel void k_SubdivideCell_B(
    int filter,
    __GLOBAL__ int* vn_header,
    __GLOBAL__ Vertex* vertex_array,
    const __GLOBAL__ uchar* to_subdivide,//) {
    __GLOBAL__ int* size) {
  // const int vi = get_global_id(0);
  // SubdivideCell_B(
  //     vi, filter, vn_header, vertex_array, to_subdivide);
  // if (filter == 0) {
  //   atomic_inc(size);
  // }
}

kernel void k_SubdivideCell_C(
    int filter,
    __GLOBAL__ int* vn_header,
    __GLOBAL__ Vertex* vertex_array,
    const __GLOBAL__ uchar* to_subdivide,//) {
    __GLOBAL__ int* size) {
  // const int vi = get_global_id(0);
  // SubdivideCell_C(
  //     vi, filter, vn_header, vertex_array, to_subdivide);
  // if (filter == 0) {
  //   atomic_inc(size);
  // }
}

// To be called after kSubdivideCell
kernel void k_ClipGeometries(
    int filter,
    __GLOBAL__ int* vn_header,
    __GLOBAL__ Vertex* vertex_array,
    const __GLOBAL__ uchar* to_subdivide,
    __GLOBAL__ int* vi2geometries_array,
    int num_gvertices, __GLOBAL__ intn* gvertices,
    int num_goffsets, __GLOBAL__ int* goffsets) {
  const int vi = get_global_id(0);
  ClipGeometriesAfterSubdivide(
      vi, filter, vn_header, vertex_array, to_subdivide,
      vi2geometries_array,
      num_gvertices, gvertices, num_goffsets, goffsets);
}

kernel void k_GetAmbiguous3_orig(
    global const Vertex* vertices,
    global const int* point_labels,
    global int* ambiguous) {
  const size_t i = get_global_id(0);
  ambiguous[i] = IsAmbiguous3_orig(i, vertices, point_labels);
}

kernel void k_ComputeNonEmptyVertexDistances(
    const int lvi,
    __GLOBAL__ int* vn_header,
    __GLOBAL__ Vertex* vertex_array,
    __GLOBAL__ GeomPoint* cpoints_array,
    __GLOBAL__ int* vi2geometries_array,
    int num_gvertices, __GLOBAL__ intn* gvertices,
    int num_goffsets, __GLOBAL__ int* goffsets) {
  const size_t vi = get_global_id(0);
  ComputeNonEmptyVertexDistances(
      vi, lvi, vn_header, vertex_array, cpoints_array, vi2geometries_array,
      num_gvertices, gvertices, num_goffsets, goffsets);
}

kernel void k_PullDistances(
    __GLOBAL__ int* vn_header,
    __GLOBAL__ Vertex* vertex_array,
    __GLOBAL__ GeomPoint* cpoints_array,
    __GLOBAL__ uchar* changed) {
  const size_t vi = get_global_id(0);
  PullDistances(vi, vn_header, vertex_array, cpoints_array, changed);
}

kernel void k_FindAmbiguous(
    __GLOBAL__ int* vn_header,
    __GLOBAL__ Vertex* vertex_array,
    __GLOBAL__ GeomPoint* cpoints_array,
    int max_level,
    __GLOBAL__ uchar* ambiguous,
    __GLOBAL__ int* count) {
  const size_t vi = get_global_id(0);
  FindAmbiguous(vi, vn_header, vertex_array, cpoints_array,
                max_level, ambiguous, count);
}
