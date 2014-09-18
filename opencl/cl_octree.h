#ifndef __CL_OCTREE_H__
#define __CL_OCTREE_H__

#include "./vertex.h"
#include "./geometry.h"
#include "./uvertex_network.h"

NAMESPACE_OCT_BEGIN

#ifndef OPEN_CL
static uchar* g_to_subdivide;
#endif

// // This function sets the subdivision array to false
// inline void InitToSubdivide(int vi, __GLOBAL__ uchar* to_subdivide) {
//   to_subdivide[vi] = false;
// }

// Checks if a cell intersects with multiple geometries
void FindToSubdivide1(
    int vi,
    // int num_vertices, __GLOBAL__ Vertex* vertices,
    __GLOBAL__ int* nv_header, __GLOBAL__ Vertex* vertices,
    // int num_cpoints, __GLOBAL__ GeomPoint* cpoints,
    __GLOBAL__ int* vi2geometries_array,
    int max_level,
    __GLOBAL__ uchar* to_subdivide);

// Checks if a cell's neighbors either have different labels
// (in order to have a buffer region) or are subdivided more
// than kGradation than this cell (in order to have bounded
// gradation).
//
// changed: write-only argument indicating whether any cell's subdivide
// was changed to true.  This is used in iteration to ensure that
// gradation is maintained at every stage.
void FindToSubdivide2(
    int vi,
    // int num_vertices, __GLOBAL__ Vertex* vertices,
    __GLOBAL__ int* nv_header, __GLOBAL__ Vertex* vertices,
    __GLOBAL__ int* vi2nbrs,
    __GLOBAL__ int* vi2geometries_array,
    int max_level,
    __GLOBAL__ uchar* to_subdivide,
    __GLOBAL__ uchar* changed);

// This is to run on a single gpu processor (not in parallel).  It is done
// on the gpu so that we don't have to transfer to_subdivide back to the cpu.
static int CountToSubdivide(
    int num_vertices, __GLOBAL__ uchar* to_subdivide) {
  int count = 0;
  for (int vi = 0; vi < num_vertices; ++vi) {
    if (to_subdivide[vi]) {
      ++count;
    }
  }
  return count;
}

// vi2geometries_array and following arguments can be null.
void SubdivideCell(
    const int vi,
    const int filter,
    __GLOBAL__ int* vn_header,
    __GLOBAL__ Vertex* vertex_array,
    const __GLOBAL__ uchar* to_subdivide);

void SubdivideCell_A(
    const int vi,
    const int filter,
    __GLOBAL__ int* vn_header,
    __GLOBAL__ Vertex* vertex_array,
    const __GLOBAL__ uchar* to_subdivide);
void SubdivideCell_B(
    const int vi,
    const int filter,
    __GLOBAL__ int* vn_header,
    __GLOBAL__ Vertex* vertex_array,
    const __GLOBAL__ uchar* to_subdivide);
void SubdivideCell_C(
    const int vi,
    const int filter,
    __GLOBAL__ int* vn_header,
    __GLOBAL__ Vertex* vertex_array,
    const __GLOBAL__ uchar* to_subdivide);

void ClipGeometriesAfterSubdivide(
    const int vi,
    const int filter,
    __GLOBAL__ int* vn_header,
    __GLOBAL__ Vertex* vertex_array,
    const __GLOBAL__ uchar* to_subdivide,
    __GLOBAL__ int* vi2geometries_array,
    int num_gvertices, __GLOBAL__ intn* gvertices,
    int num_goffsets, __GLOBAL__ int* goffsets);

void ComputeNonEmptyVertexDistances(
    const int vi,
    const int lvi,
    __GLOBAL__ int* vn_header,
    __GLOBAL__ Vertex* vertex_array,
    __GLOBAL__ GeomPoint* cpoints_array,
    __GLOBAL__ int* vi2geometries_array,
    int num_gvertices, __GLOBAL__ intn* gvertices,
    int num_goffsets, __GLOBAL__ int* goffsets);

void PullDistances(
    const int vi,
    __GLOBAL__ int* vn_header,
    __GLOBAL__ Vertex* vertex_array,
    __GLOBAL__ GeomPoint* cpoints_array,
    __GLOBAL__ uchar* changed);

void FindAmbiguous(
    int vi,
    __GLOBAL__ int* nv_header,
    __GLOBAL__ Vertex* vertex_array,
    __GLOBAL__ GeomPoint* cpoints_array,
    int max_level,
    __GLOBAL__ uchar* ambiguous,
    __GLOBAL__ int* count);

NAMESPACE_OCT_END

#endif
