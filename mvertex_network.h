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

#ifndef __MVERTEX_NETWORK_H__
#define __MVERTEX_NETWORK_H__

//----------------------------------------------------------------------
// Managed and unmanaged vertex network structs
//----------------------------------------------------------------------

#include <cstring>

#include "./opencl/uvertex_network.h"
#include "./shared_array.h"

NAMESPACE_OCT_BEGIN

//----------------------------------------
// struct MVertexNetwork
//----------------------------------------
struct MVertexNetwork {
  shared_array<int> header;
  // int _num_vertices;
  int vertex_array_capacity;
  shared_array<Vertex> vertices;

  // int _num_cpoints;
  int cpoint_array_capacity;
  shared_array<GeomPoint> cpoints;
};

//------------------------------------------------------------
// Header functions
//------------------------------------------------------------

inline int NumVertices(MVertexNetwork vn) {
  // return vn._num_vertices;
  return vn.header[0];
}

inline int NumCPoints(MVertexNetwork vn) {
  // return vn._num_cpoints;
  return vn.header[1];
}

inline void SetNumVertices(int n, MVertexNetwork* vn) {
  // vn->_num_vertices = n;
  vn->header[0] = n;
}

inline void SetNumCPoints(int n, MVertexNetwork* vn) {
  // vn->_num_cpoints = n;
  vn->header[1] = n;
}

static void mvn_reserve(int vsize, MVertexNetwork* mvn) {
  if (vsize > mvn->vertex_array_capacity) {
    Vertex* new_array = new Vertex[vsize];
    // memcpy(new_array, mvn->vertices.get(),
    //        NumVertices(*mvn) * sizeof(Vertex));
    memcpy(new_array, mvn->vertices.get(),
           mvn->vertex_array_capacity * sizeof(Vertex));
    mvn->vertices.reset(new_array);
    mvn->vertex_array_capacity = vsize;
  }
}

static void mvn_cp_reserve(int cpsize, MVertexNetwork* mvn) {
  if (cpsize > mvn->cpoint_array_capacity) {
    GeomPoint* new_array = new GeomPoint[cpsize];
    // memcpy(new_array, mvn->cpoints.get(),
    //        mvn->num_cpoints * sizeof(GeomPoint));
    if (mvn->cpoints) {
      memcpy(new_array, mvn->cpoints.get(),
             NumCPoints(*mvn) * sizeof(GeomPoint));
    }
    mvn->cpoints.reset(new_array);
    mvn->cpoint_array_capacity = cpsize;
  }
}

// // Ensures that there will be enough space if every one of the
// // num_candidates vertices are subdivided and that each subdivision
// // yields kSubAdded new octree vertices.
// static void mvn_ensure(const int num_candidates, MVertexNetwork* mvn) {
//   // if (mvn->vertex_array_capacity - mvn->num_vertices <
//   //     num_candidates*kSubAdded) {
//   //   mvn_reserve(mvn->num_vertices+num_candidates*kSubAdded*3, mvn);
//   // }
//   if (mvn->vertex_array_capacity - NumVertices(*mvn) <
//       num_candidates*kSubAdded) {
//     mvn_reserve(NumVertices(*mvn)+num_candidates*kSubAdded*3, mvn);
//   }
// }

// Ensures that there will be enough space if every one of the
// num_candidates vertices are subdivided and that each subdivision
// yields kSubAdded new octree vertices.
static void mvn_ensure(
    const int threshold, const int new_n, MVertexNetwork* mvn) {
  if (mvn->vertex_array_capacity < threshold) {
    mvn_reserve(new_n, mvn);
  }
}

//------------------------------------------------------------
// Construction
//------------------------------------------------------------

MVertexNetwork make_mvertex_network();
UVertexNetwork make_vertex_network(MVertexNetwork mvn);
void update_mvertex_network(UVertexNetwork vn, MVertexNetwork& mvn);

NAMESPACE_OCT_END

#endif
