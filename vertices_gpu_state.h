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

#ifndef __VERTICES_GPU_STATE_H__
#define __VERTICES_GPU_STATE_H__

#include "./vertex_network.h"
#include "./opencl/uvertex_network.h"
#include "./opencl.h"

namespace oct {

#ifdef __OPEN_CL_SUPPORT__
template <int D>
class VerticesGpuState {
 public:
  VerticesGpuState(cl_kernel& kernel_update_vertices_,
                   const VertexNetwork& vertices)
      : kernel_update_vertices(kernel_update_vertices_),
        _vertices(&vertices), _changes(new VerticesChangeTracker()) {
    // todo: replace
    // _vertices.SetChangeTracker(_changes);
    InitializeVertices();
    InitializeClosestPoints();
  }
  VerticesGpuState(cl_kernel& kernel_update_vertices_,
                   const UVertexNetwork& vertices)
      : kernel_update_vertices(kernel_update_vertices_), _vertices(0) {
    std::cerr << "VerticesGpuState will do nothing" << std::endl;
  }
  ~VerticesGpuState() {
    Cleanup();
  }

  // Updates GPU if necessary
  const cl_mem* VertexArray() {
    Update();
    return &cl_vertex_array;
  }
  // Updates GPU if necessary
  const cl_mem* PointLabelArray() {
    Update();
    return &cl_point_label_array;
  }

 private:
  void InitializeVertices();
  void InitializeClosestPoints();
  void Update();

  void Cleanup() {
    if (!_vertices) return;
    clReleaseMemObject(cl_vertex_array);
    clReleaseMemObject(cl_point_label_array);
  }

 private:
  cl_kernel& kernel_update_vertices;
  const VertexNetwork* _vertices;
  VerticesChangeTracker_h _changes;
  cl_mem cl_vertex_array;
  cl_mem cl_point_label_array;

  // Sizes of the vertices array on the gpu.  If the size
  // needed is greater than these then a full re-initialization needs to be
  // done.
  int _v_array_size;
};

#else

// Non-opencl implementation
template <int D>
class VerticesGpuState {
 public:
  VerticesGpuState(const VertexNetwork& vertices) {}
  ~VerticesGpuState() {}

  void Update(const VertexNetwork& vertices) {
    // no-op
  }
};

#endif //__OPEN_CL_SUPPORT__

}

#endif
