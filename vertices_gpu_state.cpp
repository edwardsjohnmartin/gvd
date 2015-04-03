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

#ifdef __OPEN_CL_SUPPORT__
#include "./vertices_gpu_state.h"

namespace oct {

using namespace std;

template <int D>
void VerticesGpuState<D>::InitializeVertices() {
  if (!_vertices) return;

  // typedef Vertex<D> Vertex_t;
  typedef Vertex Vertex_t;

  const int N = _vertices->size();
  _v_array_size = N*2;
  int error;

  // shared_array<Vertex_t> vertex_array(new Vertex_t[_v_array_size]);
  shared_array<Vertex_t> vertex_array(new Vertex_t[N]);
  _vertices->CopyVertices(vertex_array.get());

  // cl_vertex_array =
  //     clCreateBuffer(context_cl, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR,
  //                    sizeof(Vertex_t)*_v_array_size, vertex_array.get(),
  //                    &error);
  cl_vertex_array =
      clCreateBuffer(context_cl, CL_MEM_READ_WRITE,
                     sizeof(Vertex_t)*_v_array_size, 0,
                     &error);
  CheckError(error, "clCreateBuffer");

  cl_command_queue queue_cl = CreateCommandQueue(context_cl, &error);
  CheckError(error, "CreateCommandQueue");

  error = clEnqueueWriteBuffer(
      queue_cl, cl_vertex_array, CL_TRUE, 0,
      sizeof(Vertex_t)*N, vertex_array.get(), 0, nullptr, nullptr);
  CheckError(error, "clEnqueueWriteBuffer");

  error = clFinish(queue_cl);
  CheckError(error, "clFinish");
}

template <int D>
void VerticesGpuState<D>::InitializeClosestPoints() {
  if (!_vertices) return;

  // typedef Vertex<D> Vertex_t;
  typedef Vertex Vertex_t;

  const int M = _vertices->NumClosestPoints();
  int error;

  shared_array<int> point_label_array(new int[M]);
  // _vertices.CopyPointLabels(point_label_array.get());
  for (int i = 0; i < M; ++i) {
    point_label_array[i] = _vertices->LabelAtIdx(i);
  }

  // cl_point_label_array =
  //     clCreateBuffer(context_cl, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR,
  //                    sizeof(int)*M, point_label_array.get(), &error);
  cl_point_label_array =
      clCreateBuffer(context_cl, CL_MEM_READ_WRITE,
                     sizeof(int)*M, 0, &error);
  CheckError(error, "clCreateBuffer");

  cl_command_queue queue_cl = CreateCommandQueue(context_cl, &error);
  CheckError(error, "CreateCommandQueue");

  error = clEnqueueWriteBuffer(
      queue_cl, cl_point_label_array, CL_TRUE, 0,
      sizeof(int)*M, point_label_array.get(), 0, nullptr, nullptr);
  CheckError(error, "clEnqueueWriteBuffer");

  error = clFinish(queue_cl);
  CheckError(error, "clFinish");
}

template <int D>
void VerticesGpuState<D>::Update() {
  if (!_vertices) return;

  if (!_changes->Changed())
    return;

  if (_changes->NumClosestPointsAdded() > 0)
    throw logic_error("Didn't expect the number of closest points to change");

  if (_vertices->size() > _v_array_size) {
    cout << "Reinitializing gpu memory" << endl;
    InitializeVertices();
  } else {
    cout << "Updating gpu memory" << endl;
    // typedef Vertex<D> Vertex_t;
    typedef Vertex Vertex_t;
    // cl_context context = GetCLContext();
    int error;

    // Get values out of vertex network
    const int num_changed = _changes->NumVerticesChanged();
    const int num_added = _changes->NumVerticesAdded();
    const int N = num_changed+num_added;

    shared_array<int> vi_array(new int[num_changed+num_added]);
    _changes->VerticesChanged(vi_array.get());
    for (int i = 0; i < num_added; ++i) {
      vi_array[num_changed+i] = _vertices->size()-num_added+i;
    }
    shared_array<Vertex_t> v_array(new Vertex_t[N]);
    for (int i = 0; i < N; ++i) {
      v_array[i] = _vertices->GetVertex(vi_array[i]);
    }

    cout << N << " vertices are being updated" << endl;

    // Set up the opencl buffers
    cl_mem cl_vi_array =
        clCreateBuffer(context_cl, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR,
                       sizeof(int)*N, vi_array.get(), &error);
    CheckError(error, "clCreateBuffer");
    cl_mem cl_v_array =
        clCreateBuffer(context_cl, CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR,
                       sizeof(Vertex_t)*N, v_array.get(), &error);
    CheckError(error, "clCreateBuffer");

    // Set up the kernel arguments
    clSetKernelArg(
        kernel_update_vertices, 0, sizeof(cl_mem), &cl_vertex_array);
    clSetKernelArg(
        kernel_update_vertices, 1, sizeof(cl_mem), &cl_point_label_array);
    clSetKernelArg(
        kernel_update_vertices, 2, sizeof(cl_mem), &cl_vi_array);
    clSetKernelArg(
        kernel_update_vertices, 3, sizeof(cl_mem), &cl_v_array);
  
    cl_command_queue queue_cl = CreateCommandQueue(context_cl, &error);

    // Enqueue the buffer writes
    error = clEnqueueWriteBuffer(
        queue_cl, cl_vi_array, CL_TRUE, 0,
        sizeof(int)*N, vi_array.get(), 0, nullptr, nullptr);
    CheckError(error, "clEnqueueWriteBuffer");
    error = clEnqueueWriteBuffer(
        queue_cl, cl_v_array, CL_TRUE, 0,
        sizeof(Vertex_t)*N, v_array.get(), 0, nullptr, nullptr);
    CheckError(error, "clEnqueueWriteBuffer");

    // Enqueue the kernel job
    size_t globalWorkSize[] = { N };
    cl_event gpu_event;
    error = clEnqueueNDRangeKernel(
        queue_cl, kernel_update_vertices, 1, nullptr, globalWorkSize,
        nullptr, 0, nullptr, &gpu_event);
    CheckError(error, "clEnqueueNDRangeKernel");

    // Execute
    clFinish(queue_cl);

    // Cleanup
    clReleaseMemObject(cl_vi_array);
    clReleaseMemObject(cl_v_array);
    clReleaseCommandQueue(queue_cl);
  }

  _changes->Reset();

  // Cleanup();
  // InitializeVertices();
  // InitializeClosestPoints();
}

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
// Template instantiations
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

template void VerticesGpuState<2>::InitializeVertices();
template void VerticesGpuState<3>::InitializeVertices();

template void VerticesGpuState<2>::InitializeClosestPoints();
template void VerticesGpuState<3>::InitializeClosestPoints();

template void VerticesGpuState<2>::Update();
template void VerticesGpuState<3>::Update();

}

#endif //__OPEN_CL_SUPPORT__
