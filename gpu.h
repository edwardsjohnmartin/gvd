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

#ifndef __GPU_H__
#define __GPU_H__

#include "./opencl/defs.h"
#include "./opencl/vertex.h"
#include "./opencl/uvertex_network.h"
#include "./opencl/geometry.h"
#include "./opencl.h"
#include "./shared_array.h"

NAMESPACE_OCT_BEGIN

class GpuEvents {
 public:
  GpuEvents()
      : _size(0), _events() {}
  GpuEvents(cl_event e)
      : _size(1), _events(new cl_event[_size]) {
    _events[0] = e;
  }
  GpuEvents(cl_event e0, cl_event e1)
      : _size(2), _events(new cl_event[_size]) {
    _events[0] = e0;
    _events[1] = e1;
  }
  GpuEvents(cl_event e0, cl_event e1, cl_event e2)
      : _size(3), _events(new cl_event[_size]) {
    _events[0] = e0;
    _events[1] = e1;
    _events[2] = e2;
  }
  GpuEvents(GpuEvents events, cl_event e)
      : _size(events.size()+1), _events(new cl_event[_size]) {
    for (int i = 0; i < _size-1; ++i) {
      _events[i] = events[i];
    }
    _events[_size-1] = e;
  }

  int size() const { return _size; }
  cl_event* get() { return _events.get(); }
  cl_event operator[](int i) const { return _events[i]; }

 private:
  int _size;
  shared_array<cl_event> _events;
};

#ifdef __OPEN_CL_SUPPORT__

//------------------------------------------------------------------------------
// Gpu class with OpenCL support
//------------------------------------------------------------------------------

class Gpu {
 public:
  Gpu();
  ~Gpu();

  // Karras calls
  void CreateMPoints(const int n);

  void CreateToSubdivide(int n);
  void CreateVertices(int n);
  void CreateCPoints(int n);
  void CreateVi2Geometries(int n);
  void CreateGeomVertices(GeomVertices geom_vertices);

  // Ensures that the buffer in gpu memory is sufficient to hold n vertices
  void EnsureVertices(
      int threshold, int new_n, GpuEvents wait_events = GpuEvents());
  void EnsureCPoints(GpuEvents wait_events = GpuEvents());

  cl_event EnqueueWriteToSubdivide(
      int n, uchar* array, GpuEvents wait_events = GpuEvents());
  cl_event EnqueueWriteVNHeader(
      int* header, GpuEvents wait_events = GpuEvents());
  cl_event EnqueueWriteVertices(
      int n, Vertex* array, GpuEvents wait_events = GpuEvents());
  cl_event EnqueueWriteCPoints(
      int n, GeomPoint* array, GpuEvents wait_events = GpuEvents());
  cl_event EnqueueWriteVi2Geometries(
      int* array, GpuEvents wait_events = GpuEvents());
  cl_event EnqueueWriteChanged(
      uchar* value, GpuEvents wait_events = GpuEvents());
  cl_event EnqueueWriteSize(
      int* value, GpuEvents wait_events = GpuEvents());
  GpuEvents EnqueueWriteGeomVertices(
      GeomVertices geom_vertices, GpuEvents wait_events = GpuEvents());

  cl_event EnqueueReadToSubdivide(
      int n, uchar* array, GpuEvents wait_events = GpuEvents());
  cl_event EnqueueReadVNHeader(
      int* header, GpuEvents wait_events = GpuEvents());
  cl_event EnqueueReadVertices(
      int n, Vertex* array, GpuEvents wait_events = GpuEvents());
  cl_event EnqueueReadCPoints(
      int n, GeomPoint* array, GpuEvents wait_events = GpuEvents());
  cl_event EnqueueReadVi2Geometries(
      int* array, GpuEvents wait_events = GpuEvents());
  cl_event EnqueueReadChanged(
      uchar* changed, GpuEvents wait_events = GpuEvents());
  cl_event EnqueueReadCount(
      int* count, GpuEvents wait_events = GpuEvents());
  cl_event EnqueueReadSize(
      int* size, GpuEvents wait_events = GpuEvents());

  cl_event EnqueueFindToSubdivide1(
      int num_vertices, int max_level, GpuEvents wait_events = GpuEvents());
  cl_event EnqueueFindToSubdivide2(
      int num_vertices, int max_level, GpuEvents wait_events = GpuEvents());
  cl_event EnqueueCountToSubdivide(
      int num_vertices, GpuEvents wait_events = GpuEvents());
  cl_event EnqueueComputeSparseVi2GeometriesSize(
      GpuEvents wait_events = GpuEvents());
  void MakeSparseVi2Geometries(
      const int new_n, GpuEvents wait_events = GpuEvents());
  void CondenseVi2Geometries(
      GpuEvents wait_events = GpuEvents());
  GpuEvents EnqueueSubdivideCell(
      int num_vertices, GpuEvents wait_events = GpuEvents());
  GpuEvents EnqueueClipGeometries(
      int num_vertices, GpuEvents wait_events = GpuEvents());
  GpuEvents EnqueueComputeNonEmptyVertexDistances(
      int num_vertices, GpuEvents wait_events = GpuEvents());
  cl_event EnqueuePullDistances(
      int num_vertices, GpuEvents wait_events = GpuEvents());
  cl_event EnqueueFindAmbiguous(
      int num_vertices, int max_level, GpuEvents wait_events = GpuEvents());

  void Finish() {
    clFinish(queue);
  }

 private:
  cl_command_queue queue;

  // Update constructor and destructor when adding a buffer.
  
  // Karras
  cl_mem mpoints;
  int mpoints_n;

  cl_mem to_subdivide;
  cl_mem vn_header;
  cl_mem vertices;
  cl_mem vi2nbrs;
  cl_mem cpoints;
  cl_mem vi2geometries;
  cl_mem changed;
  cl_mem count;
  cl_mem size;
  cl_mem gvertices;
  cl_mem goffsets;

  // Array sizes
  int to_subdivide_n;
  int vertices_n;
  int cpoints_n;
  int vi2geometries_n;
  int gvertices_n;
  int goffsets_n;

 public:
  int VerticesSize() const { return vertices_n; }
  int CPointsSize() const { return cpoints_n; }
  int Vi2GeometriesSize() const { return vi2geometries_n; }
};

#else // __OPEN_CL_SUPPORT__

//------------------------------------------------------------------------------
// Gpu class with no OpenCL support
// Throws exceptions
//------------------------------------------------------------------------------

class gpu_error : public std::logic_error {
 public:
  gpu_error() : std::logic_error("OpenCL not supported") {}
};

class Gpu {
 public:
  Gpu() {}
  ~Gpu() {}

  void CreateToSubdivide(int n) { throw gpu_error(); }
  void CreateVertices(int n) { throw gpu_error(); }
  void CreateCPoints(int n) { throw gpu_error(); }
  void CreateVi2Geometries(int n) { throw gpu_error(); }
  void CreateGeomVertices(GeomVertices geom_vertices) { throw gpu_error(); }

  // Ensures that the buffer in gpu memory is sufficient to hold n vertices
  void EnsureVertices(
      int threshold, int new_n, GpuEvents wait_events = GpuEvents()) {
    throw gpu_error(); }
  void EnsureCPoints(GpuEvents wait_events = GpuEvents()) {
    throw gpu_error(); }

  cl_event EnqueueWriteToSubdivide(
      int n, uchar* array, GpuEvents wait_events = GpuEvents()) {
    throw gpu_error(); }
  cl_event EnqueueWriteVNHeader(
      int* header, GpuEvents wait_events = GpuEvents()) {
    throw gpu_error(); }
  cl_event EnqueueWriteVertices(
      int n, Vertex* array, GpuEvents wait_events = GpuEvents()) {
    throw gpu_error(); }
  cl_event EnqueueWriteCPoints(
      int n, GeomPoint* array, GpuEvents wait_events = GpuEvents()) {
    throw gpu_error(); }
  cl_event EnqueueWriteVi2Geometries(
      int* array, GpuEvents wait_events = GpuEvents()) {
    throw gpu_error(); }
  cl_event EnqueueWriteChanged(
      uchar* value, GpuEvents wait_events = GpuEvents()) {
    throw gpu_error(); }
  cl_event EnqueueWriteSize(
      int* value, GpuEvents wait_events = GpuEvents()) {
    throw gpu_error(); }
  GpuEvents EnqueueWriteGeomVertices(
      GeomVertices geom_vertices, GpuEvents wait_events = GpuEvents()) {
    throw gpu_error(); }

  cl_event EnqueueReadToSubdivide(
      int n, uchar* array, GpuEvents wait_events = GpuEvents()) {
    throw gpu_error(); }
  cl_event EnqueueReadVNHeader(
      int* header, GpuEvents wait_events = GpuEvents()) {
    throw gpu_error(); }
  cl_event EnqueueReadVertices(
      int n, Vertex* array, GpuEvents wait_events = GpuEvents()) {
    throw gpu_error(); }
  cl_event EnqueueReadCPoints(
      int n, GeomPoint* array, GpuEvents wait_events = GpuEvents()) {
    throw gpu_error(); }
  cl_event EnqueueReadVi2Geometries(
      int* array, GpuEvents wait_events = GpuEvents()) {
    throw gpu_error(); }
  cl_event EnqueueReadChanged(
      uchar* changed, GpuEvents wait_events = GpuEvents()) {
    throw gpu_error(); }
  cl_event EnqueueReadCount(
      int* count, GpuEvents wait_events = GpuEvents()) {
    throw gpu_error(); }
  cl_event EnqueueReadSize(
      int* size, GpuEvents wait_events = GpuEvents()) {
    throw gpu_error(); }

  cl_event EnqueueFindToSubdivide1(
      int num_vertices, int max_level, GpuEvents wait_events = GpuEvents()) {
    throw gpu_error(); }
  cl_event EnqueueFindToSubdivide2(
      int num_vertices, int max_level, GpuEvents wait_events = GpuEvents()) {
    throw gpu_error(); }
  cl_event EnqueueCountToSubdivide(
      int num_vertices, GpuEvents wait_events = GpuEvents()) {
    throw gpu_error(); }
  cl_event EnqueueComputeSparseVi2GeometriesSize(
      GpuEvents wait_events = GpuEvents()) {
    throw gpu_error(); }
  void MakeSparseVi2Geometries(
      const int new_n, GpuEvents wait_events = GpuEvents()) {
    throw gpu_error(); }
  void CondenseVi2Geometries(
      GpuEvents wait_events = GpuEvents()) {
    throw gpu_error(); }
  GpuEvents EnqueueSubdivideCell(
      int num_vertices, GpuEvents wait_events = GpuEvents()) {
    throw gpu_error(); }
  GpuEvents EnqueueClipGeometries(
      int num_vertices, GpuEvents wait_events = GpuEvents()) {
    throw gpu_error(); }
  GpuEvents EnqueueComputeNonEmptyVertexDistances(
      int num_vertices, GpuEvents wait_events = GpuEvents()) {
    throw gpu_error(); }
  cl_event EnqueuePullDistances(
      int num_vertices, GpuEvents wait_events = GpuEvents()) {
    throw gpu_error(); }
  cl_event EnqueueFindAmbiguous(
      int num_vertices, int max_level, GpuEvents wait_events = GpuEvents()) {
    throw gpu_error(); }

  void Finish() { throw gpu_error(); }

 public:
  int VerticesSize() const { throw gpu_error(); }
  int CPointsSize() const { throw gpu_error(); }
  int Vi2GeometriesSize() const { throw gpu_error(); }
};

#endif // __OPEN_CL_SUPPORT__

NAMESPACE_OCT_END

#endif

