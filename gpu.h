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

class Gpu {
 public:
  Gpu();
  ~Gpu();

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

NAMESPACE_OCT_END

#endif
