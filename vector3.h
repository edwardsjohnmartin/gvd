// This file contains functions useful for octree building from 2D vector
// structures, such as polygons.

#ifndef __OCT_VECTOR3_H__
#define __OCT_VECTOR3_H__

#include <vector>
#include <utility>

#include "./opencl/vec.h"
#include "./opencl/triangle.h"
#include "./opencl/distance3.h"
#include "./opencl/vertex.h"
#include "./shared_ptr.h"
#include "./shared_array.h"

namespace oct {

class LabeledGeometry3 {
 public:
  LabeledGeometry3() : _label(-1) {}
  LabeledGeometry3(const std::vector<Triangle> triangles, int label)
      : _label(label), _triangles(triangles) {
  }
  // LabeledGeometry3(const Triangle* triangles, int count, int label)
  //     : _label(label), _triangles(triangles), _count(count) {
  // }

  int GetLabel() const { return _label; }
  size_t size() const { return _triangles.size(); }
  // size_t size() const { return _count; }
  bool empty() const { return _triangles.empty(); }
  // bool empty() const { return _count == 0; }
  // const Triangle* GetTriangles() const { return _triangles; }
  const std::vector<Triangle>& GetTriangles() const { return _triangles; }
  const std::vector<Triangle>& GetPrimitives() const { return GetTriangles(); }

 private:
  int _label;
  std::vector<Triangle> _triangles;
  // const Triangle* _triangles;
  // const int _count;
};

typedef LabeledGeometry3 LabeledGeometry;

Distance3i distance_geom(
    const int3& p, const LabeledGeometry3& geometry,
    int3* verts);
Distance3i distance_geom(const int3& p, const int* geometry,
                    int3* verts);
PointAndLabel distance_geoms(
    const int3& p, const std::vector<LabeledGeometry3>& geometries,
    int3* verts, int* verts_offsets);
PointAndLabel distance_geoms(
    const int3& p, const int* geometries,
    int3* verts, int* verts_offsets);

void ClipGeometry(
    const LabeledGeometry3& geometry, const int3& base_point,
    const level_t level,
    std::vector<std::vector<LabeledGeometry3> >& clipped,
    int3* vertices);

}

#endif
