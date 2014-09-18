// This file contains functions useful for octree building from 2D vector
// structures, such as polygons.

#ifndef __OCT_VECTOR2_H__
#define __OCT_VECTOR2_H__

#include <vector>
#include <utility>
#include <cfloat>

#include "./opencl/defs.h"
#include "./opencl/vec.h"
#include "./opencl/edge.h"

namespace oct {

struct Range {
  Range() {
    t[0] = FLT_MAX;
    t[1] = -FLT_MAX;
  }
  Range(float t0, float t1) {
    if (t1 < t0) std::swap(t0, t1);
    t[0] = t0;
    t[1] = t1;
  }
  void operator()(const float& f) {
    t[0] = std::min(t[0], f);
    t[1] = std::max(t[1], f);
  }
  bool In(const float& f) const {
    return f >= t[0] && f <= t[1];
  }
  bool Intersects(const Range& rhs) const {
    return rhs.In(t[0]) || rhs.In(t[1]) ||
        In(rhs.t[0]) || In(rhs.t[1]);
  }
  float t[2];
};

// Data associated with a primitive (simplex in our case) that can
// be reused in intersection computation using the Separating Axis
// Theorem.
struct SATData2 {
  SATData2() {}
  // SATData2(const Edge& e, const std::vector<int2>& vertices) {
  SATData2(const Edge& e, const int2* vertices) {
    verts[0] = vertices[e.s[0]];
    verts[1] = vertices[e.s[1]];
    Init();
  }
  static float2 Ortho(const float2& v) {
    return make_float2(-v.s[1], v.s[0]);
  }
  static float2 Ortho(const int2& v) {
    return make_float2(-v.s[1], v.s[0]);
  }
  static float2 ToFloat(const int2& v) {
    return make_float2(v.s[0], v.s[1]);
  }

  const int2& operator[](int i) const { return verts[i]; }

  void Init() {
    float2 plane_t = Ortho(verts[1] - verts[0]);
    // const float norm = plane_t.norm();
    const float norm = length(plane_t);
    plane = plane_t / norm;
    // range(ToFloat(verts[0]) * plane);
    range(dot(ToFloat(verts[0]), plane));
  }
  int2 verts[2];
  Range range;
  float2 plane;
};

class LabeledGeometry2 {
 public:
  typedef Edge Primitive;
  static const int D;

 public:
  LabeledGeometry2() : _label(-1) {}
  explicit LabeledGeometry2(int label) : _label(label) {}
  LabeledGeometry2(int2* vertices, int label)
      : _vertices(vertices), _label(label) {
  }
  LabeledGeometry2(int2* vertices,
                  const std::vector<Edge> edges, int label)
      : _vertices(vertices), _edges(edges), _axes(edges.size()), _label(label) {
    for (int i = 0; i < _edges.size(); ++i) {
      SetAxis(i);
    }
  }
  LabeledGeometry2(const std::vector<Edge> edges, int label)
      : _edges(edges), /*_axes(edges.size()),*/ _label(label) {
    // for (int i = 0; i < _edges.size(); ++i) {
    //   SetAxis(i);
    // }
  }

  int GetLabel() const { return _label; }
  size_t size() const { return _edges.size(); }
  bool empty() const { return _edges.empty(); }
  int2* GetVertices() const { return _vertices; }
  const std::vector<Edge>& GetEdges() const { return _edges; }
  const std::vector<Edge>& GetPrimitives() const { return GetEdges(); }
  const std::vector<SATData2>& GetAxes() const { return _axes; }
  const Edge& operator[](int i) const { return _edges[i]; }

  void Add(const Edge& e, const SATData2& axis) {
    _edges.push_back(e);
    _axes.push_back(axis);
  }

 private:
  void SetAxis(int i) {
    // _axes[i] = SATData2(_edges[i], *_vertices);
    _axes[i] = SATData2(_edges[i], _vertices);
  }

 public:
  static int count;

 private:
  // mutable oct::shared_ptr<std::vector<int2> > _vertices;
  mutable int2* _vertices;
  std::vector<Edge> _edges;
  // Axes used in decomposition, which uses the Separating Axis Theorem (SAT)
  std::vector<SATData2> _axes;
  int _label;
};

inline std::ostream& operator<<(
    std::ostream& out, const LabeledGeometry2& geom) {
  out << "Label = " << geom.GetLabel()
      << " Vertices = ";
  for (int i = 0; i < geom.GetEdges().size(); ++i) {
    out << (*geom.GetVertices()).s[geom.GetEdges()[i].s[0]] << " ";
  }
  return out;
}

PointAndLabel distance_geoms(
    const int2& p, const std::vector<LabeledGeometry2>& geometries,
    int2* verts, int* verts_offsets);
PointAndLabel distance_geoms(
    const int2& p, const int* geometries,
    int2* verts, int* verts_offsets);

void ClipGeometry(
    const LabeledGeometry2& geometry, const int2& base_point,
    const level_t level,
    std::vector<std::vector<LabeledGeometry2> >& clipped,
    int2* vertices);

// Finds which sub-quadrants the line a-b intersects with
void Intersects(const int2& a, const int2& b,
                const int2& base, const int width,
                bool intersects[]);

// bool Intersects(const int2& base, const int width,
//                 const Triangle& t, const std::vector<int2>& vertices);
bool Intersects(const int2& base, const int width,
                const SATData2& d);
}

#endif
