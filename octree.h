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

#ifndef __JME_OCTREE_H__
#define __JME_OCTREE_H__

#include <memory.h>
#include <math.h>
#include <climits>
#include <vector>
#include <list>
#include <iostream>
#include <fstream>
#include <stdexcept>
#include <utility>
#include <set>
#include <algorithm>

#include "./opencl/defs.h"
#include "./opencl/bit.h"
#include "./opencl/vec.h"
#include "./opencl/triangle.h"
#include "./opencl/geometry.h"
#include "./opencl/uvertex_network.h"

#include "./options.h"
#include "./timer.h"
#include "./vertex_network.h"
#include "./shared_ptr.h"
#include "./bb.h"
#include "./opencl.h"
#include "./vectorn.h"

namespace oct {

// See comments to Index() in bit.h.
//
// kFaceIndices[D][dim][i] where D is the dimension of the space,
// dim is the face dimension (0 = point, 1 = edge, 2 = face, 3 = 4D face),
// and i is the index of the face.
//
// Suppose we're in 3D.  To get the 6 indices of all 2D faces,
//    const int D = 3;
//    const int dim = 2;
//    const int num_faces = kNumFaceIndices[dim];
//    for (int i = 0; i < num_faces; ++i) {
//      cout << kFaceIndices[D][dim][i] << endl;
//    }
extern const int kFaceIndices[4][4][12];
// kNumFaceIndices[D][dim]
extern const int kNumFaceIndices[4][4];

extern const int kNumSubdividedArray[4];
// Number of vertices in a subdivided cell.
//   2D: 9
//   3D: 27
// static const int kNumSubdivided = kNumSubdividedArray[DIM];

// Usage:
//  Ensure(vec, i);
//  vec[i] = a;
template <typename Vector>
void Ensure(Vector& v, int idx) {
  if (v.size() <= idx)
    v.resize(idx+1);
}

// Usage:
//  Ensure(vec, i, default_val);
//  vec[i] = a;
template <typename Vector, typename T>
void Ensure(Vector& v, int idx, const T& default_val) {
  if (v.size() <= idx)
    v.resize(idx+1, default_val);
}

//------------------------------------------------------------------------------
// Position
//
// Position of subcell lvi's 0 index vertex (current cell 0 idx at base_point
// and width 2*width)
//------------------------------------------------------------------------------
template <typename Point>
Point Position(const int lvi, const Point& base_point, const int width) {
  Point p(base_point);
  for (int j = 0; j < DIM; ++j) {
    if (lvi & (1 << j)) {
      p.s[j] += width;
    }
  }
  return p;
}


//------------------------------------------------------------------------------
// Position
//
// Position of vertex in direction d (current cell center at base_point
// and has width 2*width)
//------------------------------------------------------------------------------
template <typename Point>
Point Position(const int pos, const int neg, const Point& base_point,
               const index_t width) {
  static const int D = DIM;
  Point p(base_point);
  for (int i = 0; i < D; ++i) {
    const int mask = (1<<i);
    // if (d.Pos() & mask) {
    if (pos & mask) {
      p.s[i] += width;
      // } else if (d.Neg() & mask) {
    } else if (neg & mask) {
      p.s[i] -= width;
    }
  }
  return p;
}

//------------------------------------------------------------------------------
// HeapVertex
//------------------------------------------------------------------------------
struct HeapVertex {
  HeapVertex(const int dist_, const int vi_)
      : dist(dist_), vi(vi_) {}
  bool operator<(const HeapVertex& rhs) const {
    if (dist == rhs.dist) {
      return vi < rhs.vi;
    }
    return dist < rhs.dist;
  }
  friend std::ostream& operator<<(std::ostream& out, const HeapVertex& v) {
    out << "dist = " << v.dist << " vi = " << v.vi;
    return out;
  }
  int dist;
  int vi;
};

//------------------------------------------------------------------------------
// ComputeCornerDistances
//
// Compute distances for each corner vertex of a cell
//------------------------------------------------------------------------------
template <int D, typename LabeledGeometry>
void ComputeCornerDistances(
    const intn& base_point,
    const level_t level,
    const std::vector<LabeledGeometry>& cell_geometries,
    VertexNetwork& vertices,
    const int corner_vertices[],
    const GeomVertices& geom_vertices,
    std::multiset<HeapVertex>& heap,
    const OctreeOptions& o) {

  const level_t max_level = o.max_level;
  const index_t width = Level2CellWidth(level);

  for (int lvi = 0; lvi < (1<<D); ++lvi) {
    const intn v_point = Position(lvi, base_point, width);
    const bool simple_dist = o.simple_dist;
    static const int kSimpleThreshold = 2;
    const int vi = corner_vertices[lvi];
    intn point = make_intn();
    int label;
    int dist;
    if (simple_dist &&
        cell_geometries.size() == 1 && max_level - level < kSimpleThreshold) {
      point = v_point;
      label = cell_geometries[0].GetLabel();
      dist = 0;
    } else {
      const PointAndLabel dpair = distance_geoms(
          v_point, cell_geometries,
          // geom_vertices.all_geom_vertices.get(),
          // geom_vertices.all_geom_vertex_offsets.get());
          // geom_vertices.all_geom_vertices,
          // geom_vertices.all_geom_vertex_offsets);
          geom_vertices.vertices,
          geom_vertices.offsets);
      point = dpair.p;
      label = dpair.l;
      // dist = (point - v_point).norm();
      dist = length(point - v_point);
    }

    const int curve_pointi = vertices.CreateNewClosestPoint(point, label);
    const bool closest_changed = vertices.SetClosestPoint(vi, curve_pointi);
    if (closest_changed) {
      heap.insert(HeapVertex(dist, vi));
    }
  }
}

//------------------------------------------------------------------------------
// Geometry/array conversion routines
//------------------------------------------------------------------------------

template <typename LabeledGeometry>
size_t GeometryArraySize(std::vector<LabeledGeometry>& geoms);

template <typename LabeledGeometry>
size_t GeometryArraySize(std::vector<std::vector<LabeledGeometry> >& geoms);

template <typename LabeledGeometry>
std::vector<LabeledGeometry> ConvertArrayToGeometries(
    const int* geom_array);

template <typename LabeledGeometry>
int ConvertGeometriesToArray(
    std::vector<LabeledGeometry>& geoms, int* geom_array);

template <typename LabeledGeometry>
shared_array<int> ConvertGeometriesToArray(
    std::vector<LabeledGeometry>& geoms);

template <typename LabeledGeometry>
std::vector<std::vector<LabeledGeometry> > ConvertArrayToGeometries2(
    const int* geom_array);

template <typename LabeledGeometry>
shared_array<int> ConvertGeometriesToArray(
    std::vector<std::vector<LabeledGeometry> >& geoms);

//------------------------------------------------------------------------------
// ComputeCornerDistances
//
// Compute distances for each corner vertex of a cell
//------------------------------------------------------------------------------
template <int D>
void ComputeCornerDistancesParallel3(
    const int i, // in opencl use the work unit id
    const intn* base_points,
    const level_t* all_levels,
    const int* all_geom_array,
    const int* all_corner_vertices,
    intn* points,
    int* labels,
    const level_t max_level,
    intn* all_geom_vertices,
    int* all_geom_vertex_offsets) {

  const level_t level = all_levels[i];
  const index_t width = Level2CellWidth(level);
  const int kNumCorners = (1 << D);

  const intn base_point = base_points[i];
  const int* corner_vertices = all_corner_vertices+i*kNumCorners;
  const int base_vi = corner_vertices[0];
  const int* geom_array = all_geom_array+all_geom_array[base_vi];

  for (int lvi = 0; lvi < kNumCorners; ++lvi) {
    const intn v_point = Position(lvi, base_point, width);
    intn point = make_intn();
    int label;
    const PointAndLabel dpair = distance_geoms(
        v_point, geom_array,
        all_geom_vertices, all_geom_vertex_offsets);
    point = dpair.p;
    label = dpair.l;
    points[i*kNumCorners+lvi] = point;
    labels[i*kNumCorners+lvi] = label;
  }
}

template <typename LabeledGeometry>
bool HasMultipleLabels(const std::vector<LabeledGeometry>& cell_geometries) {
  assert (!cell_geometries.empty());
  const int const_label = cell_geometries[0].GetLabel();
  for (int i = 1; i < cell_geometries.size(); ++i) {
    if (const_label != cell_geometries[i].GetLabel()) {
      return true;
    }
  }
  return false;
}

template <int D>
vector<int> GetNeighborsToSubdivide(
    const int vi,
    const intn& base_point,
    const char level,
    VertexNetwork& vertices);

// Conflicting cells both have geometry with differing labels.
template <typename LabeledGeometry>
bool Conflict(
    const int vi0, const int vi1,
    std::vector<std::vector<LabeledGeometry> >& base2geometries) {
  if (base2geometries[vi0].empty() || base2geometries[vi1].empty())
    return false;
  if (HasMultipleLabels(base2geometries[vi0]) ||
      HasMultipleLabels(base2geometries[vi1]))
    return true;
  return
      base2geometries[vi0][0].GetLabel() != base2geometries[vi1][0].GetLabel();
}

// Conflicting cells both have geometry with differing labels.
static bool Conflict(
    const int vi0, const int vi1, Vi2Geometries vi2geometries) {
  Geometries g0 = get_geometries(vi0, vi2geometries);
  Geometries g1 = get_geometries(vi1, vi2geometries);
  if (cell_is_empty(g0) || cell_is_empty(g1))
    return false;
  if (cell_has_multiple_labels(g0) || cell_has_multiple_labels(g1))
    return true;
  Geometry gg0 = get_geometry(0, g0);
  // todo: fix this
  for (int i = 1; i < g_n(g0) && g_m(gg0) == 0; ++i) {
    gg0 = get_geometry(i, g0);
  }
  Geometry gg1 = get_geometry(0, g1);
  // todo: fix this
  for (int i = 1; i < g_n(g1) && g_m(gg1) == 0; ++i) {
    gg1 = get_geometry(i, g1);
  }
  return g_label(gg0) != g_label(gg1);
}

// Subdivides cell vi and splits geometries.  Adds all new cells to
// the output iterator new_cells.
// Returns the new subcells with lvi indexing (see bit.h).
// If slvi2vi isn't needed, pass in null.
// template <typename LabeledGeometry, typename NewCellsIter>
template <typename NewCellsIter>
shared_array<int> SubdivideCell(
    const int vi,
    std::vector<std::vector<LabeledGeometry> >& base2geometries,
    const GeomVertices& geom_vertices,
    int* slvi2vi,
    vector<set<int> >& adjacent_cells,
    NewCellsIter new_cells,
    VertexNetwork& vertices) {

  const intn base_point = vertices.Position(vi);
  // cout << "Subdividing cell. vi = " << vi << " base_point = " << base_point << endl;
  const level_t level = vertices.CellLevel(vi);
  std::vector<LabeledGeometry>& cell_geometries = base2geometries[vi];

  // Subdivide
  shared_array<int> subcells = vertices.Subdivide(vi, slvi2vi);

  // Clip geometries
  std::vector<std::vector<LabeledGeometry> > sub_geometries(1<<DIM);
  for (int i = 0; i < cell_geometries.size(); ++i) {
    // ClipGeometry(cell_geometries[i], base_point, level, sub_geometries,
    //              geom_vertices.Vertices(cell_geometries[i]));
    // ClipGeometry(cell_geometries[i], base_point, level, sub_geometries,
    //              geom_vertices.Vertices(cell_geometries[i].GetLabel()));
    ClipGeometry(
        cell_geometries[i], base_point, level, sub_geometries,
        get_geom_vertices(cell_geometries[i].GetLabel(), geom_vertices));
  }

  // Assign geometries to each sub-cell and add sub-cells to output iterator.
  for (int lvi = 0; lvi < (1<<DIM); ++lvi) {
    const int sub_vi = subcells[lvi];

    Ensure(base2geometries, sub_vi);
    base2geometries[sub_vi] = sub_geometries[lvi];

    *new_cells++ = sub_vi;
  }


  const set<int> adjacent_to_vi = adjacent_cells[vi];

  // Remove all vi -> n_vi
  adjacent_cells[vi].clear();
  // Remove all n_vi -> vi
  for (set<int>::const_iterator it = adjacent_to_vi.begin();
       it != adjacent_to_vi.end(); ++it) {
    adjacent_cells[*it].erase(vi);
  }

  for (int lvi = 0; lvi < (1<<DIM); ++lvi) {
    Ensure(adjacent_cells, subcells[lvi]);
  }

  for (int lvi = 0; lvi < (1<<DIM); ++lvi) {
    const int sub_vi = subcells[lvi];
    if (base2geometries[sub_vi].empty())
      continue;
      
    // Set adjacencies of subcells
    for (int lvi2 = lvi+1; lvi2 < (1<<DIM); ++lvi2) {
      if (Conflict(subcells[lvi], subcells[lvi2], base2geometries)) {
        adjacent_cells[subcells[lvi]].insert(subcells[lvi2]);
        adjacent_cells[subcells[lvi2]].insert(subcells[lvi]);
      }
    }

    // Update links to neighbors that were adjacent to vi
    for (set<int>::const_iterator it = adjacent_to_vi.begin();
         it != adjacent_to_vi.end(); ++it) {
      const int n_vi = *it;
      if (vertices.AreAdjacent(sub_vi, n_vi)) {
        if (Conflict(sub_vi, n_vi, base2geometries)) {
          adjacent_cells[sub_vi].insert(n_vi);
          adjacent_cells[n_vi].insert(sub_vi);
        }
      }
    }
  }

  return subcells;
}

template <typename NewCellsIter>
shared_array<int> SubdivideCellGpu(
    const int vi,
    std::vector<std::vector<LabeledGeometry> >& base2geometries,
    const GeomVertices& geom_vertices,
    int* slvi2vi,
    vector<set<int> >& adjacent_cells,
    NewCellsIter new_cells,
    UVertexNetwork& vertices) {

  // const intn base_point = vertices.Position(vi);
  const intn base_point = vertices.vertices[vi].position;
  // const level_t level = vertices.CellLevel(vi);
  const level_t level = CellLevel(vi, vertices);
  std::vector<LabeledGeometry>& cell_geometries = base2geometries[vi];

  // Subdivide
  // shared_array<int> subcells = Subdivide(vi, slvi2vi, &vertices);
  shared_array<int> subcells(new int[1<<DIM]);
  // Subdivide(vi, subcells.get(), slvi2vi, &vertices);
  Subdivide(vi, slvi2vi, &vertices);
  memcpy(subcells.get(), vertices.vertices[vi].corners, sizeof(int)*(1<<DIM));

  // Clip geometries
  std::vector<std::vector<LabeledGeometry> > sub_geometries(1<<DIM);
  for (int i = 0; i < cell_geometries.size(); ++i) {
    // ClipGeometry(cell_geometries[i], base_point, level, sub_geometries,
    //              geom_vertices.Vertices(cell_geometries[i]));
    // ClipGeometry(cell_geometries[i], base_point, level, sub_geometries,
    //              geom_vertices.Vertices(cell_geometries[i].GetLabel()));
    ClipGeometry(
        cell_geometries[i], base_point, level, sub_geometries,
        get_geom_vertices(cell_geometries[i].GetLabel(), geom_vertices));
  }

  // Assign geometries to each sub-cell and add sub-cells to output iterator.
  for (int lvi = 0; lvi < (1<<DIM); ++lvi) {
    const int sub_vi = subcells[lvi];

    Ensure(base2geometries, sub_vi);
    base2geometries[sub_vi] = sub_geometries[lvi];

    *new_cells++ = sub_vi;
  }


  const set<int> adjacent_to_vi = adjacent_cells[vi];

  // Remove all vi -> n_vi
  adjacent_cells[vi].clear();
  // Remove all n_vi -> vi
  for (set<int>::const_iterator it = adjacent_to_vi.begin();
       it != adjacent_to_vi.end(); ++it) {
    adjacent_cells[*it].erase(vi);
  }

  for (int lvi = 0; lvi < (1<<DIM); ++lvi) {
    Ensure(adjacent_cells, subcells[lvi]);
  }

  for (int lvi = 0; lvi < (1<<DIM); ++lvi) {
    const int sub_vi = subcells[lvi];
    if (base2geometries[sub_vi].empty())
      continue;
      
    // Set adjacencies of subcells
    for (int lvi2 = lvi+1; lvi2 < (1<<DIM); ++lvi2) {
      if (Conflict(subcells[lvi], subcells[lvi2], base2geometries)) {
        adjacent_cells[subcells[lvi]].insert(subcells[lvi2]);
        adjacent_cells[subcells[lvi2]].insert(subcells[lvi]);
      }
    }

    // Update links to neighbors that were adjacent to vi
    for (set<int>::const_iterator it = adjacent_to_vi.begin();
         it != adjacent_to_vi.end(); ++it) {
      const int n_vi = *it;
      // if (vertices.AreAdjacent(sub_vi, n_vi)) {
      if (AreAdjacent(sub_vi, n_vi, vertices)) {
        if (Conflict(sub_vi, n_vi, base2geometries)) {
          adjacent_cells[sub_vi].insert(n_vi);
          adjacent_cells[n_vi].insert(sub_vi);
        }
      }
    }
  }

  return subcells;
}

// // Subdivides cell vi and splits geometries.  Adds all new cells to
// // the output iterator new_cells.
// // Returns the new subcells with lvi indexing (see bit.h).
// // If slvi2vi isn't needed, pass in null.
// // template <typename LabeledGeometry, typename NewCellsIter>
// template <typename NewCellsIter>
// shared_array<int> SubdivideCellGpu(
//     const int vi,
//     Vi2Geometries& vi2geometries,
//     const GeomVertices& geom_vertices,
//     int* slvi2vi,
//     vector<set<int> >& adjacent_cells,
//     NewCellsIter new_cells,
//     UVertexNetwork& vertices);

//------------------------------------------------------------------------------
// PushDistance
//
// Pushes the closest point of vi to its neighbor in direction d.  Returns
// true if the closest point of the neighbor changed.
//------------------------------------------------------------------------------
template <int D>
HeapVertex PushDistance(
    VertexNetwork& vertices, const int vi, const Direction& d,
    const intn& v_point, const int label,
    const OctreeOptions& o) {

  const intn curve_point = vertices.ClosestPoint(vi);
  const int curve_pointi = vertices.ClosestPointIndex(vi);

  const int n_vi = vertices.Neighbor(vi, d);
  const intn n_v_point = vertices.Position(n_vi);
  const int n_dist = length(n_v_point - curve_point);

  const bool closest_changed = vertices.SetClosestPoint(n_vi, curve_pointi);
  if (closest_changed) {
    return HeapVertex(n_dist, n_vi);
  }
  return HeapVertex(-1, -1);
}

//------------------------------------------------------------------------------
// PushDistance
//
// Pushes the closest point of vi to its neighbor in direction d.  Returns
// true if the closest point of the neighbor changed.
//------------------------------------------------------------------------------
static HeapVertex PushDistance(
    UVertexNetwork& vertices, const int vi, const Direction& d,
    const intn& v_point, const int label,
    const OctreeOptions& o) {

  // const intn curve_point = vertices.ClosestPoint(vi);
  // const int curve_pointi = vertices.ClosestPointIndex(vi);
  const intn curve_point = ClosestPoint(vi, vertices);
  const int curve_pointi = ClosestPointIndex(vi, vertices);

  // const int n_vi = vertices.Neighbor(vi, d);
  const int n_vi = Neighbor(vi, d, vertices);
  // const intn n_v_point = vertices.Position(n_vi);
  const intn n_v_point = vertices.vertices[n_vi].position;
  const int n_dist = length(n_v_point - curve_point);

  // const bool closest_changed = vertices.SetClosestPoint(n_vi, curve_pointi);
  const bool closest_changed = SetClosestPoint(n_vi, curve_pointi, &vertices);
  if (closest_changed) {
    return HeapVertex(n_dist, n_vi);
  }
  return HeapVertex(-1, -1);
}

//------------------------------------------------------------------------------
// SetDistances
//
// Performs a wavefront expansion out from the cells in the heap, setting
// distances on vertices as it goes.
//------------------------------------------------------------------------------
void SetDistances(
    VertexNetwork& vertices, std::multiset<HeapVertex>& heap,
    const OctreeOptions& o);

void SetDistances(
    UVertexNetwork& vertices, std::multiset<HeapVertex>& heap,
    const OctreeOptions& o);

//------------------------------------------------------------
// LoadWavefront
//
// Given a vertex vi, find all of its neighbors and add them
// to a wavefront.
//------------------------------------------------------------
template <int D>
void LoadWavefront(const int vi, const intn& center,
                   const VertexNetwork& vertices,
                   multiset<HeapVertex>& heap) {
  // typedef Direction<D> Direction_t;
  typedef Direction Direction_t;

  for (int axis = 0; axis < D; ++axis) {
    for (int pos = 0; pos < 2; ++pos) {
      // const Direction_t d = Direction_t::FromAxis(axis, pos);
      const Direction_t d = DirectionFromAxis(axis, pos);
      const int n_vi = vertices.Neighbor(vi, d);
      if (n_vi != -1 && vertices.ClosestPointIndex(n_vi) != -1) {
        const index_t level = vertices.NeighborLevel(vi, d);
        const index_t width = Level2CellWidth(level);
        const intn v_point = Position(d.pos, d.neg, center, width);
        const index_t dist = length(vertices.ClosestPoint(n_vi)-v_point);
        heap.insert(HeapVertex(dist, n_vi));
      }
    }
  }
}

template <int D>
class VertexVisitor {
 public:
  virtual bool operator()(
      const int vi, const intn& p) = 0;
};

// Ambiguous cells snipped

template <int D, typename LabeledGeometry>
struct NonEmptyVertexCollector : public VertexVisitor<D> {
  typedef VertexNetwork VertexNetwork_t;

  NonEmptyVertexCollector(
      VertexNetwork_t& vertices_,
      std::vector<int>& cell_vertices_,
      std::vector<intn >& points_,
      std::vector<std::vector<LabeledGeometry> >& cell_geometries_)
      : vertices(vertices_), cell_vertices(cell_vertices_),
      points(points_), cell_geometries(cell_geometries_) {}

  bool operator()(const int vi, const intn& p) {
    if (vertices.IsBase(vi) && !cell_geometries[vi].empty()) {
      cell_vertices.push_back(vi);
      points.push_back(p);
    }
    return true;
  }

  VertexNetwork_t& vertices;
  std::vector<int>& cell_vertices;
  std::vector<intn >& points;
  std::vector<std::vector<LabeledGeometry> >& cell_geometries;
};

struct UNonEmptyVertexCollector : public VertexVisitor<DIM> {
  UNonEmptyVertexCollector(
      UVertexNetwork& vertices_,
      std::vector<int>& cell_vertices_,
      std::vector<intn >& points_,
      std::vector<std::vector<LabeledGeometry> >& cell_geometries_)
      : vertices(vertices_), cell_vertices(cell_vertices_),
      points(points_), cell_geometries(cell_geometries_) {}

  bool operator()(const int vi, const intn& p) {
    if (IsBase(vi, vertices) && !cell_geometries[vi].empty()) {
      cell_vertices.push_back(vi);
      points.push_back(p);
    }
    return true;
  }

  UVertexNetwork& vertices;
  std::vector<int>& cell_vertices;
  std::vector<intn >& points;
  std::vector<std::vector<LabeledGeometry> >& cell_geometries;
};

template <int D>
struct CompareVerticesVisitor {
  typedef Direction Direction_t;
  typedef Direction Dir;

  CompareVerticesVisitor() : _points(new std::vector<intn>()) {}

  bool operator()(const int base_vi, const intn& p) {
    _points->push_back(p);
    return true;
  }

  bool operator==(const CompareVerticesVisitor& rhs) {
    if (_points->size() != rhs._points->size()) return false;
    for (int i = 0; i < _points->size(); ++i) {
      if (_points->at(i) != rhs._points->at(i)) return false;
    }
    return true;
  }

  bool operator!=(const CompareVerticesVisitor& rhs) {
    return !(*this == rhs);
  }

 private:
  shared_ptr<std::vector<intn> > _points;
};

//------------------------------------------------------------
// FindPositionsVisitor
// Finds all positions of vertices whose indices are given.
// This is efficient if the number of vertices is small.
//------------------------------------------------------------
template <int D>
struct FindPositionsVisitor : public VertexVisitor<D> {
  typedef std::map<int, intn> Map;
  typedef typename Map::iterator Iterator;

  template <typename Iter>
  FindPositionsVisitor(Iter begin, Iter end)
      : _positions(new Map()) {
    for (Iter it = begin; it != end; ++it) {
      (*_positions)[*it] = make_intn();
    }
  }

  static FindPositionsVisitor<D> ForOne(const int vi) {
    FindPositionsVisitor<D> fpv;
    fpv.Add(vi);
    return fpv;
  }

  static FindPositionsVisitor<D> ForAll(const int count) {
    FindPositionsVisitor<D> fpv;
    fpv._positions.reset(0);
    fpv._array.reset(new vector<intn>(count));
    return fpv;
  }

  bool operator()(const int base_vi, const intn& p) {
    if (_positions) {
      Iterator it = _positions->find(base_vi);
      if (it != _positions->end()) {
        it->second = p;
      }
    } else {
      Ensure(*_array.get(), base_vi);
      _array->at(base_vi) = p;
    }
    return true;
  }

  const intn& operator[](const int vi) const {
    if (_positions) {
      Iterator it = _positions->find(vi);
      if (it != _positions->end()) {
        return it->second;
      }
      throw logic_error("Unknown vi in FindPositionsVisitor");
    } else {
      return _array->at(vi);
    }
  }

  bool empty() const { return _positions->empty(); }
  void clear() { _positions->clear(); }

 private:
  FindPositionsVisitor() : _positions(new Map()) {}

  void Add(const int vi) {
    (*_positions)[vi] = make_intn();
  }

 private:
  shared_ptr<Map> _positions;
  shared_ptr<std::vector<intn> > _array;
};

//------------------------------------------------------------
// IndexAndPoint
//------------------------------------------------------------
template <int D>
struct IndexAndPoint {
  IndexAndPoint()
      : vi(-1) {}
  IndexAndPoint(const int vi_, const intn& v_point_)
      : vi(vi_), v_point(v_point_) {}
  int vi;
  intn v_point;
};

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
// Function declarations
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

// void InitWaveGpu_old(
//     std::vector<int>& nonempty_vertices,
//     std::vector<intn >& nonempty_points,
//     shared_array<level_t> all_levels,
//     std::vector<vector<LabeledGeometry> >& base2geometries,
//     GeomVertices& geom_vertices,
//     UVertexNetwork& vertices,
//     multiset<HeapVertex>& heap,
//     const OctreeOptions& o);

template <int D>
std::vector<std::vector<int> > ComputeBase2Incident(
    const VertexNetwork& vertices);

//------------------------------------------------------------------------------
// BuildOctree
//
// Given polylines or triangles, build an octree.
//------------------------------------------------------------------------------
// ManagedVertexNetwork BuildOctree(
VertexNetwork BuildOctree( // TODO - why are these named differently than in the cpp file?
    const std::vector<std::vector<floatn> >& all_vertices,
    const std::vector<std::vector<Face> >& all_faces,
    const BoundingBox<floatn>& bb,
    const OctreeOptions& o);
VertexNetwork BuildExtendedOctree(
    const std::vector<std::vector<floatn> >& all_vertices,
    const std::vector<std::vector<Face> >& all_faces,
    const BoundingBox<floatn>& bb,
    const OctreeOptions& o);
}
std::vector<oct::LabeledGeometry> GenerateLabeledGeometries(
    const vector<vector<floatn> >& label2gverts,
    const vector<vector<Face> >& label2faces,
    const BoundingBox<floatn>& bb);

class LabelPair {
 public:
  LabelPair() {}
  LabelPair(int a, int b) {
    labels[0] = a;
    labels[1] = b;
  }

  int operator[](const int i) const {
    return labels[i];
  }
  int& operator[](const int i) {
    return labels[i];
  }

  bool operator<(const LabelPair& rhs) const {
    for (int i = 0; i < 2; ++i) {
      if (labels[i] != rhs[i])
        return labels[i] < rhs[i];
    }
    return false;
  }
  bool operator==(const LabelPair& rhs) const {
    for (int i = 0; i < 2; ++i) {
      if (labels[i] != rhs[i])
        return false;
    }
    return true;
  }

 private:
  int labels[2];
};

#endif
