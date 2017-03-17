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

#ifndef __VERTEX_NETWORK_H__
#define __VERTEX_NETWORK_H__

#include <algorithm>
#include <vector>
#include <set>

#include "./opencl/vec.h"
#include "./opencl/vertex.h"
#include "./opencl/geometry.h"
#include "./shared_array.h"
#include "./shared_ptr.h"

namespace oct {

struct NullVisitor {
  void operator()(const int n_vi) {}
};

//----------------------------------------
// class VerticesChangeTracker
//----------------------------------------
class VerticesChangeTracker {
 public:
  VerticesChangeTracker()
      : _num_vertices_added(0), _num_closest_points_added(0) {}

  void VertexChanged(const int vi) {
    _vertices_changed.insert(vi);
  }
  void VertexAdded() {
    ++_num_vertices_added;
  }
  void ClosestPointAdded() {
    ++_num_closest_points_added;
  }

  bool Changed() const {
    return (!_vertices_changed.empty() ||
        _num_vertices_added > 0 ||
        _num_closest_points_added > 0);
  }

  // This is linear!
  int NumVerticesChanged() const {
    return _vertices_changed.size();
  }

  template <typename Out_iter>
  void VerticesChanged(Out_iter result) const {
    copy(_vertices_changed.begin(), _vertices_changed.end(), result);
  }

  int NumVerticesAdded() const {
    return _num_vertices_added;
  }
  
  int NumClosestPointsAdded() const {
    return _num_closest_points_added;
  }

  void Reset() {
    _vertices_changed.clear();
    _num_vertices_added = 0;
    _num_closest_points_added = 0;
  }
  
 private:
  std::set<int> _vertices_changed;
  int _num_vertices_added;
  int _num_closest_points_added;
};

typedef shared_ptr<VerticesChangeTracker> VerticesChangeTracker_h;

//----------------------------------------
// struct ManagedVertexNetwork
//----------------------------------------
struct ManagedVertexNetwork {
  int num_vertices;
  int vertex_array_capacity;
  shared_array<Vertex> _vertices;

  int num_cpoints;
  int cpoint_array_capacity;
  shared_array<GeomPoint> _closest_points;
};

//----------------------------------------
// class VertexNetwork
//----------------------------------------
class VertexNetwork {
 public:
  static const int kNoLabel = -1;

 public:
  VertexNetwork();
  VertexNetwork(int num_vertices, Vertex* vertices,
                int num_cpoints, GeomPoint* cpoints);
  ~VertexNetwork();

  //------------------------------------------------------------
  // Modifiers
  // Every modifier should update the change tracker _changes.
  //------------------------------------------------------------

  void Clear();

  // Constructs 2^D new cells by subdividing the cell base_vi.
  //
  // Returns the new subcells with lvi indexing (see bit.h)
  // slvi2vi is filled with full subdivided cell indices if non-null.
  shared_array<int> Subdivide(const int base_vi, int* slvi2vi = 0);

  bool SetClosestPoint(const int vi, const int candidate) {
    // if (_changes)
    //   _changes->VertexChanged(vi);

    bool changed = false;
    const int cur = v_closest_point(_vertices[vi]);
    if (cur == -1) {
      // No closest point
      v_set_closest_point(candidate, &_vertices[vi]);
      changed = true;
    } else {
      const intn& vp = Position(vi);
      const double cur_dist = length(_closest_points[cur].p-vp);
      const double new_dist = length(_closest_points[candidate].p-vp);
      if (new_dist < cur_dist) {
        v_set_closest_point(candidate, &_vertices[vi]);
        changed = true;
      }
    }
    return changed;
  }

  // Returns the new index
  int CreateNewClosestPoint(const intn& p, const int label) {
    _closest_points.push_back(make_gpoint(p, label));
    return _closest_points.size() - 1;
  }

  //------------------------------------------------------------
  // Miscellaneous functions
  //------------------------------------------------------------

  size_t size() const { return _vertices.size(); }
  bool empty() const { return _vertices.empty(); }

  // Return a copy.  Returning a reference is dangerous since
  // the vector could be reallocated upon resize.
  Vertex GetVertex(const int vi) const {
    return _vertices[vi];
  }

  Vertex operator[](int i) const {
    return _vertices[i];
  }

  int Neighbor(const int vi, const Direction& d) const {
    const Vertex& v = _vertices[vi];
    return v_neighbor_index(d, v);
  }

  bool HasNeighbor(const int vi, const Direction& d) const {
    return (Neighbor(vi, d) != -1);
  }

  bool IsNeighbor(const int vi, const int n_vi) const {
    for (int i = 1; i < (1<<DIM); i=(i<<1)) {
      Direction d = DirectionFromPosNeg(i, 0);
      int n = Neighbor(vi, d);
      if (n == n_vi) return true;
      d = DirectionFromPosNeg(0, i);
      n = Neighbor(vi, d);
      if (n == n_vi) return true;
    }
    return false;
  }

  level_t NeighborLevel(
      const int vi, const Direction& d) const {
    const Vertex& v = _vertices[vi];
    return v_neighbor_level(d, v);
  }

  // width is the total octree width
  int NeighborDist(const int vi, const Direction& d) const {
    const Vertex& v = _vertices[vi];
    // return oct::CellWidth(v_neighbor_level(d, v));
    return Level2CellWidth(v_neighbor_level(d, v));
  }

  // Finds the distance to the nearest neighbor
  int Alpha(const int vi) const {
    int alpha = kWidth;
    for (int i = 0; i < DIM; ++i) {
      for (int pos = 0; pos < 2; ++pos) {
        const Direction d = DirectionFromAxis(i, pos);
        const int nbr = Neighbor(vi, d);
        if (nbr != -1) {
          alpha = std::min(alpha, NeighborDist(vi, d));
        }
      }
    }
    return alpha;
  }

  template <typename Visitor>
  int FindNeighbor(const int a_vi, const Direction& d,
                   const level_t level, Visitor& v) const {
    // const int goal_dist = oct::CellWidth(level);
    const int goal_dist = Level2CellWidth(level);
    int cur_dist = 0;
    int cur_vi = a_vi;
    while (cur_dist < goal_dist) {
      if (v_neighbor_level(d, _vertices[cur_vi]) == -1) return -1;
      // cur_dist += oct::CellWidth(v_neighbor_level(d, _vertices[cur_vi]));
      cur_dist += Level2CellWidth(v_neighbor_level(d, _vertices[cur_vi]));
      cur_vi = v_neighbor_index(d, _vertices[cur_vi]);
      v(cur_vi);
    }
    if (cur_dist > goal_dist) return -1;
    return cur_vi;
  }

  int FindNeighbor(const int a_vi, const Direction& d,
                   const level_t level) const {
    NullVisitor v;
    return FindNeighbor(a_vi, d, level, v);
  }

  // Starting at vi (inclusive if d is positive/exclusive otherwise),
  // finds the base vertex in direction d.
  // Returns -1 if base doesn't exists or does not lie
  // in that direction
  int FindBase(const int vi, const Direction& d,
               int& dist_traveled) const {
    if (d.pos > 0 && IsBase(vi)) {
      return vi;
    }
    int cur_vi = vi;
    dist_traveled = 0;
    do {
      const int n_vi = Neighbor(cur_vi, d);
      if (n_vi == -1) {
        return n_vi;
      }
      dist_traveled += NeighborDist(cur_vi, d);
      cur_vi = n_vi;
    } while (!IsBase(cur_vi));
    // if (dist_traveled > oct::CellWidth(CellLevel(cur_vi))) {
    if (dist_traveled > Level2CellWidth(CellLevel(cur_vi))) {
      return -1;
    }
    return cur_vi;
  }

  int FindBase(const int vi, const Direction& d) const {
    int dist_traveled = 0;
    return FindBase(vi, d, dist_traveled);
  }

  int Label(const int vi) const {
    const int p_idx = v_closest_point(_vertices[vi]);
    if (p_idx == -1) return -1;
    return _closest_points[p_idx].l;
  }

  int LabelAtIdx(const int i) const {
    return _closest_points[i].l;
  }

  //------------------------------------------------------------
  // Closest point functions
  //------------------------------------------------------------

  // Return a copy.  Returning a reference is dangerous since
  // the vector could be reallocated upon resize.
  intn ClosestPoint(const int vi) const {
    const int p_idx = v_closest_point(_vertices[vi]);
    assert(p_idx != -1);
    return _closest_points[p_idx].p;
  }

  int ClosestPointIndex(const int vi) const {
    return v_closest_point(_vertices[vi]);
  }

  // Returns the distance from the vertex to its closest geometry point.
  // Returns -1 if there is no closest geometry point.
  int Dist(const int vi) const {
    const int p_idx = v_closest_point(_vertices[vi]);
    if (p_idx == -1) return -1;
    const intn a = _closest_points[p_idx].p;
    const intn b = Position(vi);
    return length(a-b);
  }

  //------------------------------------------------------------
  // Miscellaneous functions
  //------------------------------------------------------------

  // In 2D, a corner and base are equivalent.
  // In 3D, a base is always a corner, but a corner is not always a base.
  // normal - normal axis \in { 001, 010, 100 }
  bool IsCorner(const int vi, const int normal) const {
    for (int axis = 1; axis < (1<<DIM); axis=(axis<<1)) {
      if (axis != normal && Neighbor(vi, DirectionFromPosNeg(axis, 0)) == -1)
        return false;
    }
    return true;
  }

  bool IsBase(const int vi) const {
    return v_is_base(_vertices[vi]);
  }

  // These only apply if the vertex is the base of a leaf cell
  level_t CellLevel(const int vi) const {
    return v_cell_level(_vertices[vi]);
  }

  index_t CellWidth(const int vi) const {
    // return oct::CellWidth(CellLevel(vi));
    return Level2CellWidth(CellLevel(vi));
  }

  const int* GetCorners(const int base_vi) const {
    return _vertices[base_vi].corners;
  }

  const intn Position(const int vi) const {
    return _vertices[vi].position;
  }

  intn CellCenter(const int vi) const {
    return _vertices[vi].position + CellWidth(vi) / 2;
  }

  // Returns true if cells vi0 and vi1 are incident to each other
  bool AreAdjacent(const int vi0, const int vi1) const {
    if (!IsBase(vi0) || !IsBase(vi1))
      throw std::logic_error("IsIncident - not a base");
    const index_t w0 = CellWidth(vi0);
    const index_t w1 = CellWidth(vi1);
    for (int i = 0; i < DIM; ++i) {
      const int min0 = Position(vi0).s[i];
      const int min1 = Position(vi1).s[i];
      if (min0 < min1)
        if (min0 + w0 < min1)
          return false;
      if (min1 < min0)
        if (min1 + w1 < min0)
          return false;
    }
    return true;
  }

  // Returns true if vi is incident to the cell described by base_vi.
  // Check is done by looking at the positions.
  bool IsIncident(const int vi, const int base_vi) const {
    if (!IsBase(base_vi))
      throw std::logic_error("IsIncident - not a base");
    const index_t width = CellWidth(base_vi);
    const intn& p = Position(vi);
    const intn& base_p = Position(base_vi);
    for (int i = 0; i < DIM; ++i) {
      if (p.s[i] < base_p.s[i]) return false;
      if (p.s[i] > base_p.s[i]+width) return false;
    }
    return true;
  }

  // type: Vertex
  template <typename Out_iter>
  void CopyVertices(Out_iter result) const {
    // std::copy(_vertices.get(), _vertices.get() + num_vertices, result);
    std::copy(_vertices.begin(), _vertices.end(), result);
  }

  int NumClosestPoints() const {
    // return num_cpoints;
    return _closest_points.size();
  }

  // This function is linear!
  int NumCells() const {
    int count = 0;
    for (int i = 0; i < size(); ++i) {
      if (IsBase(i)) ++count;
    }
    return count;
  }

 private:
  //------------------------------------------------------------
  // private modifiers
  //------------------------------------------------------------
  // Not copy-constructable.  This is critical for various reasons,
  // including change tracking (see vertices_gpu_state.h).
  // VertexNetwork(const VertexNetwork& copy) {}
  // void operator=(const VertexNetwork& copy) {}

  void Initialize();

  // Returns the index of the vertex
  int CreateVertex(const intn& position);

  void AddEdge(const int a_vi, const int b_vi, const Direction& d,
               const level_t level);

  int Break(const int a_vi, const int b_vi, const Direction& d,
            const level_t level);

  void Subdivide(
      int slvi2vi[], const int offset, const int axes[], const int level,
      const int dim);

  void SetIsBase(const int vi, const bool leaf) {
    v_set_is_base(leaf, &_vertices[vi]);
  }

  void SetCellLevel(const int vi, level_t level) {
    v_set_cell_level(level, &_vertices[vi]);
  }

 private:
  //------------------------------------------------------------
  // Data members
  //------------------------------------------------------------

  std::vector<Vertex> _vertices;
  std::vector<GeomPoint> _closest_points;

  // mutable VerticesChangeTracker_h _changes;
};

}

#endif
