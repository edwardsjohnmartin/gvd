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

#ifndef __SEARCH_H__
#define __SEARCH_H__

#include "./opencl/defs.h"
#include "./opencl/bit.h"
#include "./opencl/vec.h"
#include "./opencl/triangle.h"
#include "./opencl/vertex.h"

#include "./timer.h"
#include "./shared_ptr.h"
#include "./bb.h"
#include "./opencl.h"
#include "./vectorn.h"

namespace oct {

// This struct is to be used with BFS and is optimized for
// small searches.
struct VisitedSet {
  VisitedSet() {}
  VisitedSet(const size_t size) {}
  bool operator[](const int vi) const {
    return visited.find(vi) != visited.end();
  }
  void SetVisited(const int vi) {
    visited.insert(vi);
  }
  std::set<int> visited;
};

// This struct is to be used with BFS and is optimized for
// full octree searches.
struct VisitedVector {
  VisitedVector() {}
  VisitedVector(const size_t size)
      : visited(size) {}
  bool operator[](const int vi) const {
    return visited[vi];
  }
  void SetVisited(const int vi) {
    visited[vi] = true;
  }
  std::vector<bool> visited;
};

template <int D>
class CompositeVertexVisitor {
 public:
  void Add(VertexVisitor<D>* v) {
    _visitors.push_back(v);
  }

  bool operator()(const int vi, const intn& p) {
    for (int i = 0; i < _visitors.size(); ++i) {
      if (!(*_visitors[i])(vi, p)) return false;
    }
    return true;
  }
 private:
  std::vector<VertexVisitor<D>*> _visitors;
};

//------------------------------------------------------------------------------
// BFS
//------------------------------------------------------------------------------

// Default BFS visitors
template <int D>
bool DefaultVertexVisitor(const int vi, const intn& p) {
  return true;
}

template <int D>
bool DefaultEdgeVisitor(const int vi, const int n_vi,
                        const intn& p,
                        const intn& n_p,
                        // const oct::Direction<D>& d) {
                        const oct::Direction& d) {
  return true;
}

//------------------------------------------------------------------------------
// BFS
//
// Performs a breadth-first search over the vertices from the (0, 0, 0) cell
// VertexVisitor is of the form
//    bool Visit(const int vi, const int3& p);
// EdgeVisitor is of the form
//    bool Visit(const int vi, const int n_vi,
//               const int3& p, const int3& q, const Direction<D>& d);
//------------------------------------------------------------------------------
template <int D, typename VisitedStruct,
typename VertexVisitor, typename EdgeVisitor>
void BFS(const VertexNetwork& vertices,
         const intn& base_point,
         const std::vector<int>& axes,
         const int max_length, const int start_vi,
         VertexVisitor v, EdgeVisitor ev) {
  if (vertices.empty()) return;

  const int num_axes = axes.size();
  VisitedStruct visited(vertices.size());
  std::list<std::pair<int, intn> > queue;
  queue.push_back(make_pair(start_vi, base_point));
  while (!queue.empty()) {
    const int vi = queue.front().first;
    const intn p = queue.front().second;
    queue.pop_front();
    bool visit = !visited[vi];
    for (int i = 0; visit && i < D; ++i) {
      visit = (p.s[i] - base_point.s[i]) <= max_length;
    }
    if (visit) {
      visited.SetVisited(vi);
      if (vertices.Position(vi) != p) {
        cerr << "Not equal points: " << endl;
        cerr << "  vi = " << vi << endl;
        cerr << "  p = " << p << endl;
        cerr << "  Position = " << vertices.Position(vi) << endl;
        throw logic_error("Not equal points");
      }
      if (!v(vi, p)) break;

      for (int i = 0; i < num_axes; ++i) {
        const int axis = axes[i];
        // const oct::Direction<D> d = Direction::FromAxis(axis, true);
        const oct::Direction d = DirectionFromAxis(axis, true);
        const int n_vi = vertices.Neighbor(vi, d);
        if (n_vi != -1) {
          intn q = p;
          q.s[axis] += vertices.NeighborDist(vi, d);
          if (q.s[axis] - base_point.s[axis] <= max_length) {
            queue.push_back(make_pair(n_vi, q));
            if (!ev(vi, n_vi, p, q, d)) return;
          }
        }
      }
    }
  }
  return;
}

//------------------------------------------------------------
// BFS
// The following visit functions are optimized for full octree
// searches.
//------------------------------------------------------------

template <int D, typename Visitor, typename EdgeVisitor>
void BFS(const VertexNetwork& vertices, Visitor v, EdgeVisitor ev) {
  std::vector<int> axes;
  for (int i = 0; i < D; ++i) axes.push_back(i);
  BFS<D, VisitedVector>(vertices, make_intn(),
                        axes, kWidth, 0, v, ev);
}

template <int D, typename Visitor>
void VisitVerticesBFS(const VertexNetwork& vertices, Visitor v) {
  BFS<D>(vertices, v, DefaultEdgeVisitor<D>);
}

template <int D, typename Visitor>
void VisitVertices(const VertexNetwork& vertices, Visitor v) {
  // BFS<D>(vertices, v, DefaultEdgeVisitor<D>);
  VisitVerticesBFS<D>(vertices, v);
}

template <int D, typename Visitor>
void VisitEdgesBFS(const VertexNetwork& vertices, Visitor v) {
  BFS<D>(vertices, DefaultVertexVisitor<D>, v);
}

template <int D, typename Visitor>
void VisitEdges(const VertexNetwork& vertices, Visitor v) {
  // BFS<D>(vertices, DefaultVertexVisitor<D>, v);
  VisitEdgesBFS<D>(vertices, v);
}

//------------------------------------------------------------
// BFS
// The following visit functions are optimized for small
// searches.
//------------------------------------------------------------

template <int D, typename Visitor>
void VisitVertices(const VertexNetwork& vertices,
                   const intn& base_point,
                   const std::vector<int>& axes,
                   const int max_length, const int start_vi,
                   Visitor v) {
  BFS<D, VisitedSet>(vertices, base_point,
                     axes, max_length, start_vi, v, DefaultEdgeVisitor<D>);
}

template <int D, typename Visitor>
void VisitVertices(const VertexNetwork& vertices,
                   const intn& base_point,
                   const int max_length, const int start_vi,
                   Visitor v) {
  std::vector<int> axes;
  for (int i = 0; i < D; ++i) axes.push_back(i);
  BFS<D, VisitedSet>(vertices, base_point,
                     axes, max_length, start_vi, v, DefaultEdgeVisitor<D>);
}

template <int D, typename Visitor>
void VisitEdges(const VertexNetwork& vertices,
                const intn& base_point,
                const std::vector<int>& axes,
                const int max_length, const int start_vi,
                Visitor v) {
  BFS<D, VisitedSet>(vertices, base_point, axes, max_length, start_vi,
                     DefaultVertexVisitor<D>, v);
}

template <int D, typename Visitor>
void VisitEdges(const VertexNetwork& vertices,
                const intn& base_point,
                const int max_length, const int start_vi,
                Visitor v) {
  std::vector<int> axes;
  for (int i = 0; i < D; ++i) axes.push_back(i);
  BFS<D, VisitedSet>(vertices, base_point, axes, max_length, start_vi,
                     DefaultVertexVisitor<D>, v);
}

template <int D, typename Visitor>
Visitor WalkAroundFace(const int vi,
                       const intn& base_point,
                       const level_t level,
                       const int primary_axis, const int secondary_axis,
                       const VertexNetwork& vertices, Visitor v) {
  const int face_width = Level2CellWidth(level);

  int cur_vi = vi;
  intn cur_point = base_point;
  Direction dirs[4] = {
    // Direction::FromAxis(primary_axis, true),
    // Direction::FromAxis(secondary_axis, true),
    // Direction::FromAxis(primary_axis, false),
    // Direction::FromAxis(secondary_axis, false) };
    DirectionFromAxis(primary_axis, true),
    DirectionFromAxis(secondary_axis, true),
    DirectionFromAxis(primary_axis, false),
    DirectionFromAxis(secondary_axis, false) };
  for (int i = 0; i < 4; ++i) {
    const Direction d = dirs[i];
    index_t length = 0;
    do {
      const int n_vi = vertices.Neighbor(cur_vi, d);
      if (n_vi == -1) {
        cout << "Neighbor is -1. Original vertex is: " << vi << endl;
        cout << "Prime axis is:" << primary_axis <<
            " Sec Axis: " << secondary_axis << endl;
        return v;
      }
      const int sub_length = Level2CellWidth(vertices.NeighborLevel(cur_vi, d));
      length += sub_length;
      const intn n_point = Position(d.pos, d.neg, cur_point, sub_length);
      if (!v(cur_vi, n_vi, cur_point, n_point, d)) return v;
      cur_vi = n_vi;
      cur_point = n_point;
    } while (length < face_width);
  }
  return v;
}

}

#endif
