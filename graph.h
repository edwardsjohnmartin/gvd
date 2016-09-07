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

#ifndef __JME_GRAPH_H__
#define __JME_GRAPH_H__

#include <vector>
#include <limits>

#include "./opencl/vec.h"
#include "./opencl/edge.h"
#include "./edge_cpp.h"

// Only topological information.  No spatial information.
class TopoGraph {
 public:
  typedef std::pair<int, int> Pair;
  typedef std::set<int> Set;

 public:
  TopoGraph() {}

  void AddEdge(int ai, int bi) {
    Ensure(ai);
    Ensure(bi);
    _links[ai].insert(bi);
    _links[bi].insert(ai);
  }

  // a is pushing b
  double3 Dir(int ai, int bi) const {
    if (ai == bi) return double3();

    std::pair<int, int> ab(ai, bi);
    return _link2dir.find(ab)->second;

    // const int ri = WhichRing(ai, bi);
    // Set aring = Ring(ai, ri);
    // aring.insert(Ring(ai, ri-1).begin(), Ring(ai, ri-1).end());
    // const Set& bring = Ring(bi, 1);

    // // // Find the direction between ai and bi and also between every
    // // // ni where ni is a neighbor of ai.
    // // std::pair<int, int> ab(ai, bi);
    // // double3 dir = _link2dir.find(ab)->second;
    // double3 dir;
    // // const std::set<int>& a_nbrs = _links[ai];
    // // const std::set<int>& b_nbrs = _links[bi];
    // typedef std::list<int> List;
    // typedef List::const_iterator iter;

    // // Neighbors of both a and b
    // List ab_nbrs;
    // // std::set_intersection(a_nbrs.begin(), a_nbrs.end(),
    // //                       b_nbrs.begin(), b_nbrs.end(),
    // //                       std::inserter(ab_nbrs, ab_nbrs.end()));
    // std::set_intersection(aring.begin(), aring.end(),
    //                       bring.begin(), bring.end(),
    //                       std::inserter(ab_nbrs, ab_nbrs.end()));
    // for (iter it = ab_nbrs.begin(); it != ab_nbrs.end(); ++it) {
    //   const int ni = *it;
    //   std::pair<int, int> nb(ni, bi);
    //   dir += _link2dir.find(nb)->second;
    // }
    // return dir.unit();
  }

  void SetDir(int ai, int bi, const double3& dir) {
    Pair ab(ai, bi);
    _link2dir[ab] = dir;
  }

  void SetArea(int ai, int bi, double area) {
    Pair ab(ai, bi);
    _link2area[ab] = area;
  }

  double GetArea(int ai, int bi) const {
    Pair ab(ai, bi);
    return _link2area.find(ab)->second;
  }

  const std::set<int>& operator[](int i) const { return _links[i]; }
  const std::set<int>& at(int i) const { if (i >= _links.size()) return _empty_set;
    return _links[i]; }

  const std::set<int>& Ring(int i, int ring) const {
    if (_rings.empty())
      ComputeRings();
    while (_rings[i].size() <= ring) {
      _rings[i].push_back(Set());
    }
    return _rings[i][ring];
  }

  int WhichRing(int center, int i) const {
    for (int ri = 0; ri < 10000; ++ri) {
      const Set& ring = Ring(center, ri);
      if (ring.find(i) != ring.end())
        return ri;
    }
    throw logic_error("Not found in ring");
  }

 private:
  void Ensure(int i) {
    if (_links.size() <= i) {
      _links.resize(i+1);
    }
  }

  void ComputeRings() const {
    _rings.resize(_links.size());
    for (int i = 0; i < _links.size(); ++i) {
      _rings[i].push_back(Set());
      _rings[i].push_back(Set());
      _rings[i][0].insert(i);
      const std::set<int>& nbrs = _links[i];
      _rings[i][1] = nbrs;
    }
    
    // ring index
    int ri = 2;
    bool stop = false;
    while (!stop) {
      stop = true;
      for (int i = 0; i < _links.size(); ++i) {
        const Set& prev_ring = _rings[i][ri-1];
        Set used = prev_ring;
        used.insert(_rings[i][ri-2].begin(), _rings[i][ri-2].end());
        Set ring;
        for (Set::const_iterator it = prev_ring.begin();
             it != prev_ring.end(); ++it) {
          const int j = *it;
          const Set& j_ring = _rings[j][1];
          std::set_difference(j_ring.begin(), j_ring.end(),
                              used.begin(), used.end(),
                              std::inserter(ring, ring.end()));
        }
        _rings[i].push_back(Set());
        _rings[i][ri] = ring;
        if (!ring.empty())
          stop = false;
      }
      ++ri;
    }
  }

 private:
  std::vector<std::set<int> > _links;
  // Maps edges of the graph to directions
  std::map<std::pair<int, int>, double3> _link2dir;
  std::map<std::pair<int, int>, double> _link2area;
  // _rings[0][1] refers to the zeroth object's one-ring neighbors
  // _rings[3][2] refers to the third object's two-ring neighbors
  mutable std::vector<std::vector<std::set<int> > > _rings;

  std::set<int> _empty_set;
};

struct CostStruct {
  CostStruct(int vi_, double cost_)
      : vi(vi_), cost(cost_) {}
  bool operator<(const CostStruct& rhs) const {
    if (cost == rhs.cost)
      return vi < rhs.vi;
    return cost < rhs.cost;
  }
  int vi;
  double cost;
};

template <int D>
class Graph {
 public:
  Graph() {}

  int Add(const doublen& p) {
    _vertices.push_back(p);
    _links.push_back(std::vector<int>());
    return _vertices.size()-1;
  }

  void AddEdge(int ai, int bi) {
    _links[ai].push_back(bi);
    _links[bi].push_back(ai);
  }

  template <typename Visitor>
  void VisitEdges(Visitor v) const {
    for (int vi = 0; vi < _vertices.size(); ++vi) {
      for (int i = 0; i < _links[vi].size(); ++i) {
        const int n_vi = _links[vi][i];
        v(vi, _vertices[vi], n_vi, _vertices[n_vi]);
      }
    }
  }

  // Gives the path from the end to the front
  template <typename Out_iter>
  void Dijkstra(const int start, const int end, Out_iter path) {
    std::vector<int> parents(_vertices.size(), -1);
    std::set<CostStruct> heap;
    std::vector<double> cost(_vertices.size(),
                             std::numeric_limits<double>::max());
    std::vector<bool> visited(_vertices.size(), false);
    heap.insert(CostStruct(start, 0));
    cost[start] = 0;
    bool done = false;
    while (!heap.empty() && !done) {
      const int vi = heap.begin()->vi;
      heap.erase(heap.begin());
      if (vi == end) {
        done = true;
      } else if (!visited[vi]) {
        visited[vi] = true;
        for (int i = 0; i < _links[vi].size(); ++i) {
          const int n_vi = _links[vi][i];
          if (!visited[n_vi]) {
            const double d = cost[vi] + length(_vertices[vi]-_vertices[n_vi]);
            if (d < cost[n_vi]) {
              cost[n_vi] = d;
              heap.insert(CostStruct(n_vi, d));
              parents[n_vi] = vi;
            }
          }
        }
      }
    }

    if (!done) {
      // throw std::logic_error("Shortest cost path failed");
      std::cerr << "Shortest cost path failed" << std::endl;
      return;
    }

    int vi = end;
    while (vi != start) {
      *path++ = vi;
      vi = parents[vi];
    }
    *path++ = start;
  }

  const std::vector<doublen>& GetVertices() const { return _vertices; }
  std::vector<doublen>& GetVertices() { return _vertices; }
  const doublen& operator[](int i) const { return _vertices[i]; }
  doublen& operator[](int i) { return _vertices[i]; }

  void Clear() { _vertices.clear(); _links.clear(); }

 private:
  std::vector<doublen> _vertices;
  std::vector<std::vector<int> > _links;
};

template <int D>
class GraphConstructor {
 public:
  typedef typename std::map<Edge, int>::const_iterator e2i_iter;
  typedef typename std::map<doublen, int>::const_iterator c2i_iter;
  typedef typename std::map<doublen, int>::const_iterator p2i_iter;

 public:
  GraphConstructor() {}

  //------------------------------------------------------------
  // For use with 2D and 3D bisectors
  //------------------------------------------------------------
  void AddEdge(const doublen& a, const doublen& b, const Edge& a_seg) {
    int a_idx, b_idx;

    e2i_iter e2i = _edge2idx.find(a_seg);
    if (e2i != _edge2idx.end()) {
      a_idx = e2i->second;
    } else {
      a_idx = _graph.Add(a);
      _edge2idx[a_seg] = a_idx;
    }

    c2i_iter c2i = _centroid2idx.find(b);
    if (c2i != _centroid2idx.end()) {
      b_idx = c2i->second;
    } else {
      b_idx = _graph.Add(b);
      _centroid2idx[b] = b_idx;
    }
    
    _graph.AddEdge(a_idx, b_idx);
  }

  //------------------------------------------------------------
  // For use with 3D bisectors
  //------------------------------------------------------------
  void AddEdge(const doublen& a, const doublen& b) {
    int a_idx, b_idx;

    {
      p2i_iter p2i = _point2idx.find(a);
      if (p2i != _point2idx.end()) {
        a_idx = p2i->second;
      } else {
        a_idx = _graph.Add(a);
        _point2idx[a] = a_idx;
      }
    } {
      p2i_iter p2i = _point2idx.find(b);
      if (p2i != _point2idx.end()) {
        b_idx = p2i->second;
      } else {
        b_idx = _graph.Add(b);
        _point2idx[b] = b_idx;
      }
    }

    _graph.AddEdge(a_idx, b_idx);
  }

  //------------------------------------------------------------
  // For use with 3D bisectors
  //------------------------------------------------------------
  void AddTriangle(const doublen& a, const doublen& b, const doublen& c) {
    AddEdge(a, b);
    AddEdge(b, c);
    AddEdge(a, c);
  }

  void ReplacePoint(const doublen& oldp, const doublen& newp) {
    const int i = _point2idx[oldp];
    _graph[i] = newp;
    _point2idx[newp] = i;
  }

  const Graph<D>& GetGraph() const { return _graph; }

 private:
  Graph<D> _graph;
  std::map<Edge, int> _edge2idx;
  std::map<doublen, int> _centroid2idx;
  std::map<doublen, int> _point2idx;
};

typedef Graph<2> Graph2;
typedef Graph<3> Graph3;
typedef GraphConstructor<2> GraphConstructor2;
typedef GraphConstructor<3> GraphConstructor3;

#endif
