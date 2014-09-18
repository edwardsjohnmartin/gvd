#ifndef __JME_TILE3_H__
#define __JME_TILE3_H__

#include "./octree.h"
#include "./graph.h"
#include "./shared_ptr.h"
#include "./ambiguous.h"
#include "./medial.h"

struct ViNormal {
  ViNormal() {}
  ViNormal(const int vi_, const int normal_)
      : vi(vi_), normal(normal_) {}
  int vi;
  int normal;
};

// A point of intersection and two labels that are split by this intersection.
// The labels are in traversal order of the associated quadtree cell.
//
// l[0] _________________ 
//     |                 |
//     |                 |
//     |                 |
//     * p               |
//     |                 |
//     |                 |
//     |_________________|
// l[1]                   
struct LabeledIntersection {
  LabeledIntersection(
      const int3& p_, const int l1, const int l2, const int dist_)
      : p(p_), dist(dist_) {
    label[0] = l1;
    label[1] = l2;
  }
  int3 p;
  int label[2];
  int dist;
};

struct LabeledEdge {
  LabeledEdge(const Edge& e_, const int olabel_)
      : e(e_), olabel(olabel_) {}
  Edge e;
  int olabel;
};

struct TileSurfaceVisitorGpu2 {
  static const int D = 3;
  typedef oct::Direction Direction_t;
  typedef oct::FindBasesVisitor<D> FBV;
  static const int M = (1<<D);

  TileSurfaceVisitorGpu2(const oct::FindBasesVisitor<3>* fbv,
                           const oct::FindBasesVisitor<3>* fbv_2,
                           vector<vector<LabeledIntersection> >* p_cell2ints,
                           const oct::VertexNetwork* vertices,
                           const oct::OctreeOptions& o)
      : _fbv(fbv), _fbv_2(fbv_2),
    _p_cell2ints(p_cell2ints),
    _vertices(vertices), _o(o) {}

  bool operator()(const int vi, const int3& p);

  void GetSharedCells(
      const int vi, const int n_vi, const int axis, const int oaxis,
      const oct::FindBasesVisitor<3>* fbv,
      ViNormal* shared);

 private:
  const oct::FindBasesVisitor<3>* _fbv;
  const oct::FindBasesVisitor<3>* _fbv_2;
  // For a cell, stores all edge intersections incident to the cell.
  // Accessed as
  //   _cell2ints[normal_idx][base_vi]
  // where normal_idx \in {0,1,2}.
  vector<vector<LabeledIntersection> >* _p_cell2ints;
  const oct::VertexNetwork* _vertices;
  const oct::OctreeOptions _o;
};

struct TileSurfaceVisitorGpu3 {
  static const int D = 3;
  typedef oct::Direction Direction_t;
  typedef oct::FindBasesVisitor<D> FBV;
  static const int M = (1<<D);

  TileSurfaceVisitorGpu3(
      const oct::FindBasesVisitor<3>* fbv,
      const vector<vector<LabeledIntersection> >* p_cell2ints,
      vector<vector<LabeledEdge> >* cell2edges,
      vector<intn>* tri_verts,
      vector<int>* vert_dist,
      TopoGraph* gvd_graph,
      const int num_meshes,
      const oct::OctreeOptions& o)
      : _fbv(fbv),
        _p_cell2ints(p_cell2ints),
        _cell2edges(cell2edges),
        _tri_verts(tri_verts),
        _vert_dist(vert_dist),
        _gvd_graph(gvd_graph),
        _num_meshes(num_meshes),
        _o(o) {}

  bool operator()(const int p_base_vi, const int3& p);

 private:
  const oct::FindBasesVisitor<3>* _fbv;
  // For a cell, stores all edge intersections incident to the cell.
  // Accessed as _cell2ints[normal_axis][base_vi] where normal_axis \in {0,1,2}.
  const vector<vector<LabeledIntersection> >* _p_cell2ints;
  // Maps a 3D cell to all triangle edges on its boundary, one array per label
  vector<vector<LabeledEdge> >* _cell2edges;
  // Vertices of the triangles, one array per label
  vector<intn>* _tri_verts;
  // Distance each vertex is from a surface.  Corresponds to _tri_verts.
  vector<int>* _vert_dist;
  TopoGraph* _gvd_graph;
  int _num_meshes;
  const oct::OctreeOptions _o;
};

#endif
