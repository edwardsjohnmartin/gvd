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

#ifndef __Segment2_H__
#define __Segment2_H__

#include <vector>

#include "../opencl/defs.h"
#include "../opencl/vec.h"

#include "../timer.h"
#include "../opencl/vertex.h"
#include "../vertex_network.h"
#include "./common.h"
#include "../bb.h"
#include "./gl2d.h"
#include "../octree.h"
#include "../graph.h"

class GVDViewer2 : public GL2D {
 public:
  typedef oct::Timer Timer;
  typedef oct::index_t index_t;
  // typedef oct::VertexNetwork<2> VertexNetwork;

  static const int kWidth = oct::kWidth;
  static const float3 red;

 public:
  GVDViewer2(const int win_width, const int win_height);

  void PrintCommands() const;

  void AddPoint(int x, int y);

  void ReadMesh(const std::string& filename);
  void WritePolygons() const;

  virtual void ProcessArgs(int argc, char** argv);
  virtual void Mouse(int button, int state, int x, int y);
  virtual void MouseMotion(int x, int y);
  virtual void PassiveMouseMotion(int x, int y);
  virtual void Keyboard(unsigned char key, int x, int y);
  virtual void Special(int key, int x, int y);
  virtual void Display();

  float2 Obj2Oct(const float2& v) const;
  float2 Oct2Obj(const int2& v) const;
  GLfloat Oct2Obj(int dist) const;

  unsigned char outcode(const float2& v) const;

  bool InBounds(const float2& v) const;

  void Zoom(const float2& target, const float zoom);

  // Given a code, returns the corner that represents the
  // intersection of the code's line and the next code's line.
  float2 corner(unsigned char code) const;

  // Clip a line to the bounds
  // a is outside, b is inside
  std::vector<float2> clip(const float2& a0, const float2& a1,
                     const float2& b0, const float2& b1) const;

  void MakePolygon();

  void SetStartSearch(int x, int y);
  void SetEndSearch(int x, int y);
  void Search(const int start, const int end);

  // void DefaultCallback(
  //     const int3& base_point, const index_t level,
  //     const index_t max_level, const bool complete);
  bool StoreVertexLocation(const int vi, const int2& p);
  void BuildOctree();
  void glSquare(GLfloat x, GLfloat y, GLfloat w) const;
  void glSquare(const float2& p, GLfloat w) const;
  void glSquareWinWidth(const float2& p, GLfloat w) const;
  bool DrawEdge(const int vi, const int n_vi,
                const int2& p, const int2& q,
                // const oct::Direction<2>& d) const;
                const oct::Direction& d) const;
  void DrawOctree() const;
  void DrawGVD() const;
  void DrawPath() const;
  void DrawVoronoi() const;
  bool DrawVertexID(const int vi, const int2& p) const;
  void DrawVertexIDs();
  bool DrawVertexLabel(const int vi, const int2& p) const;
  void DrawVertexLabels();
  bool DrawVertexDistanceLine(const int vi, const int2& p) const;
  void DrawVertexDistanceLines();
  bool DistanceFunctionVertex(const int vi, const int2& p) const;
  void DrawDistanceField();
  void PrintStatistics() const;
  void PrintHelp() const;

  void SetShowHelp(const bool show_help_) {
    show_help = show_help_;
  }
  void SetShowOctree(const bool show_octree_) {
    show_octree = show_octree_;
  }
  void SetShowGVD(const bool show_gvd_) {
    show_gvd = show_gvd_;
  }
  void SetShowVoronoi(const bool b) {
    show_voronoi = b;
  }
  void SetShowClosestPointLine(const bool show_closest_point_line_) {
    show_closest_point_line = show_closest_point_line_;
  }
  void SetShowVertexLabel(const bool show_vertex_labels_) {
    show_vertex_labels = show_vertex_labels_;
  }
  void SetShowVertexId(const bool show_vertex_ids_) {
    show_vertex_ids = show_vertex_ids_;
  }
  void SetShowClosestPointPath(const bool show_closest_point_path_) {
    show_closest_point_path = show_closest_point_path_;
  }
  void SetShowVertDist(const bool show_vertex_distance_) {
    show_vertex_distance = show_vertex_distance_;
  }
  void SetShowStatistics(const bool b) {
    show_statistics = b;
  }

  oct::VertexNetwork& Vertices() { return vertices; }
  const oct::VertexNetwork& Vertices() const { return vertices; }
  
 private:
  void HelpString(const std::string msg, const int i) const;

 private:
  float2 mouse_obj;
  bool mouse_active;

  std::vector<float2> verts;
  std::vector<std::vector<float2> > polygons;
  BoundingBox<float2> bb;

  oct::VertexNetwork vertices;
  bool dirty;
  Graph<2> gvd_graph;

  std::vector<float2> vertex_locations;
  std::vector<bool> vertex_circle;
  std::vector<float2> middle_down;
  std::vector<float2> middle_up;

  bool show_help;
  bool show_advanced_help;
  bool show_octree;
  bool show_vertex_labels;
  bool show_cell_vertices;
  bool show_vertex_distance;
  int polygon_mode;
  bool show_closest_point_line;
  bool show_closest_point_path;
  bool show_gvd;
  bool show_path;
  bool show_voronoi;
  bool show_vertex_ids;
  bool show_statistics;
  // 0 = none
  // 1 = GVD
  // 2 = distance field, sampled
  // 3 = distance field, interpolated
  int show_distance_field;
  double max_vertex_distance;
  float max_dist_oct;
  float max_dist_obj;
  oct::OctreeOptions o;
  // 0 = continuous, closed polygons
  // 1 = continuous, open
  // 2 = click polygons
  int entry_mode;
  std::vector<int> search_path;
  double3 octree_color;
};

#endif
