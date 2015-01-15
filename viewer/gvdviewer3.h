#ifndef __GVDVIEWER3_H__
#define __GVDVIEWER3_H__

#include <vector>

#include "../opencl/defs.h"
#include "../opencl/vec.h"
#include "../opencl/vertex.h"

#include "../timer.h"
#include "./common.h"
#include "../bb.h"
#include "./gl3d.h"
#include "./mesh.h"
#include "../octree.h"
#include "../graph.h"
#include "../statistics.h"

// TODO: Clean up textures
class GVDViewer3 : public GL3D {
 public:
  typedef oct::Timer Timer;
  typedef oct::index_t index_t;
  typedef oct::level_t level_t;
  typedef oct::VertexNetwork VertexNetwork;
  typedef oct::ManagedVertexNetwork ManagedVertexNetwork;

  static const int kWidth = oct::kWidth;
  static const float3 red;

 public:
  GVDViewer3(const int win_width, const int win_height);//,
          // const float2& world_min, const float2& world_max);
  ~GVDViewer3();

  int3 Obj2Oct(const float3& v) const;
  float3 Oct2Obj(const int3& v) const;
  GLfloat Oct2Obj(int dist) const;
  bool DrawVertexDistanceLine(const int vi, const int3& p) const;
  void DrawVertexDistanceLines();
  bool DrawLabel(const int vi, const int3& p) const;
  void DrawVertexLabels();
  bool DrawID(const int vi, const int3& p) const;
  void DrawVertexIDs();
  bool DrawSeparatorVertex(const int vi, const int n_vi,
                         const int3& p, const int3& q,
                           // const oct::Direction<3>& d) const;
                           const oct::Direction& d) const;

  void DrawMesh(const Mesh& mesh, const bool face_normals = false);
  void DrawMesh(const Mesh& mesh, const bool surface, const bool wireframe,
                const bool face_normals = false);
  void DrawMeshes();

  const std::vector<Mesh>& Meshes() const { return meshes; }
  const std::vector<Mesh>& GVDMeshes() const { return gvd_meshes; }


  // void TileSurfaceCpu();
  void TileSurfaceGpu();
  void DrawGVDSeparator();
  bool DrawEdge(const int vi, const int n_vi, const int3& p, const int3& q,
                // const oct::Direction<3>& d) const;
                const oct::Direction& d) const;
  void DrawOctree();

  void SetShowMesh(const bool show_mesh_) {
    show_mesh = show_mesh_;
  }
  void SetShowOctree(const bool show_octree_) {
    show_octree = show_octree_;
  }
  void SetShowVertices(const bool show_vertices_) {
    show_vertices = show_vertices_;
  }
  void SetShowGVD(const bool show_gvd_) {
    show_gvd = show_gvd_;
  }
  void InvertGVD() {
    for (int i = 0; i < gvd_meshes.size(); ++i) {
      gvd_meshes[i].InvertOrientation();
    }
  }
  void SetGVDColor(const bool b) {
    for (int i = 0; i < gvd_meshes.size(); ++i) {
      gvd_meshes[i].SetColor(b);
    }
  }
  void SetShowVertDistLines(const bool show_vertex_distance_lines_) {
    show_vertex_distance_lines = show_vertex_distance_lines_;
  }
  void SetShowAxis(const bool show_axis_) {
    show_axis = show_axis_;
  }
  void SetShowCellId(const bool show_cell_id_) {
    show_cell_id = show_cell_id_;
  }
  void SetMeshMode(const int mesh_mode_) {
    mesh_mode = mesh_mode_;
  }
  void SetGVDMode(const int gvd_mode_) {
    gvd_mode = gvd_mode_;
  }
  int GetGVDMode() const {
    return gvd_mode;
  }
  void SetShowStatistics(const bool b) {
    show_statistics = b;
  }

  // Returns the distance the objects traveled in octree space
  int Explode(const int max_dist);

  std::string ShotBase() const { return _shot_base; }
  void SaveTransformations();
  void ReadTransformations();
  
  virtual int ProcessArgs(int argc, char** argv);
  virtual void Init();
  virtual void Mouse(int button, int state, int x, int y);
  virtual void MouseMotion(int x, int y);
  virtual void Keyboard(unsigned char key, int x, int y);
  virtual void Special(unsigned char key, int x, int y);
  virtual void Display();

  // Returns the picked point in object coordinates
  virtual float3 Pick(int x, int y, bool& hit);
  virtual void Recenter(int x, int y);

  void ResetCenter() {
    // center = bb_full.center();
    center = bb_objects.center();
  }

  void ResizeMeshes(const size_t size);
  void ReduceGVD();

  float GetExplode();
  void SetExplode(float f);
  void IncExplode(float f = 1.1);
  void DecExplode();

  void PrintStatistics() const;

  // void SetMaxDepth(int level) { maxDepth = level; }
  void ReadMesh(const std::string& filename, bool gvd = false);
  void GenerateSurface(const oct::OctreeOptions& o =
                       oct::OctreeOptions::For3D());

  const VertexNetwork& GetVertices() const { return vertices; }

  const BoundingBox3f& bbox_objects() const { return bb_objects; }
  const BoundingBox3f& bbox_full() const { return bb_full; }

  float3& RotVec() { return rot_vec; }
  GLfloat& RotAngle() { return rot_angle; }

  const oct::OctreeOptions& Options() const { return o; }

  void WriteGvdMesh();

 private:
  void PrintHelp() const;
  void HelpString(const std::string msg, const int i) const;
  void DrawAxis();
  float3 MapMouse(GLfloat x, GLfloat y);

  void UpdateVBO();

  void SetStartSearch(int x, int y);
  void SetEndSearch(int x, int y);
  void Search(const int start, const int end);

  void SetGvdGraphDirs(
      const int label,
      // const std::map<std::pair<int,int>, list<Triangle> >& directed_tris,
      const std::map<LabelPair, list<Triangle> >& directed_tris,
      const std::vector<float3>& vertices,
      const vector<int>& mesh_dist);
  void ResetGvdMeshVertexColors();
  float3 ComputeExplodeDir(
      const int i, const set<int>& prev_ring, const vector<double3>& dirs);

  int SplitEdge(const float3& a, const float3& b,
                const int ai, const int bi,
                const double dist_a, const double dist_b,
                std::vector<int>& vert_dist,
                Mesh& mesh, std::map<Edge, int>& edge2idx);

 private:
  std::vector<Mesh> meshes;
  std::vector<Mesh> gvd_meshes;
  std::vector<Mesh> gvd_meshes_orig;
  GLuint* texture_ids;
  // std::vector<float3> centroids;
  BoundingBox3f bb_objects;
  BoundingBox3f bb_full;

  VertexNetwork vertices;
  ManagedVertexNetwork mvertices;
  // int maxDepth;
  Graph<3> _gvd_graph;
  std::vector<int> search_path;

  GLfloat mouse_x, mouse_y;
  GLfloat mouse_down_x, mouse_down_y;
  bool left_down;
  bool middle_down;
  bool right_down;
  float3 center;
  float3 rot_vec;
  GLfloat rot_angle;
  GLfloat rot_matrix[16];// = {1, 0, 0, 0,
                            // 0, 1, 0, 0,
                            // 0, 0, 1, 0,
                            // 0, 0, 0, 1};

  float2 mouse_down;
  float2 strafe;
  float2 down_strafe;

  GLfloat down_zoom;// = 1;

  bool scene_lighting;// = false;


  bool show_help;
  bool show_advanced_help;

  bool show_mesh;
  bool show_octree;
  bool show_vertices;
  bool show_gvd;
  bool show_vertex_distance_lines;
  bool show_all_vertex_distance_lines;
  bool show_vertex_id;
  bool show_axis;
  bool show_cell_id;
  bool show_statistics;
  oct::OctreeOptions o;
  oct::shared_ptr<TopoGraph> gvd_graph;

  // 0 - normal
  // 1 - show only anchor and (i<=t)-ring
  // 2 - show all exploded
  // 3 - show only anchor and (i<=t)-ring exploded
  // 4 - show only anchor and i-ring
  // 5 - show anchor in red, i-ring in blue
  // 6 - show only anchor in red, i-ring in blue
  int _explode_mode;
  // 0 - GVD bisector
  // 1 - mesh centroids
  int _explode_dir;

  int _gvd_distance_factor;
  oct::shared_array<std::vector<int> > vert_dist;
  oct::statistics<double> dist_stats;

  double _gvd_reduce_factor;

  double _path_size;

  int mesh_mode;
  int gvd_mode;

  GLuint _buffer_names[4];
  bool _buffers_valid;
  int _num_octree_edges;

  float zoom;
  float3 _picked;
  int _anchor_object;
  double _exploded_factor;
  int _max_ring;
  // Basename for screenshots
  std::string _shot_base;

  bool _random_colors;

  // Vertices to be shown as debug
  // vector<float3> _debug_verts;
};

#endif
