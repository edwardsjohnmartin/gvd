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

#include "./gvdviewer3.h"

#include <map>
#include <climits>
#include <algorithm>

#include "../octree.h"
#include "../gvd.h"
#include "./io.h"
#include "./pngio.h"
#include "../vector3.h"
#include "../tile3.h"
#include "../ambiguous.h"
#include "../statistics.h"
#include "../edge_cpp.h"

typedef oct::shared_ptr<GraphConstructor3> GraphConstructor3_h;
typedef oct::shared_ptr<TopoGraph> TopoGraph_h;

GLdouble mymodelview[16];

//------------------------------------------------------------
// Color conversion
//------------------------------------------------------------
// Conversion routines from grayscale to MATLAB's jet palette (a standard
// scivis color map)
inline double jet_interpolate( double val, double y0, double x0, double y1, double x1 ) {
    return (val-x0)*(y1-y0)/(x1-x0) + y0;
}

inline double jet_base( double val ) {
    if ( val <= -0.75 )
      return 0;
    else if ( val <= -0.25 )
      return jet_interpolate( val, 0.0, -0.75, 1.0, -0.25 );
    else if ( val <= 0.25 )
      return 1.0;
    else if ( val <= 0.75 )
      return jet_interpolate( val, 1.0, 0.25, 0.0, 0.75 );
    else return 0.0;
}

inline double jet_red( double gray ) {
    return jet_base( gray - 0.5 );
}
inline double jet_green( double gray ) {
    return jet_base( gray );
}
inline double jet_blue( double gray ) {
    return jet_base( gray + 0.5 );
}



struct CPoint {
  CPoint(float3 p_) : p(p_) {}
  CPoint(float3 p_, float3 c_) : p(p_), c(c_) {}
  float3 p;
  float3 c;
};
vector<CPoint> _debug_verts;
const float3 red = make_float3(1, 0, 0);
const float3 blue = make_float3(0, 0, 1);

namespace {
struct MergePair {
  MergePair(const float3& v_, const int i_) : v(v_), i(i_) {}
  bool operator<(const MergePair& rhs) const {
    return v < rhs.v;
  }
  bool operator==(const MergePair& rhs) const {
    return v == rhs.v;
  }
  float3 v;
  int i;
};

// Merges identical vertices and returns a map from old index to new index.
void GetMergeMap(vector<int>& v2v,
                 const vector<float3>& vertices,
                 vector<float3>& new_vertices) {
  v2v.resize(vertices.size());
  set<MergePair> added;
  int idx = 0;
  for (int i = 0; i < vertices.size(); ++i) {
    const MergePair p(vertices[i], idx);
    if (added.find(p) == added.end()) {
      added.insert(p);
      v2v[i] = idx++;
    } else {
      v2v[i] = added.find(p)->i;
    }
  }

  new_vertices.resize(added.size());
  typedef set<MergePair>::const_iterator Iter;
  for (Iter it = added.begin(); it != added.end(); ++it) {
    new_vertices.at(it->i) = it->v;
  }
}

void UpdateTrianglesWithVertexMap(const vector<int>& v2v,
                                  vector<Triangle>& triangles) {
  for (int i = 0; i < triangles.size(); ++i) {
    Triangle& t = triangles[i];
    triangles[i] = make_triangle(v2v[t.s[0]], v2v[t.s[1]], v2v[t.s[2]]);
  }
}

}

struct M3Callback {
  M3Callback(const GVDViewer3* m3_) : m3(m3_) {}
  const GVDViewer3* m3;
};

struct DrawVertexDistanceLineCallback : public M3Callback {
  DrawVertexDistanceLineCallback(const GVDViewer3* m3_) : M3Callback(m3_) {}
  bool operator()(const int vi, const int3& p) {
    return m3->DrawVertexDistanceLine(vi, p);
  }
};

struct DrawLabelCallback : public M3Callback {
  DrawLabelCallback(const GVDViewer3* m3_) : M3Callback(m3_) {}
  bool operator()(const int vi, const int3& p) {
    return m3->DrawLabel(vi, p);
  }
};

struct DrawIDCallback : public M3Callback {
  DrawIDCallback(const GVDViewer3* m3_) : M3Callback(m3_) {}
  bool operator()(const int vi, const int3& p) {
    return m3->DrawID(vi, p);
  }
};

struct DrawEdgeCallback : public M3Callback {
  DrawEdgeCallback(const GVDViewer3* m3_) : M3Callback(m3_) {}
  bool operator()(const int vi, const int n_vi,
                  const int3& p, const int3& q,
                  // const oct::Direction<3>& d) {
                  const oct::Direction& d) {
    return m3->DrawEdge(vi, n_vi, p, q, d);
  }
};

struct OctreeEdgeCallback : public M3Callback {
  OctreeEdgeCallback(const GVDViewer3* m3_)
      : M3Callback(m3_), vn(m3_->GetVertices()),
        vertices(new vector<int3>()),
        edges(new vector<Edge>()),
        labels(new vector<int>()),
        alphas(new vector<int>()),
        vi2vi(new vector<int>(m3_->GetVertices().size())),
        o(m3_->Options()) {}
  bool operator()(const int vi, const int3& p) {
    if (!o.OfInterest(vi)) return true;

    vi2vi->at(vi) = vertices->size();
    vertices->push_back(p);
    labels->push_back(vn.Label(vi));
    for (int i = 0; i < 3; ++i) {
      // const oct::Direction<3> d = oct::Direction<3>::FromAxis(i, true);
      const oct::Direction d = oct::DirectionFromAxis(i, true);
      const int n_vi = vn.Neighbor(vi, d);
      if (n_vi != -1 && o.OfInterest(n_vi)) {
        if (o.level_of_interest == -1 ||
            vn.NeighborLevel(vi, d) >= o.level_of_interest) {
          edges->push_back(make_edge(vi, n_vi));
        }
      }
    }

    alphas->push_back(vn.Alpha(vi));
    return true;
  }
  void UpdateIndices() {
    for (int i = 0; i < edges->size(); ++i) {
      Edge& e = edges->at(i);
      e = make_edge(vi2vi->at(e.s[0]), vi2vi->at(e.s[1]));
    }
  }
  const oct::VertexNetwork& vn;
  oct::shared_ptr<vector<int3> > vertices;
  oct::shared_ptr<vector<Edge> > edges;
  oct::shared_ptr<vector<int> > labels;
  oct::shared_ptr<vector<int> > alphas;
  oct::shared_ptr<vector<int> > vi2vi;
  oct::OctreeOptions o;
};

GVDViewer3::GVDViewer3(const int win_width, const int win_height)//,
                 // const float2& world_min, const float2& world_max)
    : GL3D(win_width, win_height), texture_ids(0) {//, world_min, world_max) {
  mouse_x = 0;
  mouse_y = 0;
  mouse_down_x = 0;
  mouse_down_y = 0;
  left_down = false;
  middle_down = false;
  right_down = false;
  rot_vec = make_float3(1, 0, 0);
  rot_angle = 0;
  memset(rot_matrix, 0, 16 * sizeof(GLfloat));
  for (int i = 0; i < 4; ++i) rot_matrix[4*i+i] = 1;
  // rot_matrix[16] = {1, 0, 0, 0,
  //                           0, 1, 0, 0,
  //                           0, 0, 1, 0,
  //                           0, 0, 0, 1};

  // maxDepth = 5;

  down_zoom = 1;

  // scene_lighting = false;

  show_help = false;
  show_advanced_help = false;

  show_mesh = true;
  show_octree = false;
  show_vertices = false;
  show_gvd = true;

  show_vertex_distance_lines = false;
  show_vertex_id = false;
  show_axis = false;
  show_cell_id = false;
  show_statistics = true;

  o = oct::OctreeOptions::For3D();
  o.ambiguous_max_level = -1;

  mesh_mode = 2;
  // mesh lines
  gvd_mode = 3;
  // no mesh lines
  // gvd_mode = 2;

  _buffers_valid = false;

  zoom = 1.5;
  _anchor_object = 0;
  _exploded_factor = 0.01;
  _explode_dir = 0;
  _explode_mode = 0;
  _max_ring = 1;

  _gvd_distance_factor = 16;
  _gvd_reduce_factor = 1;

  _shot_base = "shot";
  _path_size = .03;

  _random_colors = true;
}

GVDViewer3::~GVDViewer3() {
  if (texture_ids) {
    delete [] texture_ids;
  }
}

int3 GVDViewer3::Obj2Oct(const float3& v) const {
  const float3 size = bb_objects.size();
  const float max_size = max(size.s[0], max(size.s[1], size.s[2]));
  return
      convert_int3(make_float3(kWidth, kWidth, kWidth) * ((v-bb_objects.min())/max_size));
}

float3 GVDViewer3::Oct2Obj(const int3& v) const {
  const float3 vf = make_float3(v.s[0], v.s[1], v.s[2]);
  const float3 size = bb_objects.size();
  const float max_size = max(size.s[0], max(size.s[1], size.s[2]));
  const GLfloat w = kWidth;
  return (vf/w)*max_size+bb_objects.min();
}

GLfloat GVDViewer3::Oct2Obj(int dist) const {
  const float3 size = bb_objects.size();
  const float max_size = max(size.s[0], max(size.s[1], size.s[2]));
  const GLfloat ow = kWidth;
  return (dist/ow)*max_size;
}

bool GVDViewer3::DrawVertexDistanceLine(const int vi, const int3& p) const {
  if (o.cell_of_interest > -1 &&
      o.cells_of_interest.find(vi) == o.cells_of_interest.end()) {
    return true;
  }

  glLineWidth(2.0);
  glBegin(GL_LINES);
  // Closest in blue
  const int3& cp = vertices.ClosestPoint(vi);
  glColor3f(0, 0, 0.7);
  const float3 p_obj = Oct2Obj(p);
  glVertex3fv(p_obj.s);
  const float3 cp_obj = Oct2Obj(cp);
  glVertex3fv(cp_obj.s);
  glEnd();

  float3 c = make_float3(0, 0, .8);
  glColor3fv(c.s);
  Material m;
  m.set_diffuse(c);
  m.set_ambient(c);
  SetMaterial(m);
  float3 obj = Oct2Obj(cp);
  const float3 win = Obj2Win(obj);
  const float3 obj_offset = Win2Obj(make_float3(win.s[0]+5, win.s[1], win.s[2]));
  const double off = length(obj_offset-obj);
  glTranslatef(obj.s[0], obj.s[1], obj.s[2]);
  glutSolidCube(off);
  glTranslatef(-obj.s[0], -obj.s[1], -obj.s[2]);

  // stringstream ss;
  // ss << vertices.ClosestPointIndex(vi);
  // BitmapString(ss.str(), Oct2Obj(cp));

  return true;
}

void GVDViewer3::DrawVertexDistanceLines() {
  glDisable(GL_LIGHTING);
  glColor3f(0, 0.7, 0);
  glLineWidth(1.0);
  // glBegin(GL_LINES);
  oct::VisitVertices<3>(vertices, DrawVertexDistanceLineCallback(this));
  // glEnd();
  glEnable(GL_LIGHTING);
}

bool GVDViewer3::DrawLabel(const int vi, const int3& p) const {
  static const double base_size = 10;
  static const int max_dist = oct::kWidth / 2;

  if (!o.OfInterest(vi)) return true;

  if (vertices.ClosestPointIndex(vi) != -1) {
    const int dist = length(vertices.ClosestPoint(vi) - p);
    const double dist_normalized = dist/static_cast<double>(max_dist);

    const int label = vertices.Label(vi);
    const float3 red = make_float3(1, 0, 0);
    // float3 c(SetColor(label, red));
    float3 c = RandomColor(label, red);
    glColor3fv(c.s);
    Material m;
    m.set_diffuse(c);
    m.set_ambient(c);
    SetMaterial(m);

    const float3 obj = Oct2Obj(p);
    const float3 win = Obj2Win(obj);
    const float3 obj_offset =
        Win2Obj(make_float3(win.s[0]+base_size*dist_normalized+5, win.s[1], win.s[2]));
    const double off = length(obj_offset-obj);
    glTranslatef(obj.s[0], obj.s[1], obj.s[2]);
    glutSolidCube(off);
    glTranslatef(-obj.s[0], -obj.s[1], -obj.s[2]);
  }

  return true;
}

void GVDViewer3::DrawVertexLabels() {
  // glColor3f(0, 0, 0);
  // VisitVertices(vertices, DrawLabelCallback(this));

  bool use_vbo = IsExtensionSupported("GL_ARB_vertex_buffer_object");
  if (!use_vbo) {
    cerr << "VBOs not supported -- can't draw octree" << endl;
    return;
  }

  if (!_buffers_valid) {
    UpdateVBO();
  }

  glLineWidth(1.0);
  glColor3f(0.0, 0.0, 0.0);

  glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
  glDisable(GL_LIGHTING);

  glBindBuffer(GL_ARRAY_BUFFER, _buffer_names[2]);
  glVertexPointer(3, GL_FLOAT, sizeof(MeshVertexC), 0);
  glColorPointer(4, GL_UNSIGNED_BYTE, sizeof(MeshVertexC),
                 (char*)0 + 6*sizeof(GLfloat));
  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, _buffer_names[3]);

  glEnableClientState(GL_VERTEX_ARRAY);
  glEnableClientState(GL_COLOR_ARRAY);
  glDisableClientState(GL_NORMAL_ARRAY);

  int n = vertices.size();
  if (o.cell_of_interest > -1)
    n = o.cells_of_interest.size();
  // glDrawElements(GL_TRIANGLES, vertices.size()*36, GL_UNSIGNED_INT, 0);
  glDrawElements(GL_TRIANGLES, n*36, GL_UNSIGNED_INT, 0);

  glDisableClientState(GL_COLOR_ARRAY);
  glDisableClientState(GL_VERTEX_ARRAY);

  // glBegin(GL_LINES);
  // for (int i = 0; i < cb.edges->size(); ++i) {
  //   glVertex3fv(Oct2Obj((*cb.vertices)[(*cb.edges)[i][0]]));
  //   glVertex3fv(Oct2Obj((*cb.vertices)[(*cb.edges)[i][1]]));
  // }
  // glEnd();

  glEnable(GL_LIGHTING);
}

bool GVDViewer3::DrawID(const int vi, const int3& p) const {
  if (o.OfInterest(vi))
    BitmapString(vi, Oct2Obj(p), kCenterJustify, kCenterJustify);
  return true;
}

void GVDViewer3::DrawVertexIDs() {
  oct::VisitVertices<3>(vertices, DrawIDCallback(this));
}

bool GVDViewer3::DrawSeparatorVertex(const int vi, const int n_vi,
                         const int3& p, const int3& q,
                         // const oct::Direction<3>& d) const {
                         const oct::Direction& d) const {
  if (!vertices.IsBase(vi)) return true;

  // vector<pair<int3, Edge> > intersections;
  vector<pair<int3, oct::LabeledSegment<3> > > intersections;
  oct::GetIntersectionsAroundFace<3>(vi, p, vertices.CellLevel(vi), 0, 1,
                                     vertices,
                                     back_inserter(intersections), o);

  for (int i = 0; i < intersections.size(); ++i) {
    float3 obj = Oct2Obj(intersections[i].first);
    const float3 win = Obj2Win(obj);
    const float3 obj_offset =
        Win2Obj(make_float3(win.s[0]+5, win.s[1], win.s[2]));
    const double off = length(obj_offset-obj);
    glTranslatef(obj.s[0], obj.s[1], obj.s[2]);
    glutSolidCube(off);
    glTranslatef(-obj.s[0], -obj.s[1], -obj.s[2]);
  }

  return true;
}

void GVDViewer3::DrawMesh(const Mesh& mesh, const bool face_normals) {
  mesh.Display(false);
}

void GVDViewer3::DrawMesh(const Mesh& mesh,
                       const bool surface, const bool wireframe,
                       const bool face_normals) {
  if (surface)
    mesh.Display(false);
  if (wireframe)
    mesh.Display(true);
  // mesh.DisplayNormals();
}

float3 GVDViewer3::ComputeExplodeDir(
    const int i, const set<int>& prev_ring, const vector<double3>& dirs) {
  typedef set<int> Set;
  typedef set<int>::const_iterator Set_iter;
  double3 dir = make_double3(0);
  const double d = _exploded_factor;
  const int ri = gvd_graph->WhichRing(_anchor_object, i);
  if (ri == 1) {
    if (_explode_dir == 0) {
      dir = gvd_graph->Dir(_anchor_object, i) * d * ri;
    } else {
      // dir = (gvd_meshes[i].Centroid()
      //        -gvd_meshes[_anchor_object].Centroid()).unit() * d * ri;
      dir = convert_double3(normalize(gvd_meshes[i].Centroid()
                                -gvd_meshes[_anchor_object].Centroid())) * d * ri;
    }
  } else {
    const Set& one = gvd_graph->Ring(i, 1);
    vector<int> parents;
    set_intersection(one.begin(), one.end(),
                     prev_ring.begin(), prev_ring.end(),
                     back_inserter(parents));
    double3 dir_sum = make_double3(0);
    for (int j = 0; j < parents.size(); ++j) {
      const int pi = parents[j];
      const double3 pdir = dirs[pi] * gvd_graph->GetArea(pi, i);
      dir_sum += pdir;
    }
    if (_explode_dir == 0) {
      // dir = dir_sum.unit() * d * ri;
      dir = normalize(dir_sum) * d * ri;
    } else {
      // dir = (gvd_meshes[i].Centroid()
      //        -gvd_meshes[_anchor_object].Centroid()).unit() * d * ri;
      dir = convert_double3(normalize(gvd_meshes[i].Centroid()
                                -gvd_meshes[_anchor_object].Centroid())) * d * ri;
    }
  }
  return convert_float3(dir);
}

void GVDViewer3::DrawMeshes() {
  typedef set<int> Set;
  typedef set<int>::const_iterator Set_iter;

  glPolygonMode(GL_FRONT, GL_FILL);
  // glPolygonMode(GL_BACK, GL_POINT);
  glPolygonMode(GL_BACK, GL_FILL);
  // const double d = _exploded_factor;
  // for (int i = 0; i < meshes.size(); ++i) {
  int outer_ri = 0;
  Set prev_ring;
  Set ring = gvd_graph->Ring(_anchor_object, outer_ri);
  vector<double3> dirs(meshes.size());
  while (!ring.empty()) {
    for (Set_iter it = ring.begin(); it != ring.end(); ++it) {
      const int i = *it;
      bool draw = (i == _anchor_object);
      if (!draw) {
        draw = (_explode_mode == 0 || _explode_mode == 2 || _explode_mode == 5);
      }
      if (!draw) {
        // const Set& one = gvd_graph->Ring(_anchor_object, _max_ring);
        // draw = (one.find(i) != one.end());
        const int ri = gvd_graph->WhichRing(_anchor_object, i);
        if (_explode_mode == 1 || _explode_mode == 3)
          draw = (ri <= _max_ring);
        else if (_explode_mode == 4)
          draw = (ri == _max_ring);
        else if (_explode_mode == 6) {
          draw = (ri == _max_ring);
        }
      }
      if (draw) {
        glPushMatrix();
        double3 dir = make_double3(0);
        if (_explode_mode > 1 && _explode_mode < 4 && i != _anchor_object) {
          dir = convert_double3(ComputeExplodeDir(i, prev_ring, dirs));
          glTranslated(dir.s[0], dir.s[1], dir.s[2]);
        }
        dirs[i] = dir;
        const Material m = meshes[i].GetMaterial();
        if (i == _anchor_object && _explode_mode > 0) {
          meshes[i].SetMaterial(Material::FromDiffuseAmbient(make_float3(1, 0, 0)));
        } else if (_explode_mode == 5 || _explode_mode == 6) {
          const int ri = gvd_graph->WhichRing(_anchor_object, i);
          if (ri == _max_ring) {
            // const float3 color(0, 0, 1);
            const float3 color = make_float3(0.974624, 0.768307, 0.418448);
            meshes[i].SetMaterial(Material::FromDiffuseAmbient(color));
          }
        }
        DrawMesh(meshes[i], mesh_mode&2, mesh_mode&1);
        // if (i == _anchor_object) {
          meshes[i].SetMaterial(m);
        // }
      }
      glPopMatrix();
    }
    ++outer_ri;
    prev_ring = ring;
    ring = gvd_graph->Ring(_anchor_object, outer_ri);
  }
}

void GVDViewer3::SetGvdGraphDirs(
    const int label,
    const std::map<LabelPair, list<Triangle> >& directed_tris,
    const vector<float3>& mesh_vertices,
    const vector<int>& mesh_dist) {
  typedef std::set<int> Set;
  typedef Set::const_iterator Set_iter;
  typedef list<Triangle> List;
  typedef List::const_iterator List_iter;
  const Set& nbrs = gvd_graph->at(label);
  for (Set_iter it = nbrs.begin(); it != nbrs.end(); ++it) {
    const int label_adj = *it;
    LabelPair e(label, label_adj);
    const List& dir_tris = directed_tris.find(e)->second;

    // Get the distance statistics
    oct::statistics<double> dist_stats;
    for (List_iter lit = dir_tris.begin(); lit != dir_tris.end(); ++lit) {
      const Triangle& t = *lit;
      const float3& a = mesh_vertices[t.s[0]];
      const float3& b = mesh_vertices[t.s[1]];
      const float3& c = mesh_vertices[t.s[2]];
      // const double3 n = convert_double3((b-a)^(c-a));
      const double3 n = convert_double3(cross(b-a, c-a));
      if (length(n) > 0) {
        int avg_dist = 1;
        if (!mesh_dist.empty()) {
          const double ad = mesh_dist[t.s[0]];
          const double bd = mesh_dist[t.s[1]];
          const double cd = mesh_dist[t.s[2]];
          avg_dist = (int)((ad+bd+cd)/3);
        }
        dist_stats(avg_dist);
      }
    }
    const int mind = dist_stats.min();
    // const int maxd = dist_stats.max();
    const int ranged = dist_stats.range();

    // Compute directions
    double3 sum = make_double3(0);
    double area_sum = 0;
    for (List_iter lit = dir_tris.begin(); lit != dir_tris.end(); ++lit) {
      const Triangle& t = *lit;
      const float3& a = mesh_vertices[t.s[0]];
      const float3& b = mesh_vertices[t.s[1]];
      const float3& c = mesh_vertices[t.s[2]];
      // const double3 n = convert_double3((b-a)^(c-a));
      const double3 n = convert_double3(cross(b-a, c-a));
      if (length(n) > 0) {
        int avg_dist = 1;
        if (!mesh_dist.empty()) {
          const double ad = mesh_dist[t.s[0]];
          const double bd = mesh_dist[t.s[1]];
          const double cd = mesh_dist[t.s[2]];
          avg_dist = (int)((ad+bd+cd)/3);
        }
        const double dist = 1-(avg_dist-mind)/(double)ranged;
        sum += (n*pow(dist, 128));
        area_sum += length(n)/2;
      }
    }
    // cout << "sum = " << sum.unit() << endl;
    // const double3 n = sum.unit();
    const double3 n = normalize(sum);
    gvd_graph->SetDir(label, label_adj, n);
    gvd_graph->SetArea(label, label_adj, area_sum);
  }
}

void GVDViewer3::TileSurfaceGpu() {
  gvd_graph.reset(new TopoGraph());

  // Get the vertex->base mapping for 3D and the 3 2D planes
  oct::CompositeVertexVisitor<3> cv;
  oct::FindBasesVisitor<3> fbv(&vertices);
  cv.Add(&fbv);

  oct::FindBasesVisitor<3> fbv_2[] = {
    oct::FindBasesVisitor<3>(&vertices, 1, &fbv),
    oct::FindBasesVisitor<3>(&vertices, 2, &fbv),
    oct::FindBasesVisitor<3>(&vertices, 4, &fbv)
  };
  for (int i = 0; i < 3; ++i) {
    cv.Add(&fbv_2[i]);
  }

  // Every search of the octree is expensive, so compile four base finding
  // searches into one composite search.
  oct::VisitVerticesBFS<3>(vertices, cv);

  // // test
  // cout << "****** Testing find bases visitor stuff" << endl;
  // cout << "****** REMOVE!" << endl;

  // // report 3D incidence errors
  // for (int i = 0; i < vertices.size(); ++i) {
  //   for (int d = 0; d < (1<<3); ++d) {
  //     if (fbv.GetBase(i, d) != vertices.GetBase(i, 3, d)) {
  //       cout << "  i = " << i << " d = " << d << endl
  //            << " correct base = " << fbv.GetBase(i, d) << endl
  //            << " incorrect base = " << vertices.GetBase(i, 3, d) << endl;
  //     }
  //   }
  // }

  // // report 2D incidence errors
  // for (int i = 0; i < vertices.size(); ++i) {
  //   for (int normal = 0; normal < 3; ++normal) {
  //     for (int d = 0; d < (1<<2); ++d) {
  //       const int old_base = fbv_2[normal].GetBase(i, d);
  //       const int new_base = vertices.GetBase(i, 2, d, normal);
  //       // Only report false negatives -- old bases appear to not have all
  //       // incidences.
  //       if (old_base != new_base && old_base > -1) {
  //         cout << "  i = " << i << " d = " << d
  //              << " normal = " << normal << endl
  //              << " correct base = " << old_base << endl
  //              << " incorrect base = " << new_base
  //              << endl;
  //       }
  //     }
  //   }
  // }

  // // for (int i = 0; i < vertices.size(); ++i) {
  // //   for (int normal = 0; normal < 3; ++normal) {
  // //     for (int d = 0; d < (1<<2); ++d) {
  // //       if (fbv_2[normal].GetBase(i, d) > -1) {
  // //         cout << "  i = " << i
  // //              << " normal = " << normal
  // //              << " d = " << d
  // //              << " base = " << fbv_2[normal].GetBase(i, d) << endl;
  // //       }
  // //     }
  // //   }
  // // }
  // // cout << "***" << endl;
  // // for (int i = 0; i < vertices.size(); ++i) {
  // //   for (int normal = 0; normal < 3; ++normal) {
  // //     for (int d = 0; d < (1<<2); ++d) {
  // //       if (vertices.GetBase(i, 2, d, normal) > -1) {
  // //         cout << "  i = " << i
  // //              << " normal = " << normal
  // //              << " d = " << d
  // //              << " base = " << vertices.GetBase(i, 2, d, normal) << endl;
  // //       }
  // //     }
  // //   }
  // // }
  // // end test

  // Find 2D cell intersecting points
  vector<vector<LabeledIntersection> >p_cell2ints[8];
  for (int i = 1; i < 8; i=(i<<1)) {
    p_cell2ints[i].resize(vertices.size());
  }
  TileSurfaceVisitorGpu2 v2(&fbv, fbv_2,
                            p_cell2ints, &vertices, o);
  oct::VisitVerticesBFS<3>(vertices, v2);

  // Find 3D cell intersecting edges
  oct::shared_array<vector<vector<LabeledEdge> > > cell2edges(
      new vector<vector<LabeledEdge> >[meshes.size()]);
  for (int i = 0; i < meshes.size(); ++i) {
    cell2edges[i].resize(vertices.size());
  }
  oct::shared_array<vector<int3> > tri_verts(new vector<int3>[meshes.size()]);
  vert_dist.reset(new vector<int>[meshes.size()]);
  TileSurfaceVisitorGpu3 v3(
      &fbv, p_cell2ints, cell2edges.get(),
      tri_verts.get(), vert_dist.get(), gvd_graph.get(), meshes.size(), o);
  oct::VisitVerticesBFS<3>(vertices, v3);

  // Create meshes
  gvd_meshes = vector<Mesh>(meshes.size());
  for (int i = 0; i < gvd_meshes.size(); ++i) {
    if (_random_colors)
      gvd_meshes[i].AddMaterial(meshes[i].material(0));
    else {
      const float3 c = RandomColor(i, make_float3(1, 0, 0));
      const Material mat = Material::FromDiffuseAmbient(c);
      gvd_meshes[i].AddMaterial(mat);
    }
  }
  
  // Get the min/max distances etc
  // oct::statistics<double> dist_stats;
  dist_stats = oct::statistics<double>();
  for (int label = 0; label < meshes.size(); ++label) {
    for (int i = 0; i < vert_dist[label].size(); ++i) {
      // 1729
      if (vert_dist[label][i] < 0) {
        cerr << "label = " << label << " i = " << i << endl;
        throw logic_error("Invalid distance");
      }
      dist_stats(vert_dist[label][i]);
    }
  }

  map<LabelPair, list<Triangle> > directed_tris;
  const int coi = o.cell_of_interest;
  // Iterate over all vertices
  for (int vi = 0; vi < vertices.size(); ++vi) {
    if (coi > -1 && vi != coi && o.restricted_surface) {
      // Depending on options, we may only be building a portion of the
      // surface for debug purposes.
      continue;
    }
    double3 sum = make_double3(0);
    double dist_sum = 0;
    int count = 0;
    // For each label, find the edges associated with the label
    // and vertex vi and add their weights to the sum.
    for (int label = 0; label < meshes.size(); ++label) {
      // const vector<Edge>& edges = cell2edges[label][vi];
      const vector<LabeledEdge>& edges = cell2edges[label][vi];
      const int n = edges.size();
      if (n > 0) {
        // Go through all edges associated with this cell
        // TODO: find the weighted average
        for (int i = 0; i < n; ++i) {
          const Edge& le = edges[i].e;
          sum += convert_double3(tri_verts[label][le.s[0]])
              + convert_double3(tri_verts[label][le.s[1]]);
          dist_sum += vert_dist[label][le.s[0]] + vert_dist[label][le.s[1]];
          ++count;
        }
      }
    }
    if (count > 0) {
      const int3 center = convert_int3(sum / (count*2));
      const int center_dist = dist_sum / (count*2);
      for (int label = 0; label < meshes.size(); ++label) {
        // const vector<Edge>& edges = cell2edges[label][vi];
        const vector<LabeledEdge>& edges = cell2edges[label][vi];
        if (!edges.empty()) {
          const int center_idx = tri_verts[label].size();
          tri_verts[label].push_back(center);
          vert_dist[label].push_back(center_dist);
          dist_stats(center_dist);
      
          for (int i = 0; i < edges.size(); ++i) {
            // const Edge& e = edges[i];
            const Edge& e = edges[i].e;
            Triangle t = make_triangle(center_idx, e.s[0], e.s[1]);
            gvd_meshes[label].AddTriangle(t);

            const int olabel = edges[i].olabel;
            // LabelPair lpair(olabel, label);
            LabelPair lpair(label, olabel);
            // maps to second label
            directed_tris[lpair].push_back(triangle_inverted(t));
          }
        }
      }
    }
  }

  // t.restart("* 6");
  for (int label = 0; label < gvd_meshes.size(); ++label) {
    Mesh& mesh = gvd_meshes[label];
    for (int i = 0; i < tri_verts[label].size(); ++i) {
      mesh.AddVertex(Oct2Obj(tri_verts[label][i]));
    }
    const vector<float3>& mesh_vertices = mesh.vertices();
    SetGvdGraphDirs(label, directed_tris, mesh_vertices, vert_dist[label]);
  }
  for (int label = 0; label < gvd_meshes.size(); ++label) {
    gvd_meshes[label].compute_normals();
  }
  ResetGvdMeshVertexColors();

  GraphConstructor<3> gc;
  for (const Mesh& mesh : gvd_meshes) {
    const vector<float3>& vertices = mesh.vertices();
    for (const Triangle& t : mesh.triangles()) {
      gc.AddTriangle(convert_double3(vertices[t.s[0]]),
                     convert_double3(vertices[t.s[1]]),
                     convert_double3(vertices[t.s[2]]));
    }
  }
  _gvd_graph = gc.GetGraph();

  // t.restart("* 7");
  bb_full = gvd_meshes[0].bb();
  if (center == float3())
    ResetCenter();
}

void GVDViewer3::ResetGvdMeshVertexColors() {
  for (int label = 0; label < gvd_meshes.size(); ++label) {
    Mesh& mesh = gvd_meshes[label];
    for (int i = 0; i < vert_dist[label].size(); ++i) {
      const double d =
          (vert_dist[label][i]-dist_stats.min())/dist_stats.range();
      const double f = _gvd_distance_factor;
      const double gray = pow(1-d, f);
      const double r = jet_red(gray);
      const double g = jet_green(gray);
      const double b = jet_blue(gray);
      const float3 color = make_float3(r, g, b);
      mesh.SetVertexColor(i, color);
    }
  }
}

void GVDViewer3::DrawGVDSeparator() {
  glPolygonMode(GL_FRONT, GL_FILL);
  // glPolygonMode(GL_BACK, GL_POINT);
  //glPolygonMode(GL_BACK, GL_POINT);
  glEnable(GL_CULL_FACE);
  glCullFace(GL_BACK);
  glPointSize(0.0f);
  glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE);
  if (show_gvd) {
    for (int i = 0; i < gvd_meshes.size(); ++i) {
      DrawMesh(gvd_meshes[i], gvd_mode&2, gvd_mode&1, true);
    }
  }
  glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_FALSE);
  glLineWidth(1);
  glDisable(GL_CULL_FACE);

  glLineWidth(3.0);
  glColor3f(0.0, 0.0, 1.0);
  if (_path_size < 0.01) {
    glDisable(GL_LIGHTING);
    glBegin(GL_LINE_STRIP);
    for (int i = 0; i < search_path.size(); ++i) {
      glVertex3dv(_gvd_graph[search_path[i]].s);
    }
    glEnd();
  } else if (!search_path.empty()) {
    glEnable(GL_LIGHTING);
    glMaterialfv(GL_FRONT_AND_BACK, GL_EMISSION, make_float4(0, 0, 0, 1).s);
    glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, make_float4(.9, .9, .9, 1).s);
    glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, make_float4(0, 1, 0, 1).s);
    glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, make_float4(0, .4, 0, 1).s);
    glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, 40);
    float3 cur = make_float3(0);
    for (int i = 0; i < search_path.size()-1; ++i) {
      const float3& p = convert_float3(_gvd_graph[search_path[i]]);
      const float3& q = convert_float3(_gvd_graph[search_path[i+1]]);
      // const float3 v = (q-p).unit();
      const float3 v = normalize(q-p);
      const float len = length(q-p);
      float f = 0;
      if (length(p-cur) < _path_size*3) {
        f = _path_size*3 - length(p-cur);
      }
      for (; f < len; f += _path_size*3) {
        const float3 r = p + v*f;
        glTranslatef(r.s[0], r.s[1], r.s[2]);
        glutSolidSphere(_path_size, 10, 10);
        glTranslatef(-r.s[0], -r.s[1], -r.s[2]);
      }
    }
  }
}

bool GVDViewer3::DrawEdge(const int vi, const int n_vi,
                       const int3& p, const int3& q,
                       // const oct::Direction<3>& d) const {
                       const oct::Direction& d) const {
  glVertex3fv(Oct2Obj(p).s);
  glVertex3fv(Oct2Obj(q).s);
  return true;
}

struct OctVertex {
  OctVertex() {}
  OctVertex(const float3& p_)
      : p(p_) {}
  float3 p;  // point
  GLbyte padding[4];  // to make it 16 bytes
  // float3 n;  // normal
  // Vec4b c;  // color
  // GLbyte padding[4];  // to make it 32 bytes
};

void GVDViewer3::UpdateVBO() {
  OctreeEdgeCallback cb(this);
  oct::VisitVertices<3>(vertices, cb);
  cb.UpdateIndices();

  // static bool initialized = false;
  // static GLuint _buffer_names[2];
  // glDeleteBuffers(4, _buffer_names);
  glGenBuffers(4, _buffer_names);
  _buffers_valid = true;

  {
    //----------------------------------------
    // Octree lines
    //----------------------------------------
    // Send vertex buffer object
    const int n = cb.vertices->size();
    OctVertex* mverts = new OctVertex[n];
    for (int i = 0; i < n; ++i) {
      mverts[i] = OctVertex(Oct2Obj(cb.vertices->at(i)));
    }
    glBindBuffer(GL_ARRAY_BUFFER, _buffer_names[0]);
    glBufferData(
        GL_ARRAY_BUFFER, n*sizeof(OctVertex),
        mverts, GL_STREAM_DRAW);
    delete [] mverts;

    // Send index buffer object
    const int m = cb.edges->size();
    _num_octree_edges = m;
    GLuint* indices = new GLuint[m*2];
    for (int i = 0; i < m; ++i) {
      for (int j = 0; j < 2; ++j) {
        indices[i*2+j] = cb.edges->at(i).s[j];
      }
    }
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, _buffer_names[1]);
    glBufferData(
        GL_ELEMENT_ARRAY_BUFFER, m*2*sizeof(GLuint),
        indices, GL_STREAM_DRAW);
    delete [] indices;
  }

  //----------------------------------------
  // Octree vertices
  //----------------------------------------
  // Send vertex buffer object
  const int n = cb.vertices->size();
  MeshVertexC* mverts = new MeshVertexC[n*8];
  const float3 norm = make_float3(0);
  float3 red = make_float3(1, 0, 0);
  for (int i = 0; i < n; ++i) {
    const int3 v = cb.vertices->at(i);
    const int d = cb.alphas->at(i) / 8;
    const int l = cb.labels->at(i);
    const float3 c = RandomColor(l, red);
    mverts[i*8+0] = MeshVertexC(Oct2Obj(v+make_int3(-d, -d, -d)), norm, c);
    mverts[i*8+1] = MeshVertexC(Oct2Obj(v+make_int3(d, -d, -d)), norm, c);
    mverts[i*8+2] = MeshVertexC(Oct2Obj(v+make_int3(d, d, -d)), norm, c);
    mverts[i*8+3] = MeshVertexC(Oct2Obj(v+make_int3(-d, d, -d)), norm, c);
    mverts[i*8+4] = MeshVertexC(Oct2Obj(v+make_int3(-d, -d, d)), norm, c);
    mverts[i*8+5] = MeshVertexC(Oct2Obj(v+make_int3(d, -d, d)), norm, c);
    mverts[i*8+6] = MeshVertexC(Oct2Obj(v+make_int3(d, d, d)), norm, c);
    mverts[i*8+7] = MeshVertexC(Oct2Obj(v+make_int3(-d, d, d)), norm, c);
  }
  glBindBuffer(GL_ARRAY_BUFFER, _buffer_names[2]);
  glBufferData(
      GL_ARRAY_BUFFER, n*8*sizeof(MeshVertexC),
      mverts, GL_STREAM_DRAW);
  delete [] mverts;

  // Send index buffer object
  GLuint* indices = new GLuint[n*36];
  static const int inds[] = { 0, 1, 3,
                              3, 1, 2,
                              1, 5, 2,
                              2, 5, 6,
                              5, 4, 6,
                              6, 4, 7,
                              4, 0, 7,
                              7, 0, 3,
                              3, 2, 7,
                              7, 2, 6,
                              1, 0, 4,
                              1, 4, 5 };
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < 36; ++j) {
      indices[i*36+j] = inds[j] + i*8;
    }
  }
  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, _buffer_names[3]);
  glBufferData(
      GL_ELEMENT_ARRAY_BUFFER, n*36*sizeof(GLuint),
      indices, GL_STREAM_DRAW);
  delete [] indices;
}

void GVDViewer3::DrawOctree() {
  bool use_vbo = IsExtensionSupported("GL_ARB_vertex_buffer_object");
  if (!use_vbo) {
    cerr << "VBOs not supported -- can't draw octree" << endl;
    return;
  }

  if (!_buffers_valid) {
    UpdateVBO();
  }

  glLineWidth(1.0);
  glColor3f(0.0, 0.0, 0.0);

  glDisable(GL_LIGHTING);

  glBindBuffer(GL_ARRAY_BUFFER, _buffer_names[0]);
  glVertexPointer(3, GL_FLOAT, sizeof(OctVertex), 0);
  // glNormalPointer(GL_FLOAT, sizeof(OctVertex),
  //                 BUFFER_OFFSET(3*sizeof(GLfloat)));
  // glColorPointer(4, GL_UNSIGNED_BYTE, sizeof(MeshVertex),
  //                BUFFER_OFFSET(6*sizeof(GLfloat)));
  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, _buffer_names[1]);

  glEnableClientState(GL_VERTEX_ARRAY);
  glDisableClientState(GL_COLOR_ARRAY);
  glDisableClientState(GL_NORMAL_ARRAY);

  // glDrawElements(GL_LINES, cb.edges->size()*2, GL_UNSIGNED_INT, 0);
  glDrawElements(GL_LINES, _num_octree_edges*2, GL_UNSIGNED_INT, 0);

  glDisableClientState(GL_VERTEX_ARRAY);

  // glBegin(GL_LINES);
  // for (int i = 0; i < cb.edges->size(); ++i) {
  //   glVertex3fv(Oct2Obj((*cb.vertices)[(*cb.edges)[i][0]]));
  //   glVertex3fv(Oct2Obj((*cb.vertices)[(*cb.edges)[i][1]]));
  // }
  // glEnd();

  glEnable(GL_LIGHTING);

  // glDisable(GL_LIGHTING);
  // glBegin(GL_LINES);
  // VisitEdges(vertices, DrawEdgeCallback(this));
  // glEnd();
  // glEnable(GL_LIGHTING);
}

void GVDViewer3::PrintStatistics() const {
  int num_cells = 0;
  oct::level_t max_level = 0;
  for (int i = 0; i < vertices.size(); ++i) {
    if (vertices.IsBase(i)) {
      ++num_cells;
      max_level = max(max_level, vertices.CellLevel(i));
    }
  }
  {
    stringstream ss;
    ss << "Max level: " << (int)max_level;
    BitmapString(ss.str(), make_int2(2, 2));
  } {
    stringstream ss;
    ss << "Octree cells: " << num_cells;
    BitmapString(ss.str(), make_int2(2, 19));
  } {
    if (!_status.empty()) {
      BitmapString(_status, make_int2(2, 19+17));
    } else if (_picked != float3()) {
      stringstream ss;
      ss << "Picked: " << _picked << " | " << Obj2Oct(_picked);
      BitmapString(ss.str(), make_int2(2, 19+17));
    }
  }
}

void GVDViewer3::HelpString(const string msg, const int i) const {
  static const int sep = 17;
  // static const int y = 2;
  static const int y = window_height-15;
  static const int x = 2;
  // BitmapString(msg, Win2Obj(make_int2(x, y+sep*i)),
  //              kLeftJustify, kTopJustify);
  BitmapString(msg, make_int2(x, y-sep*i));
}

void GVDViewer3::PrintHelp() const {
  int i = 0;
  HelpString("h - toggle help", i++);
  HelpString("w - write GVD mesh to gvd-*.obj", i++);
  if (show_help) {
    HelpString("View", i++);
    HelpString("  m - toggle GVD", i++);
    HelpString("  g - toggle objects", i++);
    HelpString("  o - toggle octree", i++);
    HelpString("  n - invert GVD mesh orientation", i++);
    HelpString("  N - toggle GVD mesh colored by distance", i++);
    HelpString("  x/X - increment/decrement GVD color distance factor", i++);
    HelpString("  G - object mesh mode (filled, mesh/filled, mesh)", i++);
    HelpString("  a - toggle axis", i++);
    HelpString("  v - toggle vertex labels", i++);
    HelpString("  l - toggle closest point line", i++);
    HelpString("  i - toggle vertex IDs (slow)", i++);
    HelpString("  p - toggle statistics", i++);
    HelpString("  left drag - rotate", i++);
    HelpString("  Shift+left drag - strafe", i++);
    HelpString("  Ctrl+left drag - zoom", i++);
    HelpString("  Alt+left click - GVD distance from objects", i++);
    HelpString("  c/C - Adjust near clip plane", i++);
    HelpString("  u/U - Adjust far clip plane", i++);
    HelpString("  t - reset view", i++);
    HelpString("  T - save current view", i++);
    HelpString("  H - toggle advanced help", i++);
    HelpString("GVD computation", i++);
    HelpString("  f/d - increment/decrement max octree level", i++);
    HelpString("  V - toggle full subdivide", i++);
    HelpString("  B - toggle buffer", i++);
    HelpString("q - quit", i++);
  } else if (show_advanced_help) {
    HelpString("H - toggle advanced help", i++);
    HelpString("View", i++);
    // HelpString("  C/c - increment/decrement near clipping plane", i++);
    // HelpString("  U/u - increment/decrement far clipping plane", i++);
    HelpString("  P - toggle min path", i++);
    HelpString("  s - simple screenshot", i++);
    HelpString("  S - screenshot", i++);
    HelpString("  y - rotate 360 degrees with screenshots", i++);
    HelpString("  R/r - increment/decrement reduce factor", i++);
    HelpString("GVD computation", i++);
    HelpString("  b/B - increment/decrement ambiguous max level", i++);
  }
}

void GVDViewer3::Display() {
  glShadeModel(GL_SMOOTH);
  // glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
  // // glPolygonMode(GL_FRONT, GL_FILL);
  // // glPolygonMode(GL_BACK, GL_LINES);
  // glPolygonMode(GL_BACK, GL_FILL);

  glEnable(GL_LIGHT0);
  // glEnable(GL_LIGHT1);
  // const GLfloat light0_position[] = { -150, 150, 300, 1 };
  const GLfloat light0_ambient[] = { 0.1, 0.1, 0.1, 1 };
  const GLfloat light0_diffuse[] = { 0.6, 0.6, 0.6, 1 };
  const GLfloat light0_specular[] = { 0.3, 0.3, 0.3, 1 };
  glLightfv(GL_LIGHT0, GL_AMBIENT, light0_ambient);
  glLightfv(GL_LIGHT0, GL_DIFFUSE, light0_diffuse);
  glLightfv(GL_LIGHT0, GL_SPECULAR, light0_specular);
  glLightfv(GL_LIGHT1, GL_AMBIENT, light0_ambient);
  glLightfv(GL_LIGHT1, GL_DIFFUSE, light0_diffuse);
  glLightfv(GL_LIGHT1, GL_SPECULAR, light0_specular);

  glEnable(GL_NORMALIZE);

  float3 size = bb_full.size();
  const float3 eye = make_float3(
      0, 0, size.s[2]/2 + ((size.s[1]/2)/tan(20*M_PI/180.0))*1.1) * zoom;
  const float3 off = -make_float3(strafe.s[0] * size.s[0],
                           strafe.s[1] * size.s[1],
                           0);
  glMatrixMode(GL_MODELVIEW);

  glLoadIdentity();
  // gluLookAt(eye[0], eye[1], eye[2],
  //           0, 0, 0,
  //           0, 1, 0);
  gluLookAt(off.s[0]+eye.s[0], off.s[1]+eye.s[1], off.s[2]+eye.s[2],
            off.s[0], off.s[1], off.s[2],
            0, 1, 0);

  // Position the light now so it is always in the same place relative
  // to the camera.
  // const float3 light0 =
  //     float3(eye[0]-size[0]/2, eye[1]+size[1]/2, eye[2]+size[2]/2);
  const float4 light0 =
      make_float4(eye.s[0]-size.s[0]/2, eye.s[1]+size.s[1]/2, eye.s[2]+size.s[2]/2, 1);
  const float4 light1 =
      make_float4(eye.s[0]+size.s[0]/2, eye.s[1]-size.s[1]/2, eye.s[2]-size.s[2]/2, 1);

  // if (!scene_lighting) {
    glLightfv(GL_LIGHT0, GL_POSITION, light0.s);
    glLightfv(GL_LIGHT1, GL_POSITION, light1.s);
  // }

  glRotatef(rot_angle*180.0/M_PI, rot_vec.s[0], rot_vec.s[1], rot_vec.s[2]);
  glMultMatrixf(rot_matrix);

  // {
  //   // Axes
  //   glDisable(GL_LIGHTING);
  //   glColor3f(1, 0, 0);
  //   glBegin(GL_LINES);
  //   glVertex3fv(Win2Obj(float3(0, 0, 0.5)));
  //   glVertex3fv(Win2Obj(float3(100, 100, 0.5)));
  //   glEnd();
  //   // End axes
  // }

  // Position light here if we want it in the same place
  // relative to the world.
  // if (scene_lighting) {
  //   glLightfv(GL_LIGHT0, GL_POSITION, light);
  // }

  glDisable(GL_LIGHTING);
  glLineWidth(4);
  if (show_axis) {
    DrawAxis();
  }
  // glEnable(GL_LIGHTING);

  glTranslatef(-center.s[0], -center.s[1], -center.s[2]);

  if (show_octree) {
    DrawOctree();
  }
  // if (show_gvd) {
    DrawGVDSeparator();
  // }
  if (show_mesh) {
    glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE);
    DrawMeshes();
    glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_FALSE);
  }
  glGetDoublev(GL_MODELVIEW_MATRIX, mymodelview);
  // glDisable(GL_LIGHTING);
  // glColor3f(0, 0, 1);
  // glBegin(GL_LINES);
  // glVertex3fv(p_);
  // glVertex3fv(p_+v_);
  // glEnd();
  // glCube(p_, 5);
  // glColor3f(0, 1, 0);
  // glCube(q_, 15);
  // glEnable(GL_LIGHTING);

  if (show_vertex_distance_lines) {
    DrawVertexDistanceLines();
  }
  if (show_vertex_id) {
    DrawVertexIDs();
  }
  if (show_vertices) {
    DrawVertexLabels();
  }
  if (show_statistics) {
    PrintStatistics();
  }
  PrintHelp();

  // Show debug vertices
  for (int i = 0; i < _debug_verts.size(); ++i) {
    const float3 obj = Oct2Obj(convert_int3(_debug_verts[i].p));
    const float3 win = Obj2Win(obj);
    const float3 obj_offset =
        Win2Obj(make_float3(win.s[0]+5, win.s[1], win.s[2]));
    const double off = length(obj_offset-obj);
    glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE,
                 make_float4(_debug_verts[i].c, 1).s);
    glTranslatef(obj.s[0], obj.s[1], obj.s[2]);
    glutSolidCube(off);
    glTranslatef(-obj.s[0], -obj.s[1], -obj.s[2]);
  }

  // glFlush();
  // glutSwapBuffers();
}

typedef oct::index_t index_t;
void ProgressBarCallback(const int3& base_point,
                     const index_t level,
                     const index_t max_level,
                     const bool complete) {
  static const index_t w2 = oct::Level2CellWidth(1);
  if (complete) {
    cout << endl;
  } else if (level == 1) {
    if (base_point == make_int3(w2, 0)) {
      cout << "25%" << ".";
    } else if (base_point == make_int3(0, w2)) {
      cout << "50%" << ".";
    } else if (base_point == make_int3(w2, w2)) {
      cout << "75%" << ".";
    }
  } else if (level == 2) {
    cout << ".";
  }
  cout.flush();
}

void GVDViewer3::ReadMesh(const string& filename, bool gvd) {
  // Parse the obj file, compute the normals, read the textures
  Mesh mesh;
  if (!ParseObj(filename, mesh)) {
    cerr << "Error: File " << filename << " does not exist. Skipping it."
         << endl;
    return;
  }
  // mesh.print_stats();

  mesh.compute_normals();

  if (texture_ids) {
    delete [] texture_ids;
  }
  texture_ids = new GLuint[mesh.num_materials()];
  glGenTextures(mesh.num_materials(), texture_ids);

  for (int i = 0; i < mesh.num_materials(); ++i) {
    Material& material = mesh.material(i);
    material.LoadTexture(texture_ids[i]);
  }

  float3 red = make_float3(1, 0, 0);
  float3 c = make_float3(1, 1, 1);
  if (_random_colors) {
    c = RandomColor(meshes.size(), red);
  }
  Material mat = Material::FromDiffuseAmbient(c);
  if (gvd) {
    mat = Material::FromDiffuseAmbient(red);
  }
  mesh.AddMaterial(mat);

  if (gvd) {
    gvd_meshes.clear();
    gvd_meshes.push_back(mesh);
  } else {
    meshes.push_back(mesh);
  }
  bb_objects += mesh.bb();
  bb_full += mesh.bb();
  if (center == float3())
    ResetCenter();
}

void GVDViewer3::GenerateSurface(const oct::OctreeOptions& o_) {
  Timer t("** Building octree");
  Timer total("*** Total GVD computation");

  // Create the vertices and triangles
  vector<vector<float3> > all_vertices(meshes.size());
  vector<vector<Triangle> > all_triangles(meshes.size());
  for (int i = 0; i < meshes.size(); ++i) {
    const Mesh& mesh = meshes[i];
    all_vertices[i] = mesh.vertices();
    all_triangles[i].insert(all_triangles[i].end(),
                            mesh.triangles().begin(), mesh.triangles().end());
  }

  //------------------
  // Initialize OpenCL
  //------------------
#ifdef __OPEN_CL_SUPPORT__
  if (o.gpu) {
    OpenCLInit(3, o, o.opencl_log);
  }
#endif

  // Create octree with all_vertices and all_triangles
  const bool lock_subdivide = (o.ambiguous_max_level == -1);
  if (lock_subdivide) o.ambiguous_max_level = o.max_level;
  // vertices = oct::BuildOctree<3, Triangle, oct::LabeledGeometry3>(
  //     all_vertices, all_triangles, bb_objects, o);
  vertices.Clear();
  // mvertices = make_vertex_network(vertices);
  // oct::BuildOctree(
  //     all_vertices, all_triangles, bb_objects, vertices, o);
  vertices = oct::BuildOctree(
      all_vertices, all_triangles, bb_objects, o);
  // vertices = make_vertex_network(mvertices);
  if (lock_subdivide) o.ambiguous_max_level = -1;

  // Set the center
  if (o.cell_of_interest != -1) {
    vector<vector<int> > base2incident;
    base2incident = oct::ComputeBase2Incident<3>(vertices);
    const int vi = o.cell_of_interest;
    o.cells_of_interest.insert(vi);
    for (int i = 0; i < base2incident[vi].size(); ++i) {
      o.cells_of_interest.insert(base2incident[vi][i]);
    }
    for (int i = 0; i < vertices.size(); ++i) {
      for (int j = 0; j < base2incident[i].size(); ++j) {
        if (base2incident[i][j] == vi && i != vi) {
          o.cells_of_interest.insert(i);
        }
      }
    }

    oct::FindPositionsVisitor<3> fpv =
        oct::FindPositionsVisitor<3>::ForOne(vi);
    oct::VisitVertices<3>(vertices, fpv);
    center = Oct2Obj(fpv[vi]);
  } else if (o.center != -1) {
    oct::FindPositionsVisitor<3> fpv =
        oct::FindPositionsVisitor<3>::ForOne(o.center);
    oct::VisitVertices<3>(vertices, fpv);
    center = Oct2Obj(fpv[o.center]);
  } else {
    const int m = oct::kWidth/2;
    center = Oct2Obj(make_int3(m, m, m));
  }

  // if (o.gpu) {
    t.restart("** Tiling surface");
    TileSurfaceGpu();
  // } else {
  //   t.restart("** Tiling surface - cpu");
  //   TileSurfaceCpu();
  // }

  t.stop();

  if (_buffers_valid) {
    glDeleteBuffers(4, _buffer_names);
  }
  _buffers_valid = false;

  bb_full = bb_full.CenteredSquare();

  //------------------
  // Cleanup OpenCL
  //------------------
#ifdef __OPEN_CL_SUPPORT__
  if (o.gpu) {
    OpenCLCleanup();
  }
#endif
}

void GVDViewer3::SaveTransformations() {
  ofstream out("view.config");
  if (out) {
    for (int i = 0; i < 16; ++i) {
      out << rot_matrix[i] << " ";
    }
    out << endl;
    for (int i = 0; i < 3; ++i) {
      out << center.s[i] << " ";
    }
    out << endl;
    for (int i = 0; i < 2; ++i) {
      out << strafe.s[i] << " ";
    }
    out << endl;
    out << zoom;
    out << endl;
    out << bb_full.min() << endl << bb_full.max() << endl;
  }
  out.close();
}

void GVDViewer3::ReadTransformations() {
  ifstream in("view.config");
  if (in) {
    for (int i = 0; i < 16; ++i) {
      in >> rot_matrix[i];
    }
    for (int i = 0; i < 3; ++i) {
      in >> center.s[i];
    }
    for (int i = 0; i < 2; ++i) {
      in >> strafe.s[i];
    }
    in >> zoom;
    float3 minp, maxp;
    in >> minp >> maxp;
    bb_full = BoundingBox3f(minp, maxp);
    in.close();
  }
}

void PrintKeyCommands();

void PrintUsage() {
  cout << "gvd-viewer2 / gvd-viewer3" << endl;
  cout << "A 2- and 3-Dimensional Generalized Voronoi Diagram approximator\n";
  cout << "Version 1.0\n";
  cout << endl;
  cout << "Copyright 2015 by John Martin Edwards\n";
  cout << "Please send bugs and comments to edwardsjohnmartin@gmail.com\n";
  cout << endl;
  cout << "There is no warranty.\n";
  cout << endl;
  cout << "gvd-viewerx constructs a quadtree/octree that resolves between "
      "objects,\n";
  cout << "propagates a distance transform over the octree vertices, then "
      "builds a GVD\n";
  cout << "approximation using the labels and distances on the vertices. For "
      "more\n";
  cout << "information, please see:\n";
  cout << endl;
  cout << "http://sci.utah.edu/~jedwards/research/gvd/index.html\n";

  cout << endl;
  cout << "Usage: ./gvd-viewer3 [-l maxlevel] [filenames.obj]" << endl;
  cout << "  -l maxlevel        - Limit octree to maxlevel+1 levels" << endl;
  // cout << "  -x              - default to not draw GVD" << endl;
  // cout << "  -f              - get object filenames from file" << endl;
  // cout << "  --cp            - corner indices on gpu" << endl;
  // cout << "  --init-wave-cpu - disable gpu wave initialization" << endl;
  // cout << "  --test-gpu      - various tests to ensure that gpu" << endl
       // << "                    implementations are working" << endl;
  cout << endl;
  cout << endl;
  cout <<
      "If you use gvd-viewerx, I would love to hear from you. A short\n"
      "email to edwardsjohnmartin@gmail.com would be greatly appreciated.\n"
      "This is also the email to use to submit bug reports and feature\n"
      "requests." << endl;
  cout << endl;
  cout << "If you use gvd-viewerx for a publication, please cite" << endl;
  cout << endl;
  cout <<
      "   John Edwards et al, \"Approximating the Generalized Voronoi Diagram\n"
      "   of closely spaced objects,\" in Computer Graphics Forum, 2015."
       << endl;
  cout << endl;
}

int GVDViewer3::ProcessArgs(int argc, char** argv) {
  if (argc < 2) {
    PrintUsage();
    exit(0);
  }

  PrintKeyCommands();

  ReadTransformations();

  string fnfilename;

  int i = 1;
  // if (argc > 1) {
  bool stop = false;
  while (i < argc && !stop) {
    stop = true;
    if (o.ProcessArg(i, argv)) {
      stop = false;
    } else if (strcmp(argv[i], "-f") == 0) {
      ++i;
      fnfilename = argv[i];
      ++i;
      stop = false;
    } else if (strcmp(argv[i], "-s") == 0) {
      ++i;
      _shot_base = argv[i];
      ++i;
      stop = false;
    } else if (strcmp(argv[i], "-x") == 0) {
      show_gvd = false;
      ++i;
      stop = false;
    } else if (strcmp(argv[i], "--uniform-colors") == 0) {
      _random_colors = false;
      ++i;
      stop = false;
    } 
  }

  if (o.help) {
    PrintUsage();
    exit(0);
  }

  vector<string> filenames;
  for (; i < argc; ++i) {
    filenames.push_back(argv[i]);
  }

  if (!fnfilename.empty()) {
    ifstream in(fnfilename.c_str());
    if (!in) {
      cerr << "Failed to open " << fnfilename << endl;
      exit(1);
    }
    while (!in.eof()) {
      string fn;
      in >> fn;
      if (!fn.empty() && fn[0] != '#') {
        filenames.push_back(fn);
      }
    }
  }

  int num_tris = 0;
  bb_objects = BoundingBox3f();
  bb_full = BoundingBox3f();
  // for (; i < argc; ++i) {
    // string filename(argv[i]);
  for (int j = 0; j < filenames.size(); ++j) {
    const string filename = filenames[j];
    cout << "Reading " << filename << endl;
    ReadMesh(filename);
    num_tris += meshes[meshes.size()-1].num_triangles();
    // cout << "  bb = " << meshes[meshes.size()-1].bb() << endl;
  }
  bb_objects = bb_objects.CenteredSquare().ScaleCentered(o.bb_scale);

  cout << "Number of objects: " << meshes.size() << endl;
  cout << "Number of triangles: " << num_tris << endl;

  // o.simple_q = true;
  o.simple_q = false;
  GenerateSurface(o);

  return 0;
}

void GVDViewer3::Init() {
  glEnable(GL_DEPTH_TEST);
  glDepthMask(GL_TRUE);
  glDepthFunc(GL_LEQUAL);
  glClear(GL_DEPTH_BUFFER_BIT | GL_COLOR_BUFFER_BIT);

  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
  // glClearColor(1.0f, 1.0f, 1.0f, 1.0f);

  // resize the window
  // window_aspect = window_width/static_cast<float>(window_height);

  // glMatrixMode(GL_PROJECTION);
  // glLoadIdentity();
  // gluPerspective(40.0, window_aspect, 1, 1500);
}

void GVDViewer3::DrawAxis() {
  const float3 c = make_float3(0, 0, 0);
  // const float L = 20;
  const float L = 0.5;
  const float3 X = make_float3(L, 0, 0);
  const float3 Y = make_float3(0, L, 0);
  const float3 Z = make_float3(0, 0, L);

  glBegin(GL_LINES);
  glColor3f(1, 0, 0);
  glVertex3fv(c.s);
  glVertex3fv((c+X).s);
  glColor3f(0, 1, 0);
  glVertex3fv(c.s);
  glVertex3fv((c+Y).s);
  glColor3f(0, 0, 1);
  glVertex3fv(c.s);
  glVertex3fv((c+Z).s);
  glEnd();

  BitmapString("x", c+X, 5, 0, GLUT_BITMAP_8_BY_13);
  BitmapString("y", c+Y, 5, 0, GLUT_BITMAP_8_BY_13);
  BitmapString("z", c+Z, 5, 0, GLUT_BITMAP_8_BY_13);
}

float3 GVDViewer3::MapMouse(GLfloat x, GLfloat y) {
  if (x*x + y*y > 1) {
    const GLfloat len = sqrt(x*x + y*y);
    x = x/len;
    y = y/len;
  }
  const GLfloat z = sqrt(max(0.0f, 1 - x*x - y*y));
  return make_float3(x, y, z);
}

// static const int ROTATE_BUTTON = GLUT_LEFT_BUTTON;
// static const int STRAFE_BUTTON = GLUT_RIGHT_BUTTON;
// static const int ZOOM_BUTTON = GLUT_MIDDLE_BUTTON;
enum Mode { MODE_INVALID, MODE_ROTATE, MODE_STRAFE, MODE_ZOOM,
            MODE_ZOOM_IN, MODE_ZOOM_OUT, // wheel mouse
            MODE_START_PATH, MODE_END_PATH, MODE_PICK,
            MODE_RECENTER };

Mode GetMode(int button, int state, int x, int y) {
  Mode mode = MODE_INVALID;
  int mod = glutGetModifiers();

  if (button == GLUT_LEFT_BUTTON) {
    if (mod == GLUT_ACTIVE_SHIFT) {
      mode = MODE_STRAFE;
    } else if (mod == GLUT_ACTIVE_CTRL) {
      mode = MODE_ZOOM;
    } else if (mod == GLUT_ACTIVE_ALT) {
      mode = MODE_PICK;
    } else if (mod == (GLUT_ACTIVE_ALT | GLUT_ACTIVE_SHIFT)) {
      mode = MODE_START_PATH;
    } else if (mod == (GLUT_ACTIVE_ALT | GLUT_ACTIVE_CTRL)) {
      mode = MODE_END_PATH;
    } else {
      mode = MODE_ROTATE;
    }
  } else if (button == GLUT_MIDDLE_BUTTON) {
    mode = MODE_ZOOM;
  } else if (button == GLUT_RIGHT_BUTTON) {
    mode = MODE_STRAFE;
  } else if (button == 3) { // && state != GLUT_UP) {
    mode = MODE_ZOOM_IN;
  } else if (button == 4) { // && state != GLUT_UP) {
    mode = MODE_ZOOM_OUT;
  }

  if (mode == MODE_INVALID)
    throw logic_error("Invalid mode");

  return mode;
}

string MouseCommands() {
  stringstream ss;
  ss << "Mouse commands:" << endl;
  ss << "  rotate - left" << endl;
  ss << "  strafe - right / shift+left" << endl;
  ss << "  zoom - middle / ctrl+left" << endl;
  ss << "  begin path - alt+shift+left" << endl;
  ss << "  end path - alt+ctrl+left" << endl;
  // ss << "  ctrl+left - pick location" << endl;
  // ss << "  ctrl+right - set pick as center of rotation" << endl;
  ss << endl;
  return ss.str();
}

// string MouseCommands() {
//   stringstream ss;
//   ss << "Mouse commands:" << endl;
//   ss << "  left - rotate" << endl;
//   ss << "  right - strafe" << endl;
//   ss << "  middle - zoom" << endl;
//   ss << "  shift+left - begin path" << endl;
//   ss << "  shift+right - end path" << endl;
//   ss << "  ctrl+left - pick location" << endl;
//   ss << "  ctrl+right - set pick as center of rotation" << endl;
//   ss << endl;
//   return ss.str();
// }

static const GLfloat zfactor = 2;

float3 GVDViewer3::Pick(int x, int y, bool& hit) {
  float2 m = make_float2(x, window_height-y);
  // In window coordinates, z==1 --> far
  //                        z==0 --> near

  glMatrixMode(GL_MODELVIEW_MATRIX);
  glLoadMatrixd(mymodelview);

  float3 p = Win2Obj(make_float3(m, 0));
  float3 v = Win2Obj(make_float3(m, 1)) - p;
  float mint = 99999;
  // float3 picked = make_float3(0);
  if (show_mesh) {
    for (int i = 0; i < meshes.size(); ++i) {
      const float t = meshes[i].pick(p, v);
      mint = std::min(mint, t);
    }
  }
  if (show_gvd) {
    for (int i = 0; i < gvd_meshes.size(); ++i) {
      const float t = gvd_meshes[i].pick(p, v);
      mint = std::min(mint, t);
    }
  }
  float3 q = make_float3(0);
  hit = (mint < 99999);
  if (hit) {
    // q = p+mint*v;
    q = p+v*mint;
  }
  return q;
}

int GVDViewer3::PickGvdVertex(int x, int y, int* label) {
  bool hit;
  const float3 p = Pick(x, y, hit);
  if (!hit) return -1;

  double min_dist = numeric_limits<double>::max();
  int vi = -1;
  // for (const Mesh& mesh : gvd_meshes) {
  for (int l = 0; l < gvd_meshes.size(); ++l) {
    const Mesh& mesh = gvd_meshes[l];
    const vector<float3>& vertices = mesh.vertices();
    for (int i = 0; i < vertices.size(); ++i) {
      const double d = length2(p-vertices[i]);
      if (d < min_dist) {
        min_dist = d;
        vi = i;
        *label = l;
      }
    }
  }
  return vi;
}

void GVDViewer3::Recenter(int x, int y) {
  bool hit;
  const float3 p = Pick(x, y, hit);
  if (hit) {
    center = p;
    cout << "Center updated to " << center << endl;
  }
}

void GVDViewer3::Search(const int start, const int end) {
  search_path.clear();
  _gvd_graph.Dijkstra(start, end, back_inserter(search_path));
  cout << "Searched: " << search_path.size() << endl;
  reverse(search_path.begin(), search_path.end());
  glutPostRedisplay();
}

void GVDViewer3::SetStartSearch(int x, int y) {
  cout << "Starting search" << endl;
  bool hit;
  const double3 p = convert_double3(Pick(x, y, hit));
  if (!hit) return;

  double min_dist = numeric_limits<double>::max();
  int start = -1;
  const std::vector<double3>& vertices = _gvd_graph.GetVertices();
  for (int i = 0; i < vertices.size(); ++i) {
    // const double d = (p-vertices[i]).norm2();
    const double d = length2(p-vertices[i]);
    if (d < min_dist) {
      min_dist = d;
      start = i;
    }
  }
  if (start > -1) {
    if (search_path.empty()) {
      search_path.push_back(start);
    } else {
      Search(start, search_path.back());
    }
  }
  // cout << "started search with " << start << endl;
}

void GVDViewer3::SetEndSearch(int x, int y) {
  cout << "Ending search" << endl;
  bool hit;
  const double3 p = convert_double3(Pick(x, y, hit));
  if (!hit) return;

  double min_dist = numeric_limits<double>::max();
  int end = -1;
  const std::vector<double3>& vertices = _gvd_graph.GetVertices();
  for (int i = 0; i < vertices.size(); ++i) {
    // const double d = (p-vertices[i]).norm2();
    const double d = length2(p-vertices[i]);
    if (d < min_dist) {
      min_dist = d;
      end = i;
    }
  }
  if (end > -1) {
    if (search_path.empty()) {
      search_path.push_back(end);
    } else {
      Search(search_path.front(), end);
    }
  }
  // cout << "ended search with " << end << endl;
}

// void GVDViewer3::Mouse(int button, int state, int x, int y) {
//   mouse_x = 2 * x / (GLfloat)window_width - 1;
//   mouse_y = 2 * (window_height-y) / (GLfloat)window_height - 1;
//   const Mode mode = GetMode(button, state, x, y);
//   if (button == ROTATE_BUTTON) {
//     if (state == GLUT_DOWN) {
//       if (!glutGetModifiers()) {
//         left_down = true;
//         mouse_down_x = mouse_x;
//         mouse_down_y = mouse_y;
//         mouse_down = float2(mouse_x, mouse_y);
//       // } else if (glutGetModifiers() == GLUT_ACTIVE_ALT) {
//       // doesn't seem to work
//       //   Recenter(x, y);
//       } else if (glutGetModifiers() == GLUT_ACTIVE_SHIFT) {
//         SetStartSearch(x, y);
//       } else if (glutGetModifiers() == GLUT_ACTIVE_CTRL) {
//         bool hit;
//         const float3 p = Pick(x, y, hit);
//         if (hit)
//           _picked = p;
//         else
//           _picked = float3();
//       }
//     } else {
//       left_down = false;
//       glMatrixMode(GL_MODELVIEW);
//       glLoadIdentity();
//       glRotatef(rot_angle*180.0/M_PI, rot_vec[0], rot_vec[1], rot_vec[2]);
//       glMultMatrixf(rot_matrix);
//       glGetFloatv(GL_MODELVIEW_MATRIX, rot_matrix);
//       rot_angle = 0;

//     }
//   } else if (button == STRAFE_BUTTON) {
//     if (state == GLUT_DOWN) {
//       if (glutGetModifiers() == GLUT_ACTIVE_SHIFT) {
//         SetEndSearch(x, y);
//       } else if (glutGetModifiers() == GLUT_ACTIVE_CTRL) {
//         Recenter(x, y);
//       } else {
//         middle_down = true;
//         mouse_down_x = mouse_x;
//         mouse_down_y = mouse_y;
//         mouse_down = float2(mouse_x, mouse_y);
//         down_strafe = strafe;
//       }
//     } else {
//       middle_down = false;
//     }
//   } else if (button == ZOOM_BUTTON) {
//     if (state == GLUT_DOWN) {
//       right_down = true;
//       mouse_down_x = mouse_x;
//       mouse_down_y = mouse_y;
//       mouse_down = float2(mouse_x, mouse_y);
//       down_zoom = zoom;
//     } else {
//       right_down = false;
//     }
//   } else if ((button == 3) || (button == 4)) {
//     // It's a wheel event
//     // Each wheel event reports like a button click, GLUT_DOWN then GLUT_UP
//     if (state == GLUT_UP) {
//       // Disregard redundant GLUT_UP events
//     } else {
//       // 3 is scroll up, 4 is scroll down
//       if (button == 3)
//         zoom = zoom * 1.1;
//       else
//         zoom = zoom * 0.9;
//     }
//   }

//   glutPostRedisplay();
// }

void GVDViewer3::Mouse(int button, int state, int x, int y) {
  mouse_x = 2 * x / (GLfloat)window_width - 1;
  mouse_y = 2 * (window_height-y) / (GLfloat)window_height - 1;
  const Mode mode = GetMode(button, state, x, y);
  if (mode == MODE_ROTATE) {
    if (state == GLUT_DOWN) {
      left_down = true;
      mouse_down_x = mouse_x;
      mouse_down_y = mouse_y;
      mouse_down = make_float2(mouse_x, mouse_y);
    } else {
      left_down = false;
      glMatrixMode(GL_MODELVIEW);
      glLoadIdentity();
      glRotatef(rot_angle*180.0/M_PI, rot_vec.s[0], rot_vec.s[1], rot_vec.s[2]);
      glMultMatrixf(rot_matrix);
      glGetFloatv(GL_MODELVIEW_MATRIX, rot_matrix);
      rot_angle = 0;
    }
  } else if (mode == MODE_STRAFE) {
    if (state == GLUT_DOWN) {
      middle_down = true;
      mouse_down_x = mouse_x;
      mouse_down_y = mouse_y;
      mouse_down = make_float2(mouse_x, mouse_y);
      down_strafe = strafe;
    } else {
      middle_down = false;
    }
  } else if (mode == MODE_ZOOM) {
    if (state == GLUT_DOWN) {
      right_down = true;
      mouse_down_x = mouse_x;
      mouse_down_y = mouse_y;
      mouse_down = make_float2(mouse_x, mouse_y);
      down_zoom = zoom;
    } else {
      right_down = false;
    }
  } else if (mode == MODE_ZOOM_IN) {
    zoom = zoom * 1.1;
  } else if (mode == MODE_ZOOM_OUT) {
    zoom = zoom * 0.9;
  } else if (mode == MODE_START_PATH) {
    SetStartSearch(x, y);
  } else if (mode == MODE_END_PATH) {
    SetEndSearch(x, y);
  } else if (mode == MODE_PICK) {
    bool hit;
    const float3 p = Pick(x, y, hit);
    if (hit) {
      _picked = p;
      int label;
      const int vi = PickGvdVertex(x, y, &label);
      stringstream ss;
      ss << "Distance = " << Oct2Obj(vert_dist[label][vi]);
      _status = ss.str();
    }
    else {
      _picked = float3();
      _status = "";
    }
  } else if (mode == MODE_RECENTER) {
    Recenter(x, y);
  }

  glutPostRedisplay();
}

void GVDViewer3::MouseMotion(int x, int y) {
  mouse_x = 2 * x / (GLfloat)window_width - 1;
  mouse_y = 2 * (window_height-y) / (GLfloat)window_height - 1;
  float2 mouse = make_float2(mouse_x, mouse_y);
  if (left_down) {
    const float3 down_v = MapMouse(mouse_down_x, mouse_down_y);
    const float3 v = MapMouse(mouse_x, mouse_y);
    // rot_vec = (down_v ^ v).unit();
    // rot_vec = normalize(down_v ^ v);
    rot_vec = normalize(cross(down_v, v));
    // rot_angle = acos((down_v * v) / length(v));
    rot_angle = acos(dot(down_v, v) / length(v));
  } else if (middle_down) {
    strafe = down_strafe + (mouse - mouse_down);
  } else if (right_down) {
    zoom = down_zoom * pow(zfactor, mouse_y - mouse_down_y);
  }

  glutPostRedisplay();
}

void PrintKeyCommands() {
}

void GVDViewer3::WriteGvdMesh() {
  cout << "Writing gvd mesh to file gvd.obj" << endl;
  for (int i = 0; i < gvd_meshes.size(); ++i) {
    {
      stringstream ss;
      ss << "gvd-" << (i+1) << ".obj";
      ofstream out(ss.str().c_str());
      WriteObj(out, gvd_meshes[i]); 
    } {
      stringstream ss;
      ss << "gvd-" << (i+1) << ".ply";
      ofstream out(ss.str().c_str());
      WritePly(out, gvd_meshes[i]); 
    }
  }
  cout << "Gvd mesh written" << endl;
}

void GVDViewer3::Keyboard(unsigned char key, int x, int y) {
  const float zoom_factor = 1.2;
  bool redisplay = true;
  switch (key) {
    case 'h':
      show_help = !show_help;
      show_advanced_help = false;
      break;
    case 'H':
      show_advanced_help = !show_advanced_help;
      show_help = false;
      break;

    // the next few are shift+[1,2,3...]
    case '!':
      ResizeMeshes(1);
      break;
    case '@': 
      ResizeMeshes(2);
      break;
    case '#':
      ResizeMeshes(3);
      break;
    case '$':
      ResizeMeshes(4);
      break;
    case '%':
      ResizeMeshes(5);
      break;
    case '^':
      ResizeMeshes(6);
      break;
    case '&':
      ResizeMeshes(7);
      break;
    case '*':
      ResizeMeshes(8);
      break;
    case '(':
      ResizeMeshes(9);
      break;
    case '0':
      _max_ring = 0;
      break;
    case '1':
      _max_ring = 1;
      break;
    case '2':
      _max_ring = 2;
      break;
    case '3':
      _max_ring = 3;
      break;
    case '4':
      _max_ring = 4;
      break;
    case '5':
      _max_ring = 5;
      break;
    case '6':
      _max_ring = 6;
      break;
    case 'o':
      show_octree = !show_octree;
      break;
    case 'g':
      show_mesh = !show_mesh;
      break;
    case 'm':
      show_gvd = !show_gvd;
      break;
    case 'G':
      mesh_mode = (mesh_mode%3) + 1;
      break;
    case 'M':
      gvd_mode = (gvd_mode%3) + 1;
      break;
    case 'z':
      zoom *= 1/zoom_factor;
      break;
    case 'Z':
      zoom *= zoom_factor;
      break;
    case 'i':
      show_vertex_id = !show_vertex_id;
      break;
    case 'a':
      show_axis = !show_axis;
      break;
    case 'n': {
      for (int i = 0; i < gvd_meshes.size(); ++i) {
        gvd_meshes[i].InvertOrientation();
        if (!gvd_meshes_orig.empty()) {
          gvd_meshes_orig[i].InvertOrientation();
        }
      }
      break;
    }
    case 'N': {
      for (int i = 0; i < gvd_meshes.size(); ++i) {
        gvd_meshes[i].SetColor(!gvd_meshes[i].IsColor());
      }
      break;
    }
    case 'x':
      _gvd_distance_factor *= 2;
      cout << "_gvd_distance_factor = " << _gvd_distance_factor << endl;
      ResetGvdMeshVertexColors();
      break;
    case 'X':
      _gvd_distance_factor /= 2;
      cout << "_gvd_distance_factor = " << _gvd_distance_factor << endl;
      ResetGvdMeshVertexColors();
      break;
    case 'l':
      show_vertex_distance_lines = !show_vertex_distance_lines;
      break;
    case 'v':
      show_vertices = !show_vertices;
      break;
    case 'Q':
      o.simple_q = !o.simple_q;
      GenerateSurface(o);
      break;
    case 'f':
      if (o.max_level < oct::kMaxLevel) {
        ++o.max_level;
      }
      cout << "Max octree level = " << static_cast<int>(o.max_level) << endl;
      GenerateSurface(o);
      break;
    case 'd':
      o.max_level = max(o.max_level-1, 0);
      cout << "Max octree level = " << static_cast<int>(o.max_level) << endl;
      GenerateSurface(o);
      break;
    case 'w':
      WriteGvdMesh();
      break;
    case 'b':
      o.ambiguous_max_level++;
      cout << "Max ambiguous subdivision level = "
           << o.ambiguous_max_level << endl;
      GenerateSurface(o);
      break;
    case 'B':
      o.ambiguous_max_level--;
      cout << "Max ambiguous subdivision level = "
           << o.ambiguous_max_level << endl;
      GenerateSurface(o);
      break;
    case 'V':
      o.full_subdivide = !o.full_subdivide;
      GenerateSurface(o);
      break;
    // case 'B':
    //   o.make_buffer = !o.make_buffer;
    //   GenerateSurface(o);
    //   break;
    case 'e':
      _explode_mode = (_explode_mode+1)%7;
      cout << "Explode mode = " << _explode_mode << endl;
      break;
    case 'E':
      _explode_dir = (_explode_dir+1)%2;
      break;
    case 'J':
      _path_size *= 1.1;
      cout << "Path size = " << _path_size << endl;
      break;
    case 'K':
      _path_size *= 0.9;
      cout << "Path size = " << _path_size << endl;
      break;
    case 'p':
      show_statistics = !show_statistics;
      break;
    case 'R':
      // ResetCenter();
      // cout << "Center reset" << endl;
      _gvd_reduce_factor *= 2;
      ReduceGVD();
      break;
    case 'r':
      _gvd_reduce_factor *= 0.5;
      ReduceGVD();
      break;
    case 't':
      ReadTransformations();
      cout << "View reset" << endl;
      break;
    case 'T':
      SaveTransformations();
      cout << "Saved current view" << endl;
      break;
    // case 'P':
    //    writePngImage("sceenshot.png", window_width, window_height);
    //    break;
    default:
      redisplay = false;
      break;

  }

  if (redisplay) {
    glutPostRedisplay();
  }
}

// v0 is good, v1 is bad
pair<float3,double> Interpolate(
    const float3& v0, const float3& v1,
    const double d0, const double d1,
    const double f, const oct::statistics<double>& dist_stats) {
  const float3 p = v0*(d1-f)/(d1-d0) + v1*(f-d0)/(d1-d0);
  double d = d0*(d1-f)/(d1-d0) + d1*(f-d0)/(d1-d0);
  d = d * dist_stats.range() + dist_stats.min();
  return make_pair(p, d);
}

// Splits a and b.
//
//      a-----------------b   -->   a--------c--------b
//
// edge2idx maps an edge to a split point.  If there is an entry, then that
// point is returned.
int GVDViewer3::SplitEdge(const float3& a, const float3& b,
              const int ai, const int bi,
              const double dist_a, const double dist_b,
              std::vector<int>& vert_dist,
              Mesh& mesh, map<Edge, int>& edge2idx) {
  const double f = _gvd_reduce_factor;
  const Edge ab = make_edge(ai, bi);
  int ci;
  if (edge2idx.find(ab) == edge2idx.end()) {
    float3 c = make_float3(0);
    double dist_c;
    oct::tie(c, dist_c) =
        Interpolate(a, b, dist_a, dist_b, f, dist_stats);
    ci = mesh.vertices().size();
    mesh.AddVertex(c);
    vert_dist.push_back(dist_c);
    edge2idx[ab] = ci;
  } else {
    ci = edge2idx[ab];
  }
  return ci;
}

// Cuts out all triangles in the mesh that are above the average distance
// from an object
void GVDViewer3::ReduceGVD() {
  // hack
  if (gvd_meshes_orig.empty()) {
    gvd_meshes_orig = gvd_meshes;
  }

  const double f = _gvd_reduce_factor;
  for (int label = 0; label < gvd_meshes_orig.size(); ++label) {
    const Mesh& orig_mesh = gvd_meshes_orig[label];
    const vector<Triangle>& orig_triangles = orig_mesh.triangles();
    const vector<float3>& orig_vertices = orig_mesh.vertices();
    vert_dist[label].resize(orig_vertices.size());
    Mesh mesh;
    mesh.AddVertices(orig_vertices.begin(), orig_vertices.end());
    map<Edge, int> edge2idx;
    for (int i = 0; i < orig_triangles.size(); ++i) {
      const Triangle& t = orig_triangles[i];
      bool good[3];
      double d[3];
      for (int j = 0; good && j < 3; ++j) {
        const int k = t.s[j];
        d[j] = (vert_dist[label][k]-dist_stats.min())/dist_stats.range();
        good[j] = (d[j] < f);
      }
      const int num_good = (good[0]?1:0) + (good[1]?1:0) + (good[2]?1:0);
      if (num_good == 3) {
        mesh.AddTriangle(t);
      } else if (num_good == 2) {
        // Cut into two triangles
        // v0 and v1 are the good ones.  v2 is bad
        int idx2;
        for (int j = 0; j < 3; ++j) {
          if (!good[j]) idx2 = j;
        }
        const int idx0 = (idx2+1)%3;
        const int idx1 = (idx2+2)%3;
        const int vi0 = t.s[idx0];
        const int vi1 = t.s[idx1];
        const double d0 = d[idx0];
        const double d1 = d[idx1];
        const double d2_old = d[idx2];
        const float3 v0 = mesh.vertices()[t.s[idx0]];
        const float3 v1 = mesh.vertices()[t.s[idx1]];

        const int vi2 = SplitEdge(
            v1, mesh.vertices()[t.s[idx2]], vi1, t.s[idx2], d1, d2_old,
            vert_dist[label], mesh, edge2idx);

        const int vi3 = SplitEdge(
            v0, mesh.vertices()[t.s[idx2]], vi0, t.s[idx2], d0, d2_old,
            vert_dist[label], mesh, edge2idx);

        mesh.AddTriangle(make_triangle(vi0, vi1, vi2));
        mesh.AddTriangle(make_triangle(vi0, vi2, vi3));
      } else if (num_good == 1) {
        // Chop triangle
        int idx0;
        for (int j = 0; j < 3; ++j) {
          if (good[j]) idx0 = j;
        }
        const int idx1 = (idx0+1)%3;
        const int idx2 = (idx0+2)%3;
        const int vi0 = t.s[idx0];
        const double d0 = d[idx0];
        const double d1_old = d[idx1];
        const double d2_old = d[idx2];
        const float3 v0 = mesh.vertices()[vi0];

        const int vi1 = SplitEdge(
            v0, mesh.vertices()[t.s[idx1]], vi0, t.s[idx1], d0, d1_old,
            vert_dist[label], mesh, edge2idx);

        const int vi2 = SplitEdge(
            v0, mesh.vertices()[t.s[idx2]], vi0, t.s[idx2], d0, d2_old,
            vert_dist[label], mesh, edge2idx);

        mesh.AddTriangle(make_triangle(vi0, vi1, vi2));
      }
    }
    mesh.compute_normals();
    mesh.AddMaterial(gvd_meshes[label].material(0));
    mesh.SetColor(gvd_meshes[label].IsColor());
    gvd_meshes[label] = mesh;
  }
  ResetGvdMeshVertexColors();
  cout << "Mesh reduced " << _gvd_reduce_factor << endl;
}

void GVDViewer3::ResizeMeshes(const size_t size) {
  meshes.resize(min(meshes.size(), size));

  const BoundingBox3f bb_full_bak = bb_full;
  const BoundingBox3f bb_objects_bak = bb_objects;
  const float3 old_center = center;

  // ReadTransformations();
  // SaveTransformations();

  GenerateSurface(o);
  center = old_center;
  bb_full = bb_full_bak;
  bb_objects = bb_objects_bak;
}

float GVDViewer3::GetExplode() {
  return _exploded_factor;
}
void GVDViewer3::SetExplode(float f) {
  _exploded_factor = f;
}

void GVDViewer3::IncExplode(float f) {
  _exploded_factor *= f;
}

void GVDViewer3::DecExplode() {
  _exploded_factor *= 0.9;
}

int GVDViewer3::Explode(const int max_dist) {
  _anchor_object = 0;
  const float3 old_center = center;
  vector<float3> dirs(meshes.size());
  int min_dist = 99999;
  for (int target = 1; target < meshes.size(); ++target) {
    set<int> prev_ring;
    vector<double3> temp_dirs;
    float3 dir = make_float3(0, 1, 0);
    if (vert_dist) {
      dir = ComputeExplodeDir(target, prev_ring, temp_dirs);
      // dir = dir.unit();
      dir = normalize(dir);
    }
    dirs[target] = dir;
    // Find the min dist
    if (vert_dist) {
      for (int i = 0; i < vert_dist[target].size(); ++i) {
        min_dist = min(min_dist, vert_dist[target][i]);
      }
    }
  }
  const int dist = std::min(max_dist, (int)(min_dist*0.7));
  for (int target = 1; target < meshes.size(); ++target) {
    float3 dir = dirs[target];
    dir *= Oct2Obj(dist);
    meshes[target].Translate(dir);
  }

  // bb_objects = BoundingBox3f();
  bb_full = BoundingBox3f();
  for (int i = 0; i < meshes.size(); ++i) {
    // bb_objects += meshes[i].bb();
    bb_full += meshes[i].bb();
  }
  // bb_objects = bb_objects.CenteredSquare();

  GenerateSurface(o);
  center = old_center;

  return dist;
}

void GVDViewer3::Special(unsigned char key, int x, int y) {
  if (key == GLUT_KEY_LEFT) {
    _anchor_object = (_anchor_object+meshes.size()-1) % meshes.size();
    cout << "target object = " << _anchor_object << endl;
    for (int i = 0; i < 4; ++i) {
      cout << "  " << i << "-ring: "
         << gvd_graph->Ring(_anchor_object, i).size() << endl;
    }
  } else if (key == GLUT_KEY_RIGHT) {
    _anchor_object = (_anchor_object+1) % meshes.size();
    for (int i = 0; i < 4; ++i) {
      cout << "  " << i << "-ring: "
         << gvd_graph->Ring(_anchor_object, i).size() << endl;
    }
  } else if (key == GLUT_KEY_DOWN) {
    DecExplode();
  } else if (key == GLUT_KEY_UP) {
    IncExplode();
  }
  glutPostRedisplay();
}
