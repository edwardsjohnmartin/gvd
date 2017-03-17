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

#include "./octviewer2.h"

#include <iostream>
#include <fstream>
#include <utility>
#include <vector>
#include <map>
#include <limits>
#include <algorithm>

#include "../opencl/edge.h"
#include "../opencl/geom.h"

#include "../octree.h"
#include "../search.h"
#include "../vector2.h"
#include "../gvd.h"
#include "../graph.h"

#include "../karras.h"
#include "../OctreeUtils.h"

#ifdef SINGLE
#define REAL float
#else /* not SINGLE */
#define REAL double
#endif /* not SINGLE */
#define VOID void
extern "C" {
#include "../jrs-triangle.h"
}

Karras::OctCell fnode;

typedef oct::LabeledSegment<2> LabeledSegment;

using namespace std;

template <typename T>
T MeanValue(const float2& p, const vector<float2>& polygon,
                         const vector<T>& color) {
  const int n = polygon.size();

  for (int i = 0; i < n; ++i) {
    const int j = (i+1)%n;
    if (p == polygon[i]) return color[i];
    // if (p == polygon[j]) return color[j];
    for (int k = 0; k < 2; ++k) {
      if (p.s[k] == polygon[i].s[k] && p.s[k] == polygon[j].s[k]) {
        if (p.s[1-k] > polygon[i].s[1-k] && p.s[1-k] < polygon[j].s[1-k]) {
          const float w =
              // (p-polygon[i]).norm() / (polygon[i]-polygon[j]).norm();
              length(p-polygon[i]) / length(polygon[i]-polygon[j]);
          return (1-w)*color[i] + w*color[j];
        }
      }
    }
  }

  vector<float> alpha(n);
  vector<float> w(n);
  for (int i = 0; i < n; ++i) {
    const int j = (i+1)%n;    
    float2 A = polygon[i]-p;
    float2 B = polygon[j]-p;
    // if (p == A) return color[i];
    // if (p == B) return color[j];
    // alpha[i] = acos((A*B)/(A.norm()*B.norm()));
    alpha[i] = acos((dot(A,B))/(length(A)*length(B)));
  }
  float sum = 0;
  for (int i = 0; i < n; ++i) {
    const int h = (i+n-1)%n;
    // w[i] = (tan(alpha[h]/2)+tan(alpha[i]/2))/(polygon[i]-p).norm();
    w[i] = (tan(alpha[h]/2)+tan(alpha[i]/2))/length(polygon[i]-p);
    sum += w[i];
  }
  T c;
  for (int i = 0; i < n; ++i) {
    w[i] /= sum;
    c += (w[i] * color[i]);
  }
  return c;
}

string glErrorString(const int error) {
  switch (error) {
    case GL_NO_ERROR:
      return "GL_NO_ERROR";
      break;
    case GL_INVALID_ENUM:
      return "GL_INVALID_ENUM";
      break;
    case GL_INVALID_VALUE:
      return "GL_INVALID_VALUE";
      break;
    case GL_INVALID_OPERATION:
      return "GL_INVALID_OPERATION";
      break;
    case GL_INVALID_FRAMEBUFFER_OPERATION:
      return "GL_INVALID_FRAMEBUFFER_OPERATION";
      break;
    case GL_OUT_OF_MEMORY:
      return "GL_OUT_OF_MEMORY";
      break;
    case GL_STACK_UNDERFLOW:
      return "GL_STACK_UNDERFLOW";
      break;
    case GL_STACK_OVERFLOW:
      return "GL_STACK_OVERFLOW";
      break;
  }
}

void CheckGlError(const string& str) {
  GLenum err = glGetError();
  while (err != GL_NO_ERROR) {
    cerr << "OpenGL error " << str.c_str() <<
                    " : " << glErrorString(err);
    err = glGetError();
  }
}

void CheckGlError(const int i) {
  GLenum err = glGetError();
  while (err != GL_NO_ERROR) {
    cerr << "OpenGL error " << i <<
                    " : " << glErrorString(err);
    err = glGetError();
  }
}

struct M2Callback {
  M2Callback(const OctViewer2* m2_) : m2(m2_) {}
  const OctViewer2* m2;
};

struct M2CallbackNoConst {
  M2CallbackNoConst(OctViewer2* m2_) : m2(m2_) {}
  OctViewer2* m2;
};

struct StoreVertexLocationCallback : public M2CallbackNoConst {
  StoreVertexLocationCallback(OctViewer2* m2_) : M2CallbackNoConst(m2_) {}
  bool operator()(const int vi, const int2& p) {
    return m2->StoreVertexLocation(vi, p);
  }
};

struct DrawVertexDistanceLineCallback : public M2Callback {
  DrawVertexDistanceLineCallback(const OctViewer2* m2_) : M2Callback(m2_) {}
  bool operator()(const int vi, const int2& p) {
    return m2->DrawVertexDistanceLine(vi, p);
  }
};

struct DrawGVDLinesCallback : public M2Callback {
  typedef oct::OctreeOptions OctreeOptions;
  DrawGVDLinesCallback(const OctViewer2* m2_, const OctreeOptions& o_)
      : M2Callback(m2_),
        graph_constructor(new GraphConstructor<2>()),
        o(o_) {}
  bool operator()(const int vi, const int2& p) {
    const oct::VertexNetwork& vertices = m2->Vertices();
    // return m2->DrawGVDLines(vi, p);
    if (!vertices.IsBase(vi)) return true;

    vector<pair<int2, LabeledSegment> > intersections;
    oct::GetIntersectionsAroundFace<2>(vi, p, vertices.CellLevel(vi),
                                       0, 1,
                                       vertices,
                                       back_inserter(intersections), o);

    if (!intersections.empty()) {
      // gvd_cells.push_back(GVDCell(vi, intersections));
      double2 center = make_double2(0);
      for (int i = 0; i < intersections.size(); ++i) {
        const int2& edge_intersection = intersections[i].first;
        center += convert_double2(m2->Oct2Obj(edge_intersection));
      }
      center /= intersections.size();
      for (int i = 0; i < intersections.size(); ++i) {
        const int2& edge_intersection = intersections[i].first;
        const Edge e = make_edge(intersections[i].second[0].Index(),
                     intersections[i].second[1].Index());
        graph_constructor->AddEdge(convert_double2(m2->Oct2Obj(edge_intersection)),
                                  center, e);
      }
    }
    return true;
  }
  const Graph<2>& GetGraph() const { return graph_constructor->GetGraph(); }
  oct::shared_ptr<GraphConstructor<2> > graph_constructor;
  OctreeOptions o;
};

struct DrawEdgeCallback : public M2Callback {
  DrawEdgeCallback(const OctViewer2* m2_) : M2Callback(m2_) {}
  bool operator()(const int vi, const int n_vi,
                  const int2& p, const int2& q,
                  // const oct::Direction<2>& d) {
                  const oct::Direction& d) {
    return m2->DrawEdge(vi, n_vi, p, q, d);
  }
};

struct VertexIDCallback : public M2Callback {
  VertexIDCallback(const OctViewer2* m2_) : M2Callback(m2_) {}
  bool operator()(const int vi, const int2& p) {
    return m2->DrawVertexID(vi, p);
  }
};

struct VertexLabelCallback : public M2Callback {
  VertexLabelCallback(const OctViewer2* m2_) : M2Callback(m2_) {}
  bool operator()(const int vi, const int2& p) {
    return m2->DrawVertexLabel(vi, p);
  }
};

struct DistanceFunctionCallback : public M2Callback {
  DistanceFunctionCallback(const OctViewer2* m2_) : M2Callback(m2_) {}
  bool operator()(const int vi, const int2& p) {
    return m2->DistanceFunctionVertex(vi, p);
  }
};

// Finds the max distance from a vertex to a site.
struct DFDistCallback : public M2Callback {
  DFDistCallback(const OctViewer2* m2_,
                 float* dist_oct_,
                 float* dist_obj_)
      : M2Callback(m2_), dist_oct(dist_oct_), dist_obj(dist_obj_) {}
  bool operator()(const int vi, const int2& p_oct) {
    const int2 q_oct = m2->Vertices().ClosestPoint(vi);
    const float d_oct = length(q_oct - p_oct);
    const float2 q_obj = m2->Oct2Obj(q_oct);
    const float d_obj = length(q_obj - m2->Oct2Obj(p_oct));
    if (d_oct > *dist_oct) {
      *dist_oct = d_oct;
      *dist_obj = d_obj;
    }
    return true;
  }
  float* dist_oct;
  float* dist_obj;
};

const float3 OctViewer2::red = make_float3(1, 0, 0);

// OctViewer2::OctViewer2(const int2& win, const float2& world_min,
//                  const float2& world_max) : GL2D(win, world_min, world_max) {
OctViewer2::OctViewer2(const int win_width, const int win_height)
    : GL2D(win_width, win_height) {
  mouse_active = false;

  dirty = false;

  show_help = false;
  show_advanced_help = false;
  show_octree = true;
  show_vertex_labels = false;
  show_cell_vertices = false;
  show_vertex_distance = false;
  polygon_mode = 1;
  show_closest_point_line = false;
  show_gvd = true;
  show_path = true;
  show_voronoi = false;
  show_poly_vertices = 0;
  show_vertex_ids = false;
  show_distance_field = 0;
  show_statistics = true;
  max_vertex_distance = oct::Level2CellWidth(0);
  max_dist_oct = 0;
  max_dist_obj = 0;
  o = oct::OctreeOptions::For2D();
  entry_mode = 0;
  octree_color = make_double3(0.7, 0.7, 0.7);
  seg_a = seg_b = make_float2(-1, -1);
}

void OctViewer2::ReadMesh(const string& filename) {
  ifstream in(filename.c_str());
  if (!in) {
    cerr << "Failed to read " << filename << endl;
    return;
  }
  double x, y;
  while (!in.eof()) {
    in >> x >> y;
    if (!in.eof()) {
      float2 v = make_float2(x, y);
      verts.push_back(v);
      bb(v);
      // bb = bb.CenteredSquare();
    }
    dirty = true;
  }
  in.close();
  MakePolygon();
}

void OctViewer2::WritePolygons() const {
  for (int i = 0; i < polygons.size(); ++i) {
    stringstream ss;
    ss << "poly-" << (i+1) << ".dat";
    ofstream out(ss.str().c_str());
    const vector<float2>& polygon = polygons[i];
    for (int j = 0; j < polygon.size(); ++j) {
      out << std::setprecision(12) << polygon[j].s[0] << " "
          << std::setprecision(12) << polygon[j].s[1] << endl;
    }
    out.close();
  }
  cout << "Wrote objects" << endl;
}

class WriteGVDVisitor {
 public:
  WriteGVDVisitor(std::ostream& out) : _out(out) {}

  void operator()(int ai, const double2& a,
                  int bi, const double2& b) {
    _out << a << " " << b << endl;
  }
 private:
  std::ostream& _out;
};

void OctViewer2::WriteGvdMesh() const {
  ofstream out("gvd.dat");
  gvd_graph.VisitEdges(WriteGVDVisitor(out));
  cout << "Wrote mesh" << endl;
}

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
  cout << "Usage: ./gvd-viewer2 [-l maxlevel] [filenames.dat]" << endl;
  cout << "  -l maxlevel        - Limit octree to maxlevel+1 levels" << endl;
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

void OctViewer2::test(const int test) {
  using namespace Karras;

  const bool oct_test = (test >= 0 && test <= 5);

  if (oct_test) {
    karras_points.clear();
    vector<intn> qpoints;
    Resln resln;
    if (test == 0) {
      resln = Resln(8);
      qpoints.push_back(z2xyz(1, resln));
      qpoints.push_back(z2xyz(2, resln));
      qpoints.push_back(z2xyz(4, resln));
      qpoints.push_back(z2xyz(5, resln));
      qpoints.push_back(z2xyz(19, resln));
      qpoints.push_back(z2xyz(24, resln));
      qpoints.push_back(z2xyz(25, resln));
      qpoints.push_back(z2xyz(30, resln));
    } else if (test == 1) {
      resln = Resln(4);
      qpoints.push_back(z2xyz(0, resln));
      qpoints.push_back(z2xyz(2, resln));
      qpoints.push_back(z2xyz(7, resln));
      qpoints.push_back(z2xyz(15, resln));
      qpoints.push_back(z2xyz(15, resln));
    } else if (test == 2) {
      resln = Resln(8);
      srand(0);
      for (int i = 0; i < 10000; ++i) {
        qpoints.push_back(z2xyz(rand()%resln.volume, resln));
      }
    } else if (test == 3) {
      resln = Resln(128);
      srand(0);
      for (int i = 0; i < 10000; ++i) {
        qpoints.push_back(z2xyz(rand()%resln.volume, resln));
      }
    } else if (test == 4) {
      resln = Resln(1<<8);
      karras_points.push_back(make_floatn(0.204444, 0.0488889));
      karras_points.push_back(make_floatn(-0.417778, -0.177778));
      karras_points.push_back(make_floatn(-0.191111, -0.333333));
      karras_points.push_back(make_floatn(-0.364444, -0.417778));
      karras_points.push_back(make_floatn(-0.155556, -0.475556));
      qpoints = Quantize(karras_points, resln);
    } else if (test == 5) {
      resln = Resln(1<<8);
      karras_points.push_back(make_floatn(0.231111, 0.0977777));
      karras_points.push_back(make_floatn(-0.422222, -0.288889));
      karras_points.push_back(make_floatn(-0.288889, -0.142222));
      karras_points.push_back(make_floatn(-0.337778, 0.12));
      karras_points.push_back(make_floatn(-0.0311111, 0.155556));
      qpoints = Quantize(karras_points, resln);
    }

    cout << "Resln = (" << resln << ")" << endl;
    vector<Morton> tmp;
    for (const intn p : qpoints) {
      cout << p << endl;
      tmp.push_back(xyz2z(p, resln));
    }
    sort(tmp.begin(), tmp.end());
    for (const Morton mp : tmp) {
      // const std::bitset<16> bmp(mp);
      // cout << std::setw(5) << mp << ", " << bmp << endl;
    }

    for (int i = 0; i < 1000; ++i) {
      octree = Karras::BuildOctree(qpoints, resln, false);
      if (octree.empty()) {
        cout << "failure on iteration " << i << endl;
        exit(0);
      }
    }

    if (karras_points.empty()) {
      // We only tested with quantized points so it makes no sense to display
      exit(0);
    }
  
    bb = BoundingBox<float2>();
    for (const floatn& p : karras_points) {
      bb(p);
    }
  } else if (o.test == 6) {
    // not oct_test
    using namespace Karras;

    intn p;
    vector<OctNode> octree;
    Resln resln;

    ifstream in("find.err");
    in >> p;
    in >> resln;
    in >> octree;
    in.close();

    cout << "p = " << p << endl;
    cout << "resln = " << resln << endl;
    cout << "octree = " << octree << endl;

    FindLeaf(p, octree, resln);
  
    exit(0);
  } else if (o.test == 7) {
    using namespace Karras;

    OctCell cell;
    floatn a, b;
    vector<OctNode> octree;
    Resln resln;

    ifstream in("mccallback.err");
    in >> cell >> a >> b >> octree >> resln;
    in.close();

    cerr << "a = " << a << endl;
    cerr << "b = " << b << endl;
    cerr << "cell = " << cell << endl;
    cerr << "resln = " << resln << endl;
    const vector<CellIntersection> intersections =
        FindIntersections(a, b, cell, octree, resln);

    if (intersections.size() > 2) {
      cerr << "More than 2 intersections in cell" << endl;
    }
    cout << "Number of intersections: " << intersections.size() << endl;
  
    exit(0);
  } else if (o.test == 8) {
    float_seg aa(make_floatn(0, 1.2), make_floatn(1, 10));
    float_seg bb(make_floatn(0, 0), make_floatn(10, 10));
    vector<floatn> samples = v_sample(aa, bb);
    for (const floatn sample : samples) {
      cout << sample << endl;
    }
    exit(0);
  } else if (o.test == 9) {
    float_seg a(make_floatn(160, 65), make_floatn(155, 65));
    float_seg b(make_floatn(160, 67), make_floatn(156, 67));
    vector<floatn> samples;
    vector<floatn> origins;
    vector<float> lengths;
    FitBoxesNoIntersection(a, b, 1, &samples, &origins, &lengths);
    exit(0);
  }
}

int OctViewer2::ProcessArgs(int argc, char** argv) {
  int i = 1;
  // if (argc > 1) {
  bool stop = false;
  while (i < argc && !stop) {
    stop = true;
    if (o.ProcessArg(i, argv)) {
      stop = false;
    }
  }
  if (o.help) {
    PrintUsage();
    exit(0);
  }

  resln = Karras::Resln(1<<o.max_level);

  if (o.test > -1)
    test(o.test);

  for (; i < argc; ++i) {
    string filename(argv[i]);
    cout << filename << endl;
    o.filenames.push_back(filename);
    // ReadMesh(filename);
  }

  if (o.filenames.empty()) {
    throw logic_error("No filenames");
  }

  for (const string& f : o.filenames) {
    ReadMesh(f);
  }

  int num_edges = 0;
  for (int i = 0; i < polygons.size(); ++i) {
    num_edges += polygons[i].size();
  }

  PrintCommands();
  cout << "Number of objects: " << polygons.size() << endl;
  cout << "Number of polygon edges: " << num_edges << endl;

  return 0;
}

void OctViewer2::PrintCommands() const {
}

void OctViewer2::Keyboard(unsigned char key, int x, int y) {
  bool redisplay = true;
  // cout << glutGetModifiers() << endl;
  // cout << key << endl;
  switch (key) {
    case 'h':
      show_help = !show_help;
      show_advanced_help = false;
      break;
    case 'H':
      show_advanced_help = !show_advanced_help;
      show_help = false;
      break;
    case 'o':
      if (glutGetModifiers() & GLUT_ACTIVE_ALT) {
        for (int i = 0; i < 3; ++i) {
          octree_color.s[i] = std::min(octree_color.s[i] * 1.1, 1.0);
        }
      } else {
        show_octree = !show_octree;
      }
      break;
    case 'O':
      if (glutGetModifiers() & GLUT_ACTIVE_ALT) {
        for (int i = 0; i < 3; ++i) {
          octree_color.s[i] = octree_color.s[i] / 1.1;
        }
      }
      break;
    case 'v':
      show_vertex_labels = !show_vertex_labels;
      break;
    case 'a':
      show_poly_vertices = (++show_poly_vertices) % 4;
      break;
    case 'l':
      show_closest_point_line = !show_closest_point_line;
      break;
    case 'm':
      show_gvd = !show_gvd;
      break;
    case 'P':
      show_path = !show_path;
      break;
    case 'g':
      polygon_mode = (polygon_mode+1)%3;
      break;
    case 'r': {
      // show_voronoi = !show_voronoi;
      const float w = window_width;
      const float h = window_height;
      const float r = w/h;
      win_obj = BB(make_float2(-r, -1), make_float2(r, 1));
      break;
    }
    case 'w':
      WriteGvdMesh();
      break;
    case 'W':
      WritePolygons();
      break;
    case 'i':
      show_vertex_ids = !show_vertex_ids;
      break;
    case 'F':
      if (o.max_level < oct::kMaxLevel) {
        ++o.max_level;
      }
      // o.max_level = min((oct::level_t)(o.max_level+1),
      //                   oct::kMaxLevel);
      cout << "Max octree level = " << static_cast<int>(o.max_level) << endl;
      dirty = true;
      break;
    case 'D':
      o.max_level = maximum(o.max_level-1, 0);
      cout << "Max octree level = " << static_cast<int>(o.max_level) << endl;
      dirty = true;
      break;
    case 'f':
      ++o.karras_iterations;
      cout << "Max Karras iterations = " << o.karras_iterations << endl;
      dirty = true;
      break;
    case 'd':
      o.karras_iterations = maximum(o.karras_iterations-1, 1);
      cout << "Max Karras iterations = " << o.karras_iterations << endl;
      dirty = true;
      break;
    case 's':
      o.ambiguous_max_level++;
      cout << "Max GVD subdivision level = "
           << o.ambiguous_max_level << endl;
      dirty = true;
      break;
    case 'S':
      o.ambiguous_max_level--;
      cout << "Max GVD subdivision level = "
           << o.ambiguous_max_level << endl;
      dirty = true;
      break;
    case 'j':
      max_vertex_distance *= 2;
      break;
    case 'k':
      max_vertex_distance /= 2;
      break;
    case 'c':
      verts.clear();
      polygons.clear();
      bb = BoundingBox<float2>();
      // vertices = VertexNetwork();
      vertices.Clear();
      break;
    case 'C':
      middle_down.clear();
      middle_up.clear();
      break;
    case 'u':
      show_distance_field = (show_distance_field+1)%4;
      break;
    case 'V':
      o.full_subdivide = !o.full_subdivide;
      dirty = true;
      break;
    case 'B':
      o.make_buffer = !o.make_buffer;
      dirty = true;
      break;
    case 'z': {
      if (!polygons.empty()) {
        polygons.resize(polygons.size()-1);
      }
      break;
    }
    case 'e':
      entry_mode = (entry_mode+1)%3;
      break;
    case 't': {
      // test
      const float2 win_min_obj = make_float2(-0.4, -0.4);
      const float2 win_max_obj = make_float2(0.4, 0.4);
      win_obj = BB(win_min_obj, win_max_obj);
      redisplay = true;
      break;
    }
    case 'p':
      show_statistics = !show_statistics;
      break;
    default:
      redisplay = false;
      break;
  }
  if (redisplay) {
    glutPostRedisplay();
  }
}

void OctViewer2::Special(int key, int x, int y) {
  // switch (key) {
  //   case GLUT_KEY_DOWN:
  //     maxDepth--;
  //     cout << "octree depth = " << maxDepth << endl;
  //     dirty = true;
  //     break;
  //   case GLUT_KEY_UP:
  //     maxDepth++;
  //     cout << "octree depth = " << maxDepth << endl;
  //     dirty = true;
  //     break;
  // }
  // glutPostRedisplay();
}

// x, y in window coordinates
void OctViewer2::AddPoint(int x, int y) {
  float2 v = Win2Obj(make_float2(x, y));
  verts.push_back(v);
  bb(v);
  // bb = bb.CenteredSquare();
  dirty = true;
}

void OctViewer2::Search(const int start, const int end) {
  search_path.clear();
  gvd_graph.Dijkstra(start, end, back_inserter(search_path));
  reverse(search_path.begin(), search_path.end());
  glutPostRedisplay();
}

void OctViewer2::SetStartSearch(int x, int y) {
  const double2 p = convert_double2(Win2Obj(make_float2(x, y)));
  double min_dist = numeric_limits<double>::maximum();
  int start = -1;
  const std::vector<double2>& vertices = gvd_graph.GetVertices();
  for (int i = 0; i < vertices.size(); ++i) {
    // const double d = (p-vertices[i]).norm2();
    const double d = length2(p-vertices[i]);
    if (d < min_dist) {
      min_dist = d;
      start = i;
    }
  }
  if (start > -1) {
    // search_path.clear();
    // search_path.push_back(start);
    if (search_path.empty()) {
      search_path.push_back(start);
    } else {
      Search(start, search_path.back());
    }
  }
}

void OctViewer2::SetEndSearch(int x, int y) {
  const double2 p = convert_double2(Win2Obj(make_float2(x, y)));
  double min_dist = numeric_limits<double>::maximum();
  int end = -1;
  const std::vector<double2>& vertices = gvd_graph.GetVertices();
  for (int i = 0; i < vertices.size(); ++i) {
    // const double d = (p-vertices[i]).norm2();
    const double d = length2(p-vertices[i]);
    if (d < min_dist) {
      min_dist = d;
      end = i;
    }
  }
  if (end > -1) {
    // if (!search_path.empty()) {
    //   Search(search_path.front(), end);
    // }
    if (search_path.empty()) {
      search_path.push_back(end);
    } else {
      Search(search_path.front(), end);
    }
  }
}

void OctViewer2::Zoom(const float2& target, const float zoom) {
  const float w = window_width;
  const float h = window_height;
  const float r = w/h;
  const BB full(make_float2(-r, -1), make_float2(r, 1));
  const BB zoomed = full * (1/zoom);
  BB centered;
  centered(target-zoomed.size()/2);
  centered(target+zoomed.size()/2);
  win_obj = centered;
}

void OctViewer2::Find(int x, int y) {
  using namespace Karras;

  floatn fv = Obj2Oct(Win2Obj(make_floatn(x, y)));
  intn v = make_intn(fv.s[0], fv.s[1]);
  fnode = Karras::FindLeaf(v, octree, resln);
  dirty = true;
  glutPostRedisplay();

  // seg_a = seg_b;
  // // seg_b = v;
  // seg_b = fv;
  
  // // if (seg_a == make_int2(-1, -1))
  // if (seg_a == make_float2(-1, -1))
  //   return;

  // const vector<CellIntersection> all_intersections = Walk(seg_a, seg_b);
  // intersections.clear();
  // for (const CellIntersection& i : all_intersections) {
  //   // const intn p = i.p;
  //   const floatn p = i.p;
  //   intersections.push_back(p);
  // }
}

// Cell walk visitor
// typedef void (*cwv)(Karras::OctCell, const intn& a, const intn& b,
typedef void (*cwv)(Karras::OctCell, const floatn& a, const floatn& b,
    const vector<OctNode>& octree, const Karras::Resln& resln, void* data);

void CellWalk(
    // const intn& a, const intn& b,
    const floatn& a, const floatn& b,
    const vector<OctNode>& octree, const Karras::Resln& resln,
    cwv v, void* data) {
  using namespace Karras;

  int dir[DIM];
  for (int i = 0; i < DIM; ++i) {
    if (b.s[i] > a.s[i]) dir[i] = 1;
    else if (b.s[i] < a.s[i]) dir[i] = -1;
    else dir[i] = 0;
  }
  CellIntersection cur(0, a);
  const Karras::OctCell bcell = Karras::FindLeaf(convert_intn(b),
                                                 octree, resln);
  vector<CellIntersection> local_intersections;
  vector<CellIntersection> all_intersections;
  bool done = false;
  int count = 0;
  Karras::OctCell cell = Karras::FindLeaf(convert_intn(cur.p), octree, resln);
  do {
    ++count;
    if (count > 10000)
      throw logic_error("Infinite loop");
    v(cell, a, b, octree, resln, data);

    const vector<CellIntersection> local_intersections =
        FindIntersections(a, b, cell, octree, resln);
    for (const CellIntersection& i : local_intersections) {
      if (i.t > cur.t) {
        cur = i;
      }
    }
    for (int i = 0; i < DIM; ++i) {
      const int end = cell.get_origin().s[i];
      if (cur.p.s[i] == end && cur.p.s[i] != 0 && dir[i] == -1) {
        --cur.p.s[i];
      }
    }

    Karras::OctCell new_cell = Karras::FindLeaf(
        convert_intn(cur.p), octree, resln);
    done = local_intersections.empty() ||
        (cell.get_origin() == bcell.get_origin()) ||
        (new_cell.get_origin() == cell.get_origin());
    cell = new_cell;
  } while (!done);
}

struct MCData {
  MCData(vector<CellIntersections>& labels_) : labels(labels_) {}
  vector<CellIntersections>& labels;
  int cur_label;
};

void write_seg(const floatn& a, const floatn& b, const string& fn) {
  ofstream out(fn);
  out << std::fixed << a << endl << b << endl;
  out.close();
}

void write_seg(const float_seg& s, const string& fn) {
  write_seg(s.a(), s.b(), fn);
}

void write_cell(const intn& o, const int& w, const string& fn) {
  ofstream out(fn);
  out << o << endl;
  out << o+make_intn(w, 0) << endl;
  out << o+make_intn(w, w) << endl;
  out << o+make_intn(0, w) << endl;
  out << o << endl;
  out.close();
}

void error_log(Karras::OctCell cell, const floatn& a, const floatn& b,
               const vector<OctNode>& octree, const Karras::Resln& resln) {
  cerr << "a = " << std::fixed << a << endl;
  cerr << "b = " << std::fixed << b << endl;
  cerr << "cell = " << cell << endl;
  cerr << "resln = " << resln << endl;
  ofstream out("mccallback.err");
  out << cell << " " << a << " " << b
      << " " << octree << " " << resln << endl;
  out.close();
  write_seg(a, b, "data1.dat");
  const intn o = cell.get_origin();
  const int w = cell.get_width();
  write_cell(o, w, "data2.dat");
  const int w2 = w/2;
  // Boundary
  write_cell(o-make_intn(w2, w2), 2*w, "data3.dat");
}

inline bool within(const float& a, const float& b, const float& epsilon) {
  return fabs(a-b) < epsilon;
}

// void MCCallback(Karras::OctCell cell, const intn& a, const intn& b,
void MCCallback(Karras::OctCell cell, const floatn& a, const floatn& b,
                  const vector<OctNode>& octree, const Karras::Resln& resln,
                  void* data) {
  // const static float EPSILON = 1e-6;

  using namespace Karras;
  MCData* d = static_cast<MCData*>(data);
  const int octant = cell.get_octant();

  // debug
  OctCell test_cell;
  stringstream ss("8796416 6561024 256 22009 3 -1");
  ss >> test_cell;
  if (cell == test_cell || cell == fnode) {
    cout << "here in MCCallback" << endl;
  }

  const vector<CellIntersection> intersections =
      FindIntersections(a, b, cell, octree, resln);
  if (intersections.size() > 2) {
    error_log(cell, a, b, octree, resln);
    cerr << "More than 2 intersections with a cell" << endl;
    for (const CellIntersection& ci : intersections) {
      cerr << "  " << ci.p << endl;
    }
    throw logic_error("More than 2 intersections with a cell");
  }
  if (intersections.size() == 2) {
    float_seg seg(intersections[0].p,
                  intersections[1].p);
    CellIntersections& cell_intersections = d->labels[cell.get_parent_idx()];
    cell_intersections.set(octant, d->cur_label, seg);
  } else if (intersections.size() == 1) {
    // Find which endpoint is in cell
    const BoundingBox<intn> bb = cell.bb();
    const floatn p = intersections[0].p;
    floatn q;
    bool degenerate = false;
    if (bb.in_half_open(convert_intn(a), EPSILON)) {
      q = a;
    } else if (bb.in_half_open(convert_intn(b), EPSILON)) {
      q = b;
    } else if ((within(p.x, bb.min().x, EPSILON) ||
                within(p.x, bb.max().x, EPSILON)) &&
               (within(p.y, bb.min().y, EPSILON) ||
                within(p.y, bb.max().y, EPSILON))) {
      degenerate = true;
    } else {
      degenerate = true;
      // error_log(cell, a, b, octree, resln);
      // throw logic_error("Only one intersection but neither endpoint is"
      //                   " in bounding box");
    }
    if (!degenerate) {
      float_seg seg(p, q);
      CellIntersections& cell_intersections = d->labels[cell.get_parent_idx()];
      cell_intersections.set(octant, d->cur_label, seg);
    }
  }
}

void OctViewer2::FindMultiCells() {
  using namespace Karras;

  cell_intersections.clear();
  cell_intersections.resize(octree.size(), CellIntersections());
  MCData data(cell_intersections);

  vector<vector<float2> > temp_polygons(polygons);
  if (!verts.empty()) {
    temp_polygons.push_back(verts);
  }

  // For each line segment in each polygon do a cell walk, updating
  // the cell_intersections structure.
  intersections.clear();
  for (int j = 0; j < temp_polygons.size(); ++j) {
    const std::vector<float2> polygon = temp_polygons[j];
    data.cur_label = j;
    for (int i = 0; i < polygon.size() - 1; ++i) {
      const floatn a = Obj2Oct(polygon[i]);
      const floatn b = Obj2Oct(polygon[i+1]);
      CellWalk(a, b, octree, resln, MCCallback, &data);
      vector<CellIntersection> local = Walk(a, b);
      for (const CellIntersection& ci : local) {
        intersections.push_back(ci.p);
      }
    }
  }

  // Find points that we want to add in order to run a more effective
  // Karras octree construction.
  extra_qpoints.clear();
  _origins.clear();
  _lengths.clear();
  for (int i = 0; i < octree.size(); ++i) {
    const CellIntersections& intersection = cell_intersections[i];
    for (int octant = 0; octant < (1<<DIM); ++octant) {
      const int num_labels = intersection.num_labels(octant);
      // cout << "cell = " << i << "/" << octant << ": " << num_labels << endl;
      // if (intersection.is_multi(octant)) {
      //   {
      for (int j = 0; j < num_labels; ++j) {
      // for (int j = 0; j < 1; ++j) {
        for (int k = j+1; k < num_labels; ++k) {
        // for (int k = j+1; k < std::min(2, num_labels); ++k) {
          // We'll compare segments at j and k.
          float_seg segs[2] = { intersection.seg(j, octant),
                                intersection.seg(k, octant) };
          vector<floatn> samples;
          vector<floatn> origins;
          vector<float> lengths;
          if (!multi_intersection(segs[0], segs[1])) {
            try {
              if (i == fnode.get_parent_idx() && octant == fnode.get_octant()) {
                cout << "here in FitBoxes" << endl;
                write_seg(segs[0], "seg0.dat");
                write_seg(segs[1], "seg1.dat");
              }
              FitBoxes(segs[0], segs[1], 1, &samples, &origins, &lengths);
            } catch(logic_error& e) {
              cerr << "segments: " << segs[0] << " " << segs[1] << endl;
              cerr << "labels: " << intersection.label(0, octant)
                   << " " << intersection.label(1, octant) << endl;
              write_seg(segs[0], "seg0.dat");
              write_seg(segs[1], "seg1.dat");
              throw e;
            }
          }
          _origins.insert(_origins.end(), origins.begin(), origins.end());
          _lengths.insert(_lengths.end(), lengths.begin(), lengths.end());
          for (const floatn& sample : samples) {
            extra_qpoints.push_back(convert_intn(sample));
          }
        }
      }
    }
  }
}

// void WalkCallback(Karras::OctCell cell, const intn& a, const intn& b,
void WalkCallback(Karras::OctCell cell, const floatn& a, const floatn& b,
                  const vector<OctNode>& octree, const Karras::Resln& resln,
                  void* data) {
  
  using namespace Karras;

  vector<CellIntersection>& all_intersections =
      *static_cast<vector<CellIntersection>*>(data);

  const vector<CellIntersection> local_intersections =
      FindIntersections(a, b, cell, octree, resln);
  for (const CellIntersection& i : local_intersections) {
    // const intn p = i.p;
    all_intersections.push_back(i);
  }
}

vector<Karras::CellIntersection> OctViewer2::Walk(
    // const intn& a, const intn& b) {
    const floatn& a, const floatn& b) {
  using namespace Karras;

  vector<CellIntersection> all_intersections;
  CellWalk(a, b, octree, resln, WalkCallback, &all_intersections);

  return all_intersections;
}

// vector<Karras::CellIntersection> OctViewer2::Walk(
//     const intn& a, const intn& b, vector<Karras::OctCell>* cells) {
//   using namespace Karras;

//   int dir[DIM];
//   for (int i = 0; i < DIM; ++i) {
//     if (b.s[i] > a.s[i]) dir[i] = 1;
//     else if (b.s[i] < a.s[i]) dir[i] = -1;
//     else dir[i] = 0;
//   }
//   CellIntersection cur(0, a);
//   const Karras::OctCell acell = Karras::FindLeaf(a, octree, resln);
//   const Karras::OctCell bcell = Karras::FindLeaf(b, octree, resln);
//   if (cells)
//     cells->push_back(acell);
//   vector<CellIntersection> local_intersections;
//   vector<CellIntersection> all_intersections;
//   // intersections.clear();
//   bool done = false;
//   // Currently does not work for a vertical line?
//   int count = 0;
//   Karras::OctCell cell = Karras::FindLeaf(cur.p, octree, resln);
//   do {
//     ++count;
//     if (count > 10000) throw logic_error("Infinite loop");
//     local_intersections = FindIntersections(a, b, cell, octree, resln);
//     for (const CellIntersection& i : local_intersections) {
//       const intn p = i.p;
//       // intersections.push_back(p);
//       all_intersections.push_back(i);
//       if (i.t > cur.t) {
//         cur = i;
//       }
//     }
//     for (int i = 0; i < DIM; ++i) {
//       const int end = cell.get_origin().s[i];
//       if (cur.p.s[i] == end && dir[i] == -1) {
//         --cur.p.s[i];
//       }
//     }
//     Karras::OctCell new_cell = Karras::FindLeaf(cur.p, octree, resln);
//     done = local_intersections.empty() ||
//         (cell.get_origin() == bcell.get_origin()) ||
//         (new_cell.get_origin() == cell.get_origin());
//     cell = new_cell;
//   } while (!done);

//   return all_intersections;
// }

float2 rubberband_start = make_float2(0);
float2 rubberband_cur = make_float2(0);
bool rubberband_mode = false;
void OctViewer2::Mouse(int button, int state, int x, int y) {
  if (button == GLUT_LEFT_BUTTON) {
    if (state == GLUT_DOWN) {
      if (glutGetModifiers() == GLUT_ACTIVE_CTRL) {
        rubberband_start = Win2Obj(make_float2(x, y));
        rubberband_mode = true;
      // } else if (glutGetModifiers() == GLUT_ACTIVE_SHIFT) {
      //   SetStartSearch(x, y);
      } else if (entry_mode == 2 && glutGetModifiers() == GLUT_ACTIVE_SHIFT) {
        AddPoint(x, y);
        verts.push_back(verts[0]);
        MakePolygon();
        dirty = true;
        glutPostRedisplay();
      } else if (glutGetModifiers() == GLUT_ACTIVE_SHIFT) {
        Find(x, y);
        glutPostRedisplay();
      } else {
        AddPoint(x, y);
        glutPostRedisplay();
        if (entry_mode < 2) {
          mouse_active = true;
        }
      }
    } else if (state == GLUT_UP) {
      if (glutGetModifiers() == GLUT_ACTIVE_CTRL) {
        rubberband_mode = false;
        const float w = window_width;
        const float h = window_height;
        const float r = w/h;
        const float2 rubberband_end = Win2Obj(make_float2(x, y));
        win_obj = BB();
        win_obj(rubberband_start);
        win_obj(rubberband_end);
        const float bbw = win_obj.size().s[0];
        const float bbh = win_obj.size().s[1];
        const float bbr = bbw / bbh;
        if (bbr > r) {
          float2 new_pt =
              make_float2(win_obj.min().s[0], win_obj.max().s[1] + (h/w)*bbw - bbh);
          win_obj(new_pt);
        } else {
          float2 new_pt =
              make_float2(win_obj.max().s[0] + (w/h)*bbh - bbw, win_obj.min().s[1]);
          win_obj(new_pt);
        }
        rubberband_start = float2();
        rubberband_cur = float2();
        glutPostRedisplay();
      } else {
        if (entry_mode < 2 && mouse_active) {
          mouse_active = false;
          // todo replace this when moving back to polygons
          // if (entry_mode == 0)
          //   verts.push_back(verts[0]);
          MakePolygon();
          // todo replace this when moving back to polygons
          // dirty = true;
          // glutPostRedisplay();
        }
      }
    }
  } else if (button == GLUT_RIGHT_BUTTON) {
    if (state == GLUT_DOWN) {
      if (glutGetModifiers() == GLUT_ACTIVE_SHIFT) {
        SetEndSearch(x, y);
      } else if (entry_mode == 2) {
        AddPoint(x, y);
        verts.push_back(verts[0]);
        MakePolygon();
        dirty = true;
        glutPostRedisplay();
      }
    }
  } else if (button == GLUT_MIDDLE_BUTTON) {
    if (state == GLUT_DOWN) {
      middle_down.push_back(Win2Obj(make_float2(x, y)));
      mouse_obj = Win2Obj(make_float2(x, y));
      glutPostRedisplay();
    } else {
      middle_up.push_back(Win2Obj(make_float2(x, y)));
    }
    glutPostRedisplay();
  }
}

void OctViewer2::MouseMotion(int x, int y) {
  // if (glutGetModifiers() == GLUT_ACTIVE_CTRL) {
  if (rubberband_mode) {
    rubberband_cur = Win2Obj(make_float2(x, y));
  } else if (mouse_active) {
    AddPoint(x, y);
  }
  mouse_obj = Win2Obj(make_float2(x, y));
  glutPostRedisplay();
}

void OctViewer2::PassiveMouseMotion(int x, int y) {
  mouse_obj = Win2Obj(make_float2(x, y));
  glutPostRedisplay();
}

float2 OctViewer2::Obj2Oct(const float2& v) const {
  const GLfloat bbw = bb.max_size();
  float2 oct = ((v-bb.min())/bbw) * (float)resln.width;
  for (int i = 0; i < DIM; ++i) {
    if (oct.s[i] >= resln.width)
      oct.s[i] = resln.width - 0.0001;
    if (oct.s[i] < 0)
      throw logic_error("Obj2Oct cannot have coord less than zero");
  }
  return oct;
}

float2 OctViewer2::Oct2Obj(const int2& v) const {
  const float2 vf = make_float2(v.s[0], v.s[1]);
  const GLfloat bbw = bb.max_size();
  // return (vf/kWidth)*bbw+bb.min();
  return (vf/resln.width)*bbw+bb.min();
}

GLfloat OctViewer2::Oct2Obj(int dist) const {
  const GLfloat bbw = bb.max_size();
  // const GLfloat ow = kWidth;
  const GLfloat ow = resln.width;
  return (dist/ow)*bbw;
}

unsigned char OctViewer2::outcode(const float2& v) const {
  unsigned char code = 0;
  if (v.s[0] < -1) {
    code = code | 0x01;
  } else if (v.s[0] > 1) {
    code = code | 0x04;
  }
  if (v.s[1] < -1) {
    code = code | 0x02;
  } else if (v.s[1] > 1) {
    code = code | 0x08;
  }
  return code;
}

bool OctViewer2::InBounds(const float2& v) const {
  return (outcode(v) == 0);
}

// Given a code, returns the corner that represents the
// intersection of the code's line and the next code's line.
float2 OctViewer2::corner(unsigned char code) const {
  if (code == 0x0001) {
    return make_float2(-1, -1);
  }
  if (code == 0x0002) {
    return make_float2(1, -1);
  }
  if (code == 0x0004) {
    return make_float2(1, 1);
  }
  if (code == 0x0008) {
    return make_float2(-1, 1);
  }
  throw logic_error("Unknown code");
}

// Clip a line to the bounds
// a is outside, b is inside
vector<float2> OctViewer2::clip(const float2& a0, const float2& a1,
                            const float2& b0, const float2& b1) const {
  // x = x0 + t*xd
  // y = y0 + t*yd
  // t = (x - x0)/xd
  // t = (y - y0)/yd
  const float ax0 = a0.s[0];
  const float ay0 = a0.s[1];
  const float axd = a1.s[0] - ax0;
  const float ayd = a1.s[1] - ay0;
  const float bx0 = b0.s[0];
  const float by0 = b0.s[1];
  const float bxd = b1.s[0] - bx0;
  const float byd = b1.s[1] - by0;
  unsigned char acode = outcode(a0);
  unsigned char bcode = outcode(b0);
  float at, bt;
  if (acode & 0x0001) {
    // clip to left
    at = (-1 - ax0) / axd;
  } else if (acode & 0x0002) {
    // clip to bottom
    at = (-1 - ay0) / ayd;
  } else if (acode & 0x0004) {
    // clip to right
    at = (1 - ax0) / axd;
  } else if (acode & 0x0008) {
    // clip to top
    at = (1 - ay0) / ayd;
  }
  if (bcode & 0x0001) {
    // clip to left
    bt = (-1 - bx0) / bxd;
  } else if (bcode & 0x0002) {
    // clip to bottom
    bt = (-1 - by0) / byd;
  } else if (bcode & 0x0004) {
    // clip to right
    bt = (1 - bx0) / bxd;
  } else if (bcode & 0x0008) {
    // clip to top
    bt = (1 - by0) / byd;
  }
  float2 anew = (a0 + (a1-a0) * at);
  float2 bnew = (b0 + (b1-b0) * bt);
  vector<float2> ret;
  ret.push_back(bnew);
  while (bcode != acode) {
    ret.push_back(corner(bcode));
    bcode *= 2;
    if (bcode > 0x0008) {
      bcode = 0x0001;
    }
  }
  ret.push_back(anew);

  return ret;
}

void OctViewer2::MakePolygon() {
  const int n = verts.size();
  polygons.push_back(verts);
  verts.clear();
}

bool OctViewer2::StoreVertexLocation(const int vi, const int2& p) {
  vertex_locations[vi] = Oct2Obj(p);
  return true;
}

void DefaultCallback(const int3& base_point,
                     const oct::index_t level,
                     const oct::index_t max_level,
                     const bool complete) {
}

void OctViewer2::BuildOctree() {
  //------------------
  // Initialize OpenCL
  //------------------
// #ifdef __OPEN_CL_SUPPORT__
//   static bool initialized = false;
//   if (o.gpu && !initialized) {
//     OpenCLInit(2, o, o.opencl_log);
//     initialized = true;
//   }
// #endif

  vector<vector<float2> > temp_polygons(polygons);
  if (!verts.empty()) {
    temp_polygons.push_back(verts);
  }
  // vector<float2> vvertices;
  karras_points.clear();
  vector<vector<float2> > all_vertices(temp_polygons.size());
  vector<vector<Edge> > all_edges(temp_polygons.size());
  for (int i = 0; i < temp_polygons.size(); ++i) {
    const vector<float2>& polygon = temp_polygons[i];
    all_vertices[i] = polygon;
    for (int j = 0; j < polygon.size()-1; ++j) {
      all_edges[i].push_back(make_edge(j, j+1));
      karras_points.push_back(temp_polygons[i][j]);
    }
    karras_points.push_back(temp_polygons[i].back());
  }

  // Find bounding box of vertices
  bb = BoundingBox<float2>();
  for (int j = 0; j < all_vertices.size(); ++j) {
    const std::vector<float2>& vertices = all_vertices[j];
    for (int i = 0; i < vertices.size(); ++i) {
      bb(vertices[i]);
    }
  }
  // bb = bb.CenteredSquare();

  // Karras iterations
  extra_qpoints.clear();
  vector<intn> qpoints = Karras::Quantize(karras_points, resln);
  int iterations = 0;
  do {
    qpoints.insert(qpoints.end(), extra_qpoints.begin(), extra_qpoints.end());
    for (const intn& qp : extra_qpoints) {
      karras_points.push_back(Oct2Obj(qp));
    }
    extra_qpoints.clear();
    if (qpoints.size() > 1) {
      octree = Karras::BuildOctree(qpoints, resln, false);
    } else {
      octree.clear();
    }
    FindMultiCells();

    ++iterations;
  } while (iterations < o.karras_iterations && !extra_qpoints.empty());
  cout << "Karras iterations: " << iterations << endl;

  // Count the number of cells with multiple intersections
  int count = 0;
  for (int i = 0; i < cell_intersections.size(); ++i) {
    for (int j = 0; j < 4; ++j) {
      if (cell_intersections[i].is_multi(j)) {
        ++count;
      }
    }
  }
  cout << "Number of multi-intersection cells: " << count << endl;
  
  //------------------
  // Cleanup OpenCL
  //------------------
// #ifdef __OPEN_CL_SUPPORT__
//   if (o.gpu) {
//     OpenCLCleanup();
//   }
// #endif
}

void OctViewer2::glSquare(GLfloat x, GLfloat y, GLfloat w) const {
  glBegin(GL_POLYGON);
  glVertex2f(x, y);
  glVertex2f(x+w, y);
  glVertex2f(x+w, y+w);
  glVertex2f(x, y+w);
  glEnd();
}

void OctViewer2::glSquare(const float2& p, GLfloat w) const {
  glSquare(p.s[0], p.s[1], w);
}

void OctViewer2::glSquareCentered(const float2& p, GLfloat w) const {
  glSquare(p.s[0]-w/2, p.s[1]-w/2, w);
}

void OctViewer2::glSquareWinWidth(const float2& p, GLfloat w) const {
  const float2 q = p - Win2Obj(w/2);
  glSquare(q, Win2Obj(w));
}

bool OctViewer2::DrawEdge(const int vi, const int n_vi,
                       const int2& p, const int2& q,
                       // const oct::Direction<2>& d) const {
                       const oct::Direction& d) const {
  glVertex2fv(Oct2Obj(p).s);
  glVertex2fv(Oct2Obj(q).s);
  return true;
}

void OctViewer2::DrawNode(
    const OctNode& parent, const int parent_idx,
    const intn origin, const int length) const {
  CheckGlError(1);

  for (int i = 0; i < 4; ++i) {
    intn o = origin;
    if (i % 2 == 1)
      o += make_intn(length/2, 0);
    if (i / 2 == 1)
      o += make_intn(0, length/2);
    if (!parent.is_leaf(i)) {
      DrawNode(octree[parent[i]], parent[i], o, length/2);
    } else {
      if (&parent == fnode.get_parent() && i == fnode.get_octant()) {
        glLineWidth(3);
        glColor3d(0.5, 0, 0.5);
        glSquare(Oct2Obj(o), Oct2Obj(length/2));
      }
      if (cell_intersections[parent_idx].is_multi(i)) {
        glLineWidth(5);
        glColor3d(1, 0, 0);
      } else {
        glLineWidth(1);
        glColor3dv(octree_color.s);
      }
      glSquare(Oct2Obj(o), Oct2Obj(length/2));
      if (show_vertex_ids) {
        stringstream ss;
        // ss << o;
        ss << parent_idx << "/" << i << " -- " << o;
        // ss << o << " " << w << " " << parent_idx << " " << i << " -1";
        BitmapString(ss.str(), Oct2Obj(o), 1, 3);
      }
    }
  }
}

void OctViewer2::DrawOctree() const {
  // cout << "DrawOctree " << octree.size() << endl;

  glLineWidth(1.0);
  glColor3dv(octree_color.s);

  if (!octree.empty()) {
    glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
    DrawNode(octree[0], 0, make_intn(0), resln.width);
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
  }

  // // Draw intersections
  // glColor3f(0.7, 0.0, 0.7);
  // for (const intn& i : intersections) {
  //   glSquareCentered(Oct2Obj(i), Win2Obj(7));
  // }

  // // Draw line segment
  // glColor3f(0.7, 0.0, 0.7);
  // glBegin(GL_LINES);
  // glVertex3fv(Oct2Obj(seg_a).s);
  // glVertex3fv(Oct2Obj(seg_b).s);
  // glEnd();
}

void DrawGVDVisitor(int ai, const double2& a,
                           int bi, const double2& b) {
  glVertex2dv(a.s);
  glVertex2dv(b.s);
}

void OctViewer2::DrawGVD() const {
  glColor3fv(red.s);
  glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
  glLineWidth(1.5);

  glBegin(GL_LINES);
  gvd_graph.VisitEdges(DrawGVDVisitor);
  glEnd();
}

void OctViewer2::DrawPath() const {
  glLineWidth(4.0);
  glColor3f(0.0, 0.0, 1.0);
  glBegin(GL_LINE_STRIP);
  for (int i = 0; i < search_path.size(); ++i) {
    glVertex2dv(gvd_graph[search_path[i]].s);
  }
  glEnd();
}

void OctViewer2::DrawVoronoi() const {
  // struct triangulateio in, mid, out, vorout;

  // set<float2> point_set;
  // gvd_cells.clear();
  // VisitVertices(vertices, DrawGVDLinesCallback(this));
  // // Update centers
  // for (int i = 0; i < gvd_cells.size(); ++i) {
  //   const GVDCell& mc = gvd_cells[i];
  //   if (mc.edge_intersections.size() > 2)
  //     point_set.insert(Oct2Obj(mc.center));
  //   for (int i = 0; i < mc.edge_intersections.size(); ++i) {
  //     point_set.insert(Oct2Obj(mc.edge_intersections[i]));
  //     glColor3f(0, 0, 1);
  //     glBegin(GL_LINES);
  //     // glVertex2fv(vertex_locations[mc.edges[i][0]]);
  //     // glVertex2fv(vertex_locations[mc.edges[i][1]]);
  //     glVertex2fv(Oct2Obj(mc.edges[i][0].Point()));
  //     glVertex2fv(Oct2Obj(mc.edges[i][1].Point()));
  //     glEnd();
  //   }
  // }

  // // const GLfloat w = Win2Obj(4);
  // // glColor3f(0, 1, 0);
  // // vector<float2> points = vertex_locations;
  // vector<float2> points(point_set.begin(), point_set.end());
  // if (points.size() < 3) return;
  // // for (int i = 0; i < vertex_locations.size(); ++i) {
  // //   const float2 p = vertex_locations[i];
  // //   // glSquare(p, w);
  // // }

  // // Define input points
  // in.numberofpoints = points.size();
  // in.numberofpointattributes = 1;
  // in.pointlist = new REAL[in.numberofpoints*2];
  // for (int i = 0; i < points.size(); ++i) {
  //   in.pointlist[i*2] = points[i][0];
  //   in.pointlist[i*2+1] = points[i][1];
  // }
  // in.pointattributelist =
  //     new REAL[in.numberofpoints*in.numberofpointattributes];
  // for (int i = 0; i < points.size(); ++i) {
  //   in.pointattributelist[i] = i;
  // }
  // in.pointmarkerlist = new int[in.numberofpoints];

  // in.numberofsegments = 0;
  // in.numberofholes = 0;
  // in.numberofregions = 1;
  // in.regionlist = new REAL[in.numberofregions*4];

  // // Return voronoi diagram in vorout

  // mid.pointlist = (REAL *) NULL;            /* Not needed if -N switch used. */
  // /* Not needed if -N switch used or number of point attributes is zero: */
  // mid.pointattributelist = (REAL *) NULL;
  // mid.pointmarkerlist = (int *) NULL; /* Not needed if -N or -B switch used. */
  // mid.trianglelist = (int *) NULL;          /* Not needed if -E switch used. */
  // /* Not needed if -E switch used or number of triangle attributes is zero: */
  // mid.triangleattributelist = (REAL *) NULL;
  // mid.neighborlist = (int *) NULL;         /* Needed only if -n switch used. */
  // /* Needed only if segments are output (-p or -c) and -P not used: */
  // mid.segmentlist = (int *) NULL;
  // /* Needed only if segments are output (-p or -c) and -P and -B not used: */
  // mid.segmentmarkerlist = (int *) NULL;
  // mid.edgelist = (int *) NULL;             /* Needed only if -e switch used. */
  // mid.edgemarkerlist = (int *) NULL;   /* Needed if -e used and -B not used. */

  // vorout.pointlist = (REAL *) NULL;        /* Needed only if -v switch used. */
  // /* Needed only if -v switch used and number of attributes is not zero: */
  // vorout.pointattributelist = (REAL *) NULL;
  // vorout.edgelist = (int *) NULL;          /* Needed only if -v switch used. */
  // vorout.normlist = (REAL *) NULL;         /* Needed only if -v switch used. */

  // /* Triangulate the points.  Switches are chosen to read and write a  */
  // /*   PSLG (p), preserve the convex hull (c), number everything from  */
  // /*   zero (z), assign a regional attribute to each element (A), and  */
  // /*   produce an edge list (e), a Voronoi diagram (v), and a triangle */
  // /*   neighbor list (n).                                              */

  // // triangulate("pczAevn", &in, &mid, &vorout);
  // char t_ops[] = "pczAevnQ";
  // triangulate(t_ops, &in, &mid, &vorout);

  // // printf("Initial triangulation:\n\n");
  // // report(&mid, 1, 1, 1, 1, 1, 0);
  // // printf("Initial Voronoi diagram:\n\n");
  // // report(&vorout, 0, 0, 0, 0, 1, 1);

  // triangulateio* io = &vorout;
  // const int n = io->numberofpoints;
  // vector<float2> vpoints(n);
  // vector<list<int> > edges(n);
  // for (int i = 0; i < n; i++) {
  //   vpoints[i] = float2(io->pointlist[i*2+0], io->pointlist[i*2+1]);
  // }
  // glColor3f(0, 1, 0);
  // for (int i = 0; i < io->numberofedges; i++) {
  //   const int ai = io->edgelist[i*2+0];
  //   int bi = io->edgelist[i*2+1];
  //   float2 a = vpoints[ai];
  //   float2 b;
  //   float2 v(io->normlist[i*2+0], io->normlist[i*2+1]);
  //   if (bi != -1) {
  //     b = vpoints[bi];
  //     glBegin(GL_LINES);
  //     glVertex2fv(a);
  //     glVertex2fv(b);
  //     glEnd();
  //   } else {
  //   }
  // }
  // //     // b = Intersect(a, v);
  // //     // a + tv = 0
  // //     float ta[] = { -1, -1, -1, -1 };
  // //     if (abs(v[0]) > 0) {
  // //       ta[0] = (buf-a[0])/v[0];
  // //       ta[1] = ((window_width-buf)-a[0])/v[0];
  // //     }
  // //     if (abs(v[1]) > 0) {
  // //       ta[2] = (buf-a[1])/v[1];
  // //       ta[3] = ((window_height-buf)-a[1])/v[1];
  // //     }
  // //     float t = 999999;
  // //     for (int i = 0; i < 4; ++i) {
  // //       if (ta[i] > 0)
  // //         t = min(t, ta[i]);
  // //     }
  // //     // b = a + v * t;
  // //     b = a + v * (t*2);
  // //     // b = a + v * window_width;
  // //     bi = vpoints.size();
  // //     vpoints.push_back(b);
  // //     edges.push_back(list<int>());
  // //   }
  // //   edges[ai].push_back(bi);
  // //   edges[bi].push_back(ai);

  // //   if (!show_delaunay) {
  // //     glLineWidth(1.0);
  // //     glBegin(GL_LINES);
  // //     glVertex2fv(a);
  // //     glVertex2fv(b);
  // //     glEnd();
  // //   }
  // // }

  // // if (show_delaunay) {
  // //   triangulateio* io = &mid;
  // //   const int n = io->numberofpoints;
  // //   vector<float2> tpoints(n);
  // //   for (int i = 0; i < n; i++) {
  // //     tpoints[i] = float2(io->pointlist[i*2+0], io->pointlist[i*2+1]);
  // //   }
  // //   for (int i = 0; i < io->numberofedges; i++) {
  // //     const int ai = io->edgelist[i*2+0];
  // //     int bi = io->edgelist[i*2+1];
  // //     float2 a = tpoints[ai];
  // //     float2 b = tpoints[bi];
  // //     glLineWidth(1.0);
  // //     glBegin(GL_LINES);
  // //     glVertex2fv(a);
  // //     glVertex2fv(b);
  // //     glEnd();
  // //   }
  // // }

  // // vector<float2> centroids = GetCentroids(vpoints, edges);
  // // if (!show_delaunay) {
  // //   glColor3f(1, 0, 0);
  // //   glBegin(GL_POINTS);
  // //   for (int i = 0; i < centroids.size(); ++i) {
  // //     glVertex2fv(centroids[i]);
  // //   }
  // //   glEnd();
  // // }

  // delete [] in.pointlist;
  // delete [] in.pointattributelist;
  // delete [] in.pointmarkerlist;
  // delete [] in.regionlist;

  // // return centroids;
}

bool OctViewer2::DrawVertexID(const int vi, const int2& p) const {
  BitmapString(vi, Oct2Obj(p),
               GL2D::kCenterJustify, GL2D::kCenterJustify);
  return true;
}

void OctViewer2::DrawVertexIDs() {
  oct::VisitVertices<2>(vertices, VertexIDCallback(this));
}

bool OctViewer2::DrawVertexLabel(const int vi, const int2& p) const {
  const double d = length(vertices.ClosestPoint(vi) - p);
  if (d < max_vertex_distance) {
    SetColor(vertices.Label(vi), red);
    if (vertices.Label(vi) == -1)
      glColor3fv(red.s);
    // glSquareWinWidth(Oct2Obj(p), pow(100,(d/(2*oct::kWidthd))));
    // glSquareWinWidth(Oct2Obj(p), 100*(d/(2*oct::kWidthd)));

    // const double w = Oct2Obj(0.3*vertices.Alpha(vi));
    const double w = Oct2Obj(0.07*d) + Win2Obj(4);
    glSquare(Oct2Obj(p)-w/2, w);
  }
  return true;
}

void OctViewer2::DrawVertexLabels() {
  glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
  oct::VisitVertices<2>(vertices, VertexLabelCallback(this));
}

bool OctViewer2::DrawVertexDistanceLine(const int vi, const int2& p) const {
  if (length(vertices.ClosestPoint(vi) - p) < max_vertex_distance) {
    const int2& cp = vertices.ClosestPoint(vi);
    glVertex2fv(Oct2Obj(p).s);
    glVertex2fv(Oct2Obj(cp).s);
  }
  return true;
}

void OctViewer2::DrawVertexDistanceLines() {
  glColor3f(0, 0.7, 0);
  glLineWidth(1.0);
  glBegin(GL_LINES);
  oct::VisitVertices<2>(vertices, DrawVertexDistanceLineCallback(this));
  glEnd();
}

struct DFWalkCallback : public M2Callback {
  DFWalkCallback(const OctViewer2* m2_,
                 const int2& base_,
                 const int width_,
                 vector<float3>* colors_,
                 vector<float>* distances_,
                 vector<float2>* points_,
                 vector<float2>* cpoints_, // closest points
                 vector<int>* labels_,
                 vector<float>* corner_distances_,
                 vector<float2>* corner_points_,
                 const float max_dist_)
      : M2Callback(m2_),
        base(base_), width(width_),
        colors(colors_), distances(distances_),
        points(points_), cpoints(cpoints_), labels(labels_),
        corner_distances(corner_distances_), corner_points(corner_points_),
        max_dist(max_dist_) {}
  bool operator()(const int vi, const int n_vi,
                  const int2& p, const int2& q,
                  // const oct::Direction<2>& d) {
                  const oct::Direction& d) {
    const float3 red = make_float3(1, 0, 0);
    const float norm_dist =
        std::min(1.0f, length(m2->Vertices().ClosestPoint(vi) - p)/(float)max_dist);
    float3 c =
        m2->RandomColor(m2->Vertices().Label(vi), red);
    c = c*(1.0f - norm_dist);
    colors->push_back(c);
    distances->push_back(norm_dist);
    points->push_back(m2->Oct2Obj(p));
    cpoints->push_back(m2->Oct2Obj(m2->Vertices().ClosestPoint(vi)));
    labels->push_back(m2->Vertices().Label(vi));

    if ((p.s[0]-base.s[0]) % width == 0 && (p.s[1]-base.s[1]) % width == 0) {
      corner_distances->push_back(norm_dist);
      corner_points->push_back(m2->Oct2Obj(p));
    }

    return true;
  }
  int2 base;
  int width;
  vector<float3>* colors;
  vector<float>* distances;
  vector<float2>* points;
  vector<float2>* cpoints;
  vector<int>* labels;

  vector<float>* corner_distances;
  vector<float2>* corner_points;
  float max_dist;
};

bool tex_init[] = { false, false, false, false, false, false };
GLuint tex_id[6];
unsigned char* texture = 0;

bool OctViewer2::DistanceFunctionVertex(const int vi, const int2& p) const {
  typedef oct::level_t level_t;
  // typedef oct::Direction<2> Direction_t;
  typedef oct::Direction Direction_t;

  if (!vertices.IsBase(vi)) return true;

  const level_t level = vertices.CellLevel(vi);
  const float w = oct::Level2CellWidth(level);
  // VisitVertices(vertices, DFDistCallback(this, &max_dist_oct, &max_dist_obj));
  // const float max_dist = oct::CellWidth(0) / 4;

  vector<float3> colors;
  vector<float> distances;
  vector<float2> points;
  vector<float2> cpoints;
  vector<float> corner_distances;
  vector<float2> corner_points;
  vector<int> labels;
  // float max_dist;
  DFWalkCallback dfwc(this, p, w, &colors, &distances,
                      &points, &cpoints, &labels,
                      &corner_distances, &corner_points,
                      max_dist_oct);
  dfwc = oct::WalkAroundFace<2>(vi, p, level, 0, 1, vertices, dfwc);

  const int2 p0 = convert_int2(Obj2Win(Oct2Obj(p+make_int2(0, 0))));
  const int2 p1 = convert_int2(Obj2Win(Oct2Obj(p+make_int2(w, w))));

  for (int y = maximum(p1.s[1], 0); y < min(p0.s[1], window_height); ++y) {
    for (int x = maximum(p0.s[0], 0); x < min(p1.s[0], window_width); ++x) {
      // Find the visible vertex with the closest point
      const float2 p0 = Win2Obj(make_float2(x, y));
      float min_dist = numeric_limits<float>::maximum();
      int min_label;
      for (int i = 0; i < cpoints.size(); ++i) {
        // const float d = (p0-cpoints[i]).norm2();
        const float d = length2(p0-cpoints[i]);
        if (d < min_dist) {
          min_dist = d;
          min_label = labels[i];
        }
      }

      float3 cf = make_float3(0);
      if (show_distance_field == 1) {
        cf = RandomColor(min_label, red);
        // cf = cf*(1.0f - sqrt(min_dist)/max_dist_obj);
      } else if (show_distance_field == 2) {
        const double tmp = std::min(1.0f, sqrt(min_dist)/max_dist_obj);
        cf = make_float3(tmp, tmp, tmp);
      } else {
        // cf = MeanValue(Win2Obj(float2(x, y)), points, colors);
        const double tmp =
            MeanValue(Win2Obj(make_float2(x, y)), points, distances);
        cf = make_float3(tmp, tmp, tmp);
      }
      const bool3 c = make_bool3(cf.s[0]*255, cf.s[1]*255, cf.s[2]*255);
      memcpy(&texture[(y*window_width+x)*3], c.s,
             sizeof(unsigned char)*3);
    }
  }
  return true;
}

void OctViewer2::DrawDistanceField() {
  const int width = window_width;
  const int height = window_height;
  const float2 p = Win2Obj(make_float2(0, 0));;
  const float2 q = Win2Obj(make_float2(width, height));
  // char* texture = gvd_texture;
  // if (show_distance_field == 2) {
  //   texture = df_texture;
  // }
  // if (!texture) {
  if (!tex_init[show_distance_field]) {
    texture = (unsigned char *)malloc(width * height * 3);
    memset(texture, 255, width*height*3*sizeof(unsigned char));

    oct::VisitVertices<2>(vertices, DistanceFunctionCallback(this));
  
    GLuint id;
    glGenTextures(1, &id);
    glBindTexture(GL_TEXTURE_2D, id);
    glTexEnvf( GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_DECAL );
    glTexParameterf( GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameterf( GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameterf( GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT );
    glTexParameterf( GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT );
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, width, height, 0, GL_RGB,
                 GL_UNSIGNED_BYTE, texture);
    free(texture);
    tex_init[show_distance_field] = true;
    tex_id[show_distance_field] = id;
  }

  glBindTexture(GL_TEXTURE_2D, tex_id[show_distance_field]);

  glColor3f(0, 1, 0);
  glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
  glEnable(GL_TEXTURE_2D);
  glBegin(GL_QUADS);
  glTexCoord2d(0, 0);
  glVertex2fv(p.s);
  glTexCoord2d(1, 0);
  glVertex2fv(make_float2(q.s[0], p.s[1]).s);
  glTexCoord2d(1, 1);
  glVertex2fv(q.s);
  glTexCoord2d(0, 1);
  glVertex2fv(make_float2(p.s[0], q.s[1]).s);
  glEnd();
  glDisable(GL_TEXTURE_2D);

  // glLineWidth(1.0);
  // glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
  // // glBegin(GL_QUADS);
  // VisitVertices(vertices, DistanceFunctionCallback(this));
  // // glEnd();
}

void OctViewer2::PrintStatistics() const {
  int num_cells = 0;
  oct::level_t max_level = 0;
  for (int i = 0; i < vertices.size(); ++i) {
    if (vertices.IsBase(i)) {
      ++num_cells;
      max_level = maximum(max_level, vertices.CellLevel(i));
    }
  }

  vector<string> lines;
  {
    stringstream ss;
    ss << "Max level: " << (int)max_level;
    lines.push_back(ss.str());
  } {
    stringstream ss;
    ss << "Quadtree cells: " << num_cells;
    lines.push_back(ss.str());
  } {
    // const int2 mouse_oct = convert_int2(Obj2Oct(mouse_obj));
    // stringstream ss;
    // // ss << mouse_oct;
    // ss << mouse_obj;
    // // lines.push_back(ss.str());
  } {
    if (show_poly_vertices == 0) {
      lines.push_back("No vertices");
    } else if (show_poly_vertices == 1) {
      lines.push_back("Latest Karras points");
    } else if (show_poly_vertices == 2) {
      lines.push_back("Next Karras points");
    } else if (show_poly_vertices == 3) {
      lines.push_back("Polygon / cell intersections");
    }
  }

  int offset = 2;
  for (const string& s : lines) {
    BitmapString(s, Win2Obj(make_float2(2, window_height-offset)),
                 kLeftJustify, kBottomJustify);
    offset += 17;
  }
}

void OctViewer2::HelpString(const string msg, const int i) const {
  static const int sep = 17;
  static const int y = 2;
  static const int x = 2;
  BitmapString(msg, Win2Obj(make_float2(x, y+sep*i)),
               kLeftJustify, kTopJustify);
}

void OctViewer2::PrintHelp() const {
  int i = 0;
  HelpString("h - toggle help", i++);
  if (show_help) {
    HelpString("View", i++);
    HelpString("  m - toggle GVD", i++);
    HelpString("  o - toggle quadtree", i++);
    HelpString("  v - toggle vertex labels", i++);
    HelpString("  l - toggle closest point line", i++);
    HelpString("  i - toggle vertex IDs", i++);
    HelpString("  p - toggle statistics", i++);
    HelpString("  Ctrl+left drag - zoom", i++);
    HelpString("  r - reset view", i++);
    HelpString("  H - toggle advanced help", i++);
    HelpString("Polygon entry", i++);
    HelpString("  z - clear last polygon drawn", i++);
    HelpString("  c - clear all polygons", i++);
    HelpString("  e - entry mode", i++); // ????
    HelpString("  W - write polygons to poly-*.dat", i++);
    HelpString("GVD computation", i++);
    HelpString("  f/d - increment/decrement max octree level", i++);
    HelpString("  V - toggle full subdivide", i++);
    HelpString("  B - toggle buffer", i++);
    // HelpString("C - ???", i++);
    HelpString("q - quit", i++);
    // HelpString("", i++);
    // HelpString("", i++);
    // HelpString("", i++);
  } else if (show_advanced_help) {
    HelpString("view", i++);
    HelpString("  P - toggle min path", i++); // show_path
    HelpString("  Alt+O/o - increase/decrease octree color", i++); // show_path
    HelpString("  H - toggle advanced help", i++);
    HelpString("misc", i++);
    HelpString("g - polygon mode", i++); // polygon_mode
    HelpString("s/S - increment/decrement ambiguous max level", i++);
    HelpString("j/k - increment/decrement max vertex distance", i++);
    // HelpString("u - show distance field", i++); // ????
    HelpString("t - test function", i++);
  }
}

void OctViewer2::Display() {
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  gluOrtho2D(win_obj.min().s[0], win_obj.max().s[0],
             win_obj.min().s[1], win_obj.max().s[1]);

  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
  // gluLookAt(0, 0, 100,
  //           0, 0, 0,
  //           0, 1, 0);
  // glTranslatef(0, 0, 100-2.7475);

  // Presentation mode
  // glTranslatef(view_obj.min()[0]/2, view_obj.min()[1]/2, 0);
  // glScalef(.5*view_obj.size()[0],
  //          .5*view_obj.size()[1],
  //          1);
  // glTranslatef(-view_obj.min()[0]/2, -view_obj.min()[1]/2, 0);

  glEnable(GL_LINE_SMOOTH);
  glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);
  glEnable(GL_BLEND);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

  // build octree
  if (dirty) {
    Timer t("Building octree");
    BuildOctree();
    t.stop();

    // for (int i = 0; i < 6; ++i) {
    //   tex_init[i] = false;
    // }
    dirty = false;
  }

  // Draw octree
  if (show_octree) {
    DrawOctree();
  }

  if (show_poly_vertices == 1) {
    // Show in blue the points that were used in the latest Karras run
    for (const floatn& p : karras_points) {
      glColor3f(0.0, 0.0, 0.7);
      glSquareCentered(p, Win2Obj(4));
    }
  } else if (show_poly_vertices == 2) {
    // Show fit boxes
    glColor3f(0.7, 0.7, 0.7);
    glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
    for (int i = 0; i < _origins.size(); ++i) {
      const floatn& o = _origins[i];
      const float& d = _lengths[i];
      glSquare(Oct2Obj(convert_intn(o)), Oct2Obj(d));
    }
    // Show in red any points that will be added in the next Karras run
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
    glColor3f(1.0, 0.0, 0.0);
    for (const intn& p : extra_qpoints) {
      glSquareCentered(Oct2Obj(p), Win2Obj(4));
    }
  } else if (show_poly_vertices == 3) {
    // Show in green anywhere the polygons intersect existing cells
    for (const floatn& pf : intersections) {
      const intn p = convert_intn(pf);
      glColor3f(1.0, 0.0, 0.0);
      glSquareCentered(Oct2Obj(p), Win2Obj(4));
    }
  }

  // if (show_distance_field > 0) {
  //   DrawDistanceField();
  // }

  // // draw octree
  // // glColor3f(0.0, 0.0, 0.0);
  // glColor3f(0.7, 0.7, 0.7);
  // glLineWidth(1.0);
  // if (show_octree) {
  //   DrawOctree();
  // }

  if (polygon_mode > 0) {
    glColor3f(0.9, 0.9, 0.9);
    // draw the polygons
    if (polygon_mode == 1) {
      glLineWidth(2.0);
    } else {
      glLineWidth(2);
    }
    glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
    for (int j = 0; j < polygons.size(); ++j) {
      if (polygon_mode == 1) {
        SetColor(j, red);
      }
      glBegin(GL_LINE_STRIP);
      for (int i = 0; i < polygons[j].size(); ++i) {
        glVertex2fv(polygons[j][i].s);
      }
      glEnd();
    }
  }

  SetColor(polygons.size(), red);
  glBegin(GL_LINE_STRIP);
  for (int i = 0; i < verts.size(); ++i) {
    glVertex2fv(verts[i].s);
  }
  glEnd();

  // // Draw other stuff
  // if (show_gvd) {
  //   DrawGVD();
  // }
  // if (show_path) {
  //   DrawPath();
  // }
  // if (show_voronoi) {
  //   DrawVoronoi();
  // }
  // if (show_vertex_ids) {
  //   DrawVertexIDs();
  // }
  // if (show_closest_point_line) {
  //   DrawVertexDistanceLines();
  // }
  // if (show_vertex_labels) {
  //   DrawVertexLabels();
  // }
  if (show_statistics) {
    PrintStatistics();
  }
  // PrintHelp();

  // // Debug circles
  // glLineWidth(1.0);
  // for (int vi = 0; vi < vertex_locations.size(); ++vi) {
  //   if (vertex_circle[vi]) {
  //     const float2& v = vertex_locations[vi];
  //     const float2 cp = Oct2Obj(vertices.ClosestPoint(vi));
  //     const float d = length(v-cp);
  //     glColor3f(0, 0, 1);
  //     glSquareWinWidth(v, 5);
  //     glBegin(GL_LINE_LOOP);
  //     for (int i = 0; i < 360; ++i) {
  //       const float t = i*M_PI/180.0;
  //       glVertex2fv((v + make_float2(d*cos(t), d*sin(t))).s);
  //     }
  //     glEnd();
  //     glColor3f(0, 1, 0);
  //     glBegin(GL_LINES);
  //     glVertex2fv(v.s);
  //     glVertex2fv(cp.s);
  //     glEnd();
  //   }
  // }

  // glLineWidth(2.0);
  // vector<float2> up(middle_up);
  // if (middle_down.size() > middle_up.size()) {
  //   up.push_back(mouse_obj);
  // }
  // // if (middle_down != float2()) {
  // for (int i = 0; i < middle_down.size(); ++i) {
  //   glColor3f(0, .3, .3);
  //   glBegin(GL_LINES);
  //   glVertex2fv(middle_down[i].s);
  //   glVertex2fv(up[i].s);
  //   glEnd();

  //   const float d = length(up[i]-middle_down[i]);
  //   glBegin(GL_LINE_LOOP);
  //   for (int j = 0; j < 360; ++j) {
  //     const float t = j*M_PI/180.0;
  //     glVertex2fv((middle_down[i] + make_float2(d*cos(t), d*sin(t))).s);
  //   }
  //   glEnd();
  // }

  // Draw rubberband
  if (rubberband_start != float2() && rubberband_cur != float2()) {
    glColor3f(0.5, 0.5, 0.5);
    glLineWidth(1.0);
    glBegin(GL_LINE_STRIP);
    glVertex2fv(rubberband_start.s);
    glVertex2f(rubberband_start.s[0], rubberband_cur.s[1]);
    glVertex2fv(rubberband_cur.s);
    glVertex2f(rubberband_cur.s[0], rubberband_start.s[1]);
    glVertex2fv(rubberband_start.s);
    glEnd();
  }
}
