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

#include <cstdlib>
#include <string>
#include <vector>
#include <sstream>
#include "./io.h"

using namespace std;

namespace {

string Dir(const string& fn) {
  int i = fn.rfind('/');
  if (i == string::npos) {
    i = fn.find('\\');
  }
  if (i != string::npos) {
    return fn.substr(0, i);
  }
  return ".";
}

float ToFloat(const string& s) {
  return ::atof(s.c_str());
}

int ToInt(const string& s) {
  return ::atoi(s.c_str());
}

string Trim(const string& s) {
  int i = 0;
  int j = s.size()-1;
  while (i < s.size() && (s[i] == ' ' || s[i] == '\t')) {
    i++;
  }
  while (j > i && (s[j] == ' ' || s[j] == '\t')) {
    j--;
  }
  if (i == j) return "";
  return s.substr(i, j+1-i);
}

string GetLine(istream* pin) {
  istream& in = *pin;
  string line;
  do {
    getline(in, line);
    line = Trim(line);
  } while (!in.eof() && line[0] == '#');
  return line;
}

void Parse(const string& line, string* t, vector<string>* values) {
  stringstream ss(line);
  ss >> *t;
  while (!ss.eof()) {
    string s;
    ss >> s;
    values->push_back(s);
  }
}

void ParseVertex(const string& s, int* vi, int* ti, int* ni) {
  // Find the offsets of the texture and normal indices
  int t_offset = s.find('/');
  int n_offset = s.size();
  if (t_offset != string::npos) {
    n_offset = s.find('/', t_offset+1);
    if (n_offset == string::npos) {
      n_offset = s.size();
    }
  } else {
    t_offset = s.size();
  }

  *vi = ToInt(s.substr(0, t_offset)) - 1;
  string ts, ns;
  if (t_offset < n_offset-1) {
    ts = s.substr(t_offset+1, n_offset-(t_offset+1));
    if (n_offset < s.size()-1) {
      ns = s.substr(n_offset+1);
    }
  }
  if (!ts.empty()) {
    *ti = ToInt(ts) - 1;
  }
  if (!ns.empty()) {
    *ni = ToInt(ns) - 1;
  }
}

void ParseMaterial(const string& fn, Mesh& mesh) {
  const string dir = Dir(fn);
  ifstream in(fn.c_str());
  string line = GetLine(&in);
  int cnt = -1, line_cnt = 0;
  while (!in.eof()) {
    string t;
    vector<string> values;
    Parse(line, &t, &values);
    if (line.empty()) {
      // do nothing
    } else if (t == "newmtl") {
      mesh.new_material(++cnt, values[0]);
    } else if (t == "Ka") {
      mesh.set_ambient(cnt,
          make_float3(ToFloat(values[0]), ToFloat(values[1]), ToFloat(values[2])));
    } else if (t == "Kd") {
      mesh.set_diffuse(cnt,
          make_float3(ToFloat(values[0]), ToFloat(values[1]), ToFloat(values[2])));
    } else if (t == "Ks") {
      mesh.set_specular(cnt,
          make_float3(ToFloat(values[0]), ToFloat(values[1]), ToFloat(values[2])));
    } else if (t == "Ns") {
      mesh.set_specular_coeff(cnt, ToFloat(values[0]));
    } else if (t == "map_Ka" || t == "map_Kd") {
      string tfn = values[0];
      if (!dir.empty()) {
        tfn = dir + "/" + tfn;
      }
      mesh.set_texture(cnt, tfn);
    } else {
      // cout << "Unrecognized material line: " << line
      //      << " " << line_cnt << endl;
    }
    line = GetLine(&in);
    line_cnt++;
  }
}

}  // end empty namespace

bool ParseObj(const string& fn, Mesh& mesh) {
  const string dir = Dir(fn);
  ifstream in(fn.c_str());
  if (!in.good()) {
    return false;
  }
  string line = GetLine(&in);
  int cnt = 0;
  while (!in.eof()) {
    if (line.empty()) {
      // do nothing
    } else if (line.find("mtllib") != string::npos) {
      string mfn = fn;
      if (!dir.empty()) {
        mfn = dir + "/" + line.substr(7);
      }
      ParseMaterial(mfn, mesh);
    } else if (line.find("usemtl") != string::npos) {
      const string fn = line.substr(7);
      mesh.set_cur_material(fn);
    } else if (line[0] == 'g') {
      // cout << "group " << line.substr(2) << endl;
    } else if (line[0] == 'o') {
      // cout << "object " << line.substr(2) << endl;
    } else if (line[0] == 's') {
      // cout << "smooth shading " << line.substr(2) << endl;
    } else {
      string t;
      vector<string> values;
      Parse(line, &t, &values);
      if (t == "v") {
        mesh.AddVertex(make_float3(
            ToFloat(values[0]), ToFloat(values[1]), ToFloat(values[2])));
      } else if (t == "vt") {
        float v[3];
        for (int i = 0; i < 3; ++i) {
          v[i] = ToFloat(values[i]);
        }
        // mesh.AddTextureVertex(float3(v[0], v[1], v[2]));
      } else if (t == "f") {
        vector<int> vertices(values.size());
        vector<int> texture(values.size(), -1);
        vector<int> normals(values.size(), -1);
        if (values.size() != 3)
          throw logic_error(
              "This software supports only triangulated surfaces");
        for (int i = 0; i < values.size(); ++i) {
          ParseVertex(values[i], &vertices[i], &texture[i], &normals[i]);
        }
        // mesh.AddPolygon(vertices, texture);
        mesh.AddTriangle(make_triangle(vertices[0], vertices[1], vertices[2]));
      } else if (t == "vn") {
        // Ignore normals
      } else {
        cout << "Unknown type " << t << " on line " << cnt+1 << endl;
      }
    }
    line = GetLine(&in);
    cnt++;
  }
  // mesh.shrink_to_fit();
  return true;
}

void WriteObj(std::ostream& out, const Mesh& mesh) {
  const vector<float3>& vertices = mesh.vertices();
  const vector<float3>& colors = mesh.colors();
  const vector<Triangle>& triangles = mesh.triangles();

  for (int i = 0; i < vertices.size(); ++i) {
    const float3& v = vertices[i];
    const float3& c = colors[i];
    out << "v " 
        << std::setprecision(12) << v.s[0] << " " 
        << std::setprecision(12) << v.s[1] << " " 
        << std::setprecision(12) << v.s[2] << " "
        << std::setprecision(12) << c.s[0] << " "
        << std::setprecision(12) << c.s[1] << " "
        << std::setprecision(12) << c.s[2] << " "
        << endl;
  }
  for (int i = 0; i < triangles.size(); ++i) {
    const Triangle& t = triangles[i];
    if (t.s[0] != t.s[1] && t.s[0] != t.s[2] && t.s[1] != t.s[2]) {
      out << "f ";
      for (int j = 0; j < 3; ++j) {
        out << (t.s[j]+1) << " ";
      }
      out << endl;
    }
  }
}

void WritePly(std::ostream& out, const Mesh& mesh) {
  const vector<float3>& vertices = mesh.vertices();
  const vector<float3>& colors = mesh.colors();
  const vector<Triangle>& triangles = mesh.triangles();

  out << "ply" << endl;
  out << "format ascii 1.0" << endl;
  // out << "element vertex 8" << endl;
  out << "element vertex " << vertices.size() << endl;
  out << "property float x" << endl;
  out << "property float y" << endl;
  out << "property float z" << endl;
  // out << "property uchar red" << endl;
  // out << "property uchar green" << endl;
  // out << "property uchar blue" << endl;
  out << "property float red" << endl;
  out << "property float green" << endl;
  out << "property float blue" << endl;
  // out << "element face 7" << endl; // number of triangles
  out << "element face " << triangles.size() << endl;
  out << "property list uchar int vertex_indices" << endl; // # verts per face
  // out << "element edge 5" << endl; // 5 edges
  // out << "property int vertex1" << endl; // index to first vertex of edge
  // out << "property int vertex2" << endl; // index to second vertex
  // out << "property uchar red" << endl; // start of edge color
  // out << "property uchar green" << endl;
  // out << "property uchar blue" << endl;
  out << "end_header" << endl;

  for (int i = 0; i < vertices.size(); ++i) {
    const float3& v = vertices[i];
    const float3& c = colors[i];
    out << std::setprecision(12) << v.s[0] << " " 
        << std::setprecision(12) << v.s[1] << " " 
        << std::setprecision(12) << v.s[2] << " "
        << std::setprecision(12) << c.s[0] << " "
        << std::setprecision(12) << c.s[1] << " "
        << std::setprecision(12) << c.s[2] << " "
        << endl;
  }
  for (int i = 0; i < triangles.size(); ++i) {
    const Triangle& t = triangles[i];
    if (t.s[0] != t.s[1] && t.s[0] != t.s[2] && t.s[1] != t.s[2]) {
      out << "3 ";
      for (int j = 0; j < 3; ++j) {
        out << t.s[j] << " ";
      }
      out << endl;
    }
  }
}

void WriteDensity(std::ostream& out, const Mesh& mesh) {
  // const vector<float3>& vertices = mesh.vertices();
  // const vector<std::vector<int> >& polygons = mesh.polygons();

  // vector<Triangle> triangles;
  // for (int i = 0; i < polygons.size(); ++i) {
  //   const vector<int>& polygon = polygons[i];
  //   if (polygon.size() > 3)
  //     throw logic_error("WriteDensity does not support arbitrary polygons");
  //   triangles.push_back(Triangle(polygon[0], polygon[1], polygon[2]));
  // }

  // const vector<vector<Triangle> > v2t =
  //     build_v2t(triangles.begin(), triangles.end());

  // out << vertices.size() << endl;
  // for (int i = 0; i < vertices.size(); ++i) {
  //   const vector<Triangle>& tris = v2t[i];
  //   double area = 0;
  //   for (int j = 0; j < tris.size(); ++j) {
  //     const Triangle& t = tris[j];
  //     const float3& p = vertices[t[0]];
  //     const float3& q = vertices[t[1]];
  //     const float3& r = vertices[t[2]];
  //     area += ((r-p)^(q-p)).norm();
  //   }
  //   if (area > 0)
  //     out << 1/(area*area) << endl;
  //   else
  //     out << 0 << endl;
  // }
}
