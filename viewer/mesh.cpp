#include <cstring>
#include <iostream>
#include <vector>
#include <set>
#include <climits>

#include "./mesh.h"
#include "../orientation.h"

using namespace std;

namespace {
#define BUFFER_OFFSET(i) ((char*)0 + (i))

float3 Normal(const float3& a, const float3& b, const float3& c) {
  // return (b-a)^(c-a);
  return cross(b-a, c-a);
}

}

// Taken from
// http://nehe.gamedev.net/tutorial/vertex_buffer_objects/22002/
bool IsExtensionSupported(const char* szTargetExtension)
{
  const unsigned char *pszExtensions = NULL;
  const unsigned char *pszStart;
  unsigned char *pszWhere, *pszTerminator;

  // Extension names should not have spaces
  pszWhere = (unsigned char *) strchr( szTargetExtension, ' ' );
  if( pszWhere || *szTargetExtension == '\0' )
    return false;

  // Get Extensions String
  pszExtensions = glGetString( GL_EXTENSIONS );
  if (!pszExtensions) {
    return false;
    cout << "no extensions" << endl;
  }

  // Search The Extensions String For An Exact Copy
  pszStart = pszExtensions;
  for(;;)
  {
    pszWhere = (unsigned char *) strstr(
        (const char *) pszStart, szTargetExtension);
    if( !pszWhere )
      break;
    pszTerminator = pszWhere + strlen( szTargetExtension );
    if( pszWhere == pszStart || *( pszWhere - 1 ) == ' ' )
      if( *pszTerminator == ' ' || *pszTerminator == '\0' )
        return true;
    pszStart = pszTerminator;
  }
  return false;
}

Mesh::Mesh() : _color(false), _num_triangles(0) {
  _vbo_initialized = false;
  _cur_mtl = -1;

  _buffer_names[0] = 0;
  _buffer_names[1] = 0;
}

Mesh::Mesh(const Mesh& rhs) {
  Copy(rhs);
}

Mesh::~Mesh() {
  if (_vbo_initialized) {
    glDeleteBuffers(2, _buffer_names);
  }
}

void Mesh::operator=(const Mesh& rhs) {
  Copy(rhs);
}

void Mesh::Copy(const Mesh& rhs) {
  _vertices = rhs._vertices;
  _normals = rhs._normals;
  _colors = rhs._colors;
  _triangles = rhs._triangles;
  _materials = rhs._materials;
  _color = rhs._color;
  _cur_mtl = rhs._cur_mtl;
  _bb = rhs._bb;
  _num_triangles = rhs._num_triangles;

  _vbo_initialized = false;
  _buffer_names[0] = 0;
  _buffer_names[1] = 0;
  
  _centroid = rhs._centroid;
}

void Mesh::InitVBO() const {
  _vbo_initialized = true;
  bool use_vbo = IsExtensionSupported("GL_ARB_vertex_buffer_object");
  if (!use_vbo) {
    cerr << "**********************" << endl;
    cerr << " VBO not supported " << endl;
    cerr << "**********************" << endl;
    return;
  }

  glDeleteBuffers(2, _buffer_names);
  glGenBuffers(2, _buffer_names);

  // Send vertex buffer object
  if (_color) {
    MeshVertexC* mverts = new MeshVertexC[_vertices.size()];
    for (int i = 0; i < _vertices.size(); ++i) {
      mverts[i] = MeshVertexC(_vertices[i], normal(i), _colors[i]);
    }
    glBindBuffer(GL_ARRAY_BUFFER, _buffer_names[0]);
    glBufferData(
        GL_ARRAY_BUFFER, _vertices.size()*sizeof(MeshVertexC),
        mverts, GL_STREAM_DRAW);
    delete [] mverts;
  } else {
    MeshVertex* mverts = new MeshVertex[_vertices.size()];
    for (int i = 0; i < _vertices.size(); ++i) {
      mverts[i] = MeshVertex(_vertices[i], normal(i));
    }
    glBindBuffer(GL_ARRAY_BUFFER, _buffer_names[0]);
    glBufferData(
        GL_ARRAY_BUFFER, _vertices.size()*sizeof(MeshVertex),
        mverts, GL_STREAM_DRAW);
    delete [] mverts;
  }

  // Send index buffer object
  GLuint* indices = new GLuint[_triangles.size()*3];
  for (int i = 0; i < _triangles.size(); ++i) {
    for (int j = 0; j < 3; ++j) {
      indices[i*3+j] = _triangles[i].s[j];
    }
  }
  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, _buffer_names[1]);
  glBufferData(
      GL_ELEMENT_ARRAY_BUFFER, _triangles.size()*3*sizeof(GLuint),
      indices, GL_STREAM_DRAW);
  delete [] indices;
}

void Mesh::Display(bool wireframe) const {
  if (!_vbo_initialized) {
    InitVBO();
  }

  // glBindBuffer(GL_ARRAY_BUFFER, _buffer_names[0]);
  // glVertexPointer(3, GL_FLOAT, sizeof(MeshVertex), 0);
  // glNormalPointer(GL_FLOAT, sizeof(MeshVertex),
  //                 BUFFER_OFFSET(3*sizeof(GLfloat)));
  // glColorPointer(4, GL_UNSIGNED_BYTE, sizeof(MeshVertex),
  //                BUFFER_OFFSET(6*sizeof(GLfloat)));

  glBindBuffer(GL_ARRAY_BUFFER, _buffer_names[0]);
  // glVertexPointer(3, GL_FLOAT, sizeof(MeshVertex), 0);
  // glNormalPointer(GL_FLOAT, sizeof(MeshVertex),
  //                 BUFFER_OFFSET(3*sizeof(GLfloat)));
  if (_color) {
    glVertexPointer(3, GL_FLOAT, sizeof(MeshVertexC), 0);
    glNormalPointer(GL_FLOAT, sizeof(MeshVertexC),
                    BUFFER_OFFSET(3*sizeof(GLfloat)));
    glColorPointer(4, GL_UNSIGNED_BYTE, sizeof(MeshVertexC),
                   BUFFER_OFFSET(6*sizeof(GLfloat)));
  } else {
    glVertexPointer(3, GL_FLOAT, sizeof(MeshVertex), 0);
    glNormalPointer(GL_FLOAT, sizeof(MeshVertex),
                    BUFFER_OFFSET(3*sizeof(GLfloat)));
  }

  // glBindBuffer(GL_ARRAY_BUFFER, _vertices->buffer_name());
  // glVertexPointer(3, GL_FLOAT, sizeof(MeshVertex), 0);

  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, _buffer_names[1]);

  // Rendering
  glEnableClientState(GL_VERTEX_ARRAY);
  glEnableClientState(GL_NORMAL_ARRAY);
  if (_color)
    glEnableClientState(GL_COLOR_ARRAY);
  else
    glDisableClientState(GL_COLOR_ARRAY);

  if (wireframe) {
    glLineWidth(1.0f);
    glDisable(GL_LIGHTING);
    glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
    glColor3f(.2f, .2f, .2f);
    glDisableClientState(GL_COLOR_ARRAY);
    glDrawElements(GL_TRIANGLES, _num_triangles*3, GL_UNSIGNED_INT, 0);
    // glEnableClientState(GL_COLOR_ARRAY);
  } else {
    // Material properties
    if (!_materials.empty() && !_color) {
      const Material& m = _materials[0];
      glMaterialfv(GL_FRONT_AND_BACK, GL_EMISSION,
                   make_float4(0, 0, 0, 1).s);
      glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR,
                   make_float4(m.specular(), 1).s);
      glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE,
                   make_float4(m.diffuse(), 1).s);
      glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT,
                   make_float4(m.ambient(), 1).s);
      glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS,
                  m.specular_coeff());
    }

    // Surface
    if (_color)
      glEnable(GL_COLOR_MATERIAL);
    glEnable(GL_LIGHTING);
    glEnable(GL_LIGHT0);
    glPolygonMode(GL_FRONT, GL_FILL);
    glShadeModel(GL_SMOOTH);
    glEnable(GL_POLYGON_OFFSET_FILL);
    glPolygonOffset(1.0, 1.0);
    glDrawElements(GL_TRIANGLES, _num_triangles*3, GL_UNSIGNED_INT, 0);
    glDisable(GL_POLYGON_OFFSET_FILL);
    if (_color)
      glDisable(GL_COLOR_MATERIAL);
  }

  if (_color)
    glDisableClientState(GL_COLOR_ARRAY);
  glDisableClientState(GL_NORMAL_ARRAY);
  glDisableClientState(GL_VERTEX_ARRAY);
}

void Mesh::DisplayNormals() const {
  // glDisable(GL_LIGHTING);
  // glLineWidth(2);
  // glColor3f(1, 0, 0);
  // glBegin(GL_LINES);
  // for (int i = 0; i < _vertices->size(); ++i) {
  //   const float3& vertex = _vertices->at(i);
  //   float3 n = normal(i) / 2;
  //   glVertex3fv(vertex);
  //   glVertex3fv(vertex+n);
  // }
  // glEnd();
  // glEnable(GL_LIGHTING);
}

void Mesh::AddVertex(const float3& v, float3 color) {
  _vertices.push_back(v);
  _colors.push_back(color);
  _bb(v);
}

void Mesh::SetVertexColor(const int i, const float3& color) {
  if (i >= _colors.size()) {
    cerr << "i = " << i << " _colors.size() = " << _colors.size() << endl;
    throw logic_error("Invalid index in SetVertexColor");
  }
  _colors[i] = color;
  _vbo_initialized = false;
}

void Mesh::Translate(const float3& t) {
  _bb = BoundingBox3f();
  for (int i = 0; i < _vertices.size(); ++i) {
    _vertices[i] += t;
    _bb(_vertices[i]);
  }
  _vbo_initialized = false;
}

void Mesh::SetAndCollapse(const vector<float3>& vertices,
                          const vector<Triangle>& triangles) {
  _vertices.clear();
  vector<int> v2v(vertices.size(), -1);
  for (int i = 0; i < triangles.size(); ++i) {
    const Triangle& t = triangles[i];
    for (int j = 0; j < 3; ++j) {
      const int old_idx = t.s[j];
      int new_idx = v2v[old_idx];
      if (new_idx == -1) {
        new_idx = _vertices.size();
        v2v[old_idx] = new_idx;
        _vertices.push_back(vertices[old_idx]);
      }
    }
    _triangles.push_back(make_triangle(v2v[t.s[0]], v2v[t.s[1]], v2v[t.s[2]]));
  }
  _num_triangles = _triangles.size();
  _vbo_initialized = false;
}

void Mesh::AddTriangle(const Triangle& t) {
  _triangles.push_back(t);
  _num_triangles++;
  _vbo_initialized = false;
}

struct MyPair {
  MyPair(const float3& v_, const int i_) : v(v_), i(i_) {}
  bool operator<(const MyPair& rhs) const {
    return v < rhs.v;
  }
  bool operator==(const MyPair& rhs) const {
    return v == rhs.v;
  }
  float3 v;
  int i;
};

void Mesh::compute_normals() {
  vector<vector<int> > v2p(_vertices.size());
  for (int i = 0; i < _triangles.size(); ++i) {
    const Triangle& t = _triangles[i];
    for (int j = 0; j < 3; ++j) {
      v2p[t.s[j]].push_back(i);
    }
  }

  _normals.resize(_vertices.size());
  for (int i = 0; i < _vertices.size(); ++i) {
    // Loop through each adjacent polygon and find
    // its weighted normal
    // float3 n = float3::zero();
    float3 n = make_float3(0);
    for (int j = 0; j < v2p[i].size(); ++j) {
      const Triangle& p = _triangles[v2p[i][j]];
      n += Normal(
          _vertices.at(p.s[0]), _vertices.at(p.s[1]), _vertices.at(p.s[2]));
    }
    // _normals[i] = n.unit();
    _normals[i] = normalize(n);
  }
}

void Mesh::OrientPolygons() {
  orient_lean(_triangles);
  if (_num_triangles != _triangles.size()) {
    cerr << "num_triangles changed in orient: " << _num_triangles
         << " " << _triangles.size() << endl;
  }
  _num_triangles = _triangles.size();
}

void Mesh::InvertOrientation() {
  // check_lean();
  for (int i = 0; i < _triangles.size(); ++i) {
    _triangles[i] = triangle_inverted(_triangles[i]);
  }
  for (int i = 0; i < _normals.size(); ++i) {
    _normals[i] = _normals[i] * -1;
  }
  _vbo_initialized = false;
}

float Mesh::pick(const float3& p, const float3& v) {
  // check_lean();
  float mint = 999999999;
  // for (int i = 0; i < _polygons.size(); ++i) {
  for (int i = 0; i < _triangles.size(); ++i) {
    // Get the plane
    const float3 P = _vertices.at(_triangles[i].s[0]);
    const float3 Q = _vertices.at(_triangles[i].s[1]);
    const float3 R = _vertices.at(_triangles[i].s[2]);
    // const float3 n = ((Q-P)^(R-P)).unit();
    // const float3 n = normalize((Q-P)^(R-P));
    const float3 n = normalize(cross(Q-P, R-P));
    // const float D = -n*P;
    const float D = -dot(n, P);
    // assert(-n*Q == D);
    // if (fabs(-n*Q - D) > 0.0001) throw logic_error("nqd");
    if (fabs(-dot(n, Q) - D) > 0.0001) throw logic_error("nqd");
    // (p+tv)*n+D=0
    // t = -(p*n+D)/v*n
    // if (v*n != 0) {
    if (dot(v, n) != 0) {
      // const float t = -(p*n+D)/(v*n);
      const float t = -(dot(p, n)+D)/dot(v, n);
      // const float3 C = p+t*v;
      const float3 C = p+v*t;
      // if (((Q-P).unit()^(C-P).unit())*n > 0 &&
      //     ((R-Q).unit()^(C-Q).unit())*n > 0 &&
      //     ((P-R).unit()^(C-R).unit())*n > 0) {
      const float3 QP = normalize(Q-P);
      const float3 CP = normalize(C-P);
      const float3 RQ = normalize(R-Q);
      const float3 CQ = normalize(C-Q);
      const float3 PR = normalize(P-R);
      const float3 CR = normalize(C-R);
      // if (dot((normalize(Q-P)^normalize(C-P)), n) > 0 &&
      //     dot((normalize(R-Q)^normalize(C-Q)), n) > 0 &&
      //     dot((normalize(P-R)^normalize(C-R)), n) > 0) {
      if (dot(cross(QP,CP), n) > 0 &&
          dot(cross(RQ,CQ), n) > 0 &&
          dot(cross(PR,CR), n) > 0) {
        if (t < mint) {
          mint = t;
        }
      }
    }
  }
  return mint;
}

void Mesh::print_stats() const {
  const int vu = (_vertices.capacity()*sizeof(float3));
  const int nu = (_normals.capacity()*sizeof(float3));
  const int tu = (_triangles.capacity()*sizeof(Triangle));
  const int mu = (_materials.capacity()*sizeof(Material));

  cout << "vertices size = " << _vertices.size() << endl;
  cout << "vertices capacity = " << _vertices.capacity() << endl;
  cout << "vertices usage = " << vu << endl;

  cout << "normals size = " << _normals.size() << endl;
  cout << "normals capacity = " << _normals.capacity() << endl;
  cout << "normals usage = " << nu << endl;

  cout << "triangles size = " << _triangles.size() << endl;
  cout << "triangles capacity = " << _triangles.capacity() << endl;
  cout << "triangles usage = " << tu << endl;

  cout << "materials usage = " << mu << endl;

  cout << "total usage = " << (vu+nu+tu+mu) << endl;
}

float3 Mesh::Centroid() const {
  if (_centroid == float3()) {
    float3 c_sum = make_float3(0);
    float area_sum = 0;
    for (int i = 0; i < _triangles.size(); ++i) {
      const float3 P = _vertices.at(_triangles[i].s[0]);
      const float3 Q = _vertices.at(_triangles[i].s[1]);
      const float3 R = _vertices.at(_triangles[i].s[2]);
      const float3 c = (P + Q + R)/3;
      // const float area = length((Q-P)^(R-P))/2;
      const float area = length(cross(Q-P, R-P))/2;
      c_sum += (c * area);
      area_sum += area;
    }
    _centroid = c_sum / area_sum;
  }
  return _centroid;
}
