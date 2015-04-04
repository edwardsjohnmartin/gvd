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

#ifndef __MESH_H__
#define __MESH_H__

#include <vector>
#include <string>

#include "../opencl/vec.h"
#include "../opencl/triangle.h"

#include "../bb.h"
#include "./material.h"
#include "./common.h"
#include "../shared_ptr.h"

struct GLfloat3 {
  GLfloat s[3];
};

static GLfloat3 make_glfloat3(const float3& f) {
  GLfloat3 glf = { f.x, f.y, f.z };
  return glf;
}

struct MeshVertexC {
  MeshVertexC() {}
  MeshVertexC(const float3& p_, const float3& n_, const bool3& c_)
      : p(make_glfloat3(p_)), n(make_glfloat3(n_)) {
    c = make_bool4(c_, 255);
  }
  MeshVertexC(const float3& p_, const float3& n_, const float3& c_)
      : p(make_glfloat3(p_)), n(make_glfloat3(n_)) {
    for (int i = 0; i < 3; ++i)
      c.s[i] = static_cast<GLbyte>(c_.s[i]*255);
    c.s[3] = 255;
  }
  // float3 p;  // point - offset 0
  // float3 n;  // normal - offset 12
  GLfloat3 p;  // point - offset 0
  GLfloat3 n;  // normal - offset 12
  bool4 c;  // color - offset 24
  GLbyte padding[4];  // to make it 32 bytes
};

struct MeshVertex
{
  MeshVertex() {}
  MeshVertex(const float3& p_, const float3& n_)
      : p(make_glfloat3(p_)), n(make_glfloat3(n_)) {}
  // float3 p;  // point
  // float3 n;  // normal
  GLfloat3 p;  // point
  GLfloat3 n;  // normal
};

class Mesh {
 public:
  Mesh();
  Mesh(const Mesh& rhs);
  ~Mesh();

  void operator=(const Mesh& rhs);

  void InitVBO() const;
  void Display(bool wireframe) const;
  void DisplayNormals() const;

  void AddVertex(const float3& v, float3 color = float3());
  template <typename Iter>
  void AddVertices(Iter begin, Iter end) {
    for (Iter i = begin; i != end; ++i) {
      AddVertex(*i);
    }
  }
  void SetVertexColor(const int i, const float3& color);
  void Translate(const float3& t);
  void AddTriangle(const Triangle& t);
  template <typename Iter>
  void AddTriangles(Iter begin, Iter end) {
    for (Iter i = begin; i != end; ++i) {
      AddTriangle(*i);
    }
  }
  void ClearTriangles() {
    _triangles.clear();
    _num_triangles = 0;
    _vbo_initialized = false;
  }

  void AddMaterial(const Material& m) {
    _materials.push_back(m);
    if (_cur_mtl == -1) _cur_mtl = _materials.size()-1;
  }

  float3 Centroid() const;

  const Material& GetMaterial() const { return _materials[0]; }
  void SetMaterial(const Material& m) { _materials[0] = m; }

  void new_material(int material_idx, const std::string& name) {
    _materials.push_back(Material(name));
  }

  void set_cur_material(const string& name) {
    for (int i = 0; i < _materials.size(); ++i) {
      if (_materials[i].name() == name) {
        set_cur_material(i);
      }
    }
  }

  void set_cur_material(int cur_mtl) {
    _cur_mtl = cur_mtl;
  }

  void set_ambient(int material_idx, const float3& ambient) {
    _materials[material_idx].set_ambient(ambient);
  }
  void set_diffuse(int material_idx, const float3& diffuse) {
    _materials[material_idx].set_diffuse(diffuse);
  }
  void set_specular(int material_idx, const float3& specular) {
    _materials[material_idx].set_specular(specular);
  }
  void set_specular_coeff(int material_idx, const float& coeff) {
    _materials[material_idx].set_specular_coeff(coeff);
  }
  void set_texture(int material_idx, const string& texture) {
    _materials[material_idx].set_texture(texture);
  }

  const std::vector<float3>& vertices() const {
    return _vertices;
  }
  const std::vector<float3>& colors() const {
    return _colors;
  }

  // Sets the triangles and vertices.  It is possible that there
  // are more vertices than needed, so collapse down to including only
  // those necessary.
  void SetAndCollapse(const vector<float3>& vertices,
                      const vector<Triangle>& triangles);

  const std::vector<Triangle>& triangles() const {
    return _triangles;
  }

  int num_triangles() const {
    return _num_triangles;
  }

  void OrientPolygons();
  void InvertOrientation();

  const Material& material(int i) const { return _materials[i]; }
  Material& material(int i) { return _materials[i]; }

  void SetColor(bool b) { _color = b; _vbo_initialized = false; }
  bool IsColor() const { return _color; }

  const BoundingBox3f& bb() const { return _bb; }
  int num_materials() const { return _materials.size(); }

  void compute_normals();
  const float3& normal(int i) const { return _normals[i]; }

  float pick(const float3& p, const float3& v);

  void print_stats() const;

 private:
  void Copy(const Mesh& rhs);

 private:
  // **** Update assignment operator if adding variables **** //

  std::vector<float3> _vertices;
  std::vector<float3> _normals;
  std::vector<float3> _colors;
  std::vector<Triangle> _triangles;
  std::vector<Material> _materials;
  // True if color is attached to vertices
  bool _color;
  int _cur_mtl;
  BoundingBox3f _bb;
  int _num_triangles;
  mutable bool _vbo_initialized;
  mutable GLuint _buffer_names[2];
  
  mutable float3 _centroid;
};

#endif
