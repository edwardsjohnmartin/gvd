#ifndef __MATERIAL_H__
#define __MATERIAL_H__

#include <string>
#include <vector>
#include "../opencl/vec.h"

class Material {
 public:
  Material();
  explicit Material(const std::string& name);

  static Material FromDiffuseAmbient(const float3& diff_amb,
                                     float coeff = 40);

  void set_ambient(const float3& ambient);
  void set_diffuse(const float3& diffuse);
  void set_specular(const float3& specular);
  void set_specular_coeff(const float& coeff);
  void set_texture(const std::string& texture);
  void LoadTexture(int texture_id);

  const std::string& name() const { return _name; }
  const float3& ambient() const { return _ambient; }
  const float3& diffuse() const { return _diffuse; }
  const float3& specular() const { return _specular; }
  const float& specular_coeff() const { return _specular_coeff; }
  const std::string& texture() const { return _texture; }
  int texture_id() const { return _texture_id; }

  friend std::ostream& operator<<(std::ostream& out, const Material& m) {
    out << "[" << m._name << ", ambient=" << m._ambient << ", diffuse="
        << m._diffuse << ", specular=" << m._specular << "]";
    return out;
  }

 private:
  std::string _name;
  float3 _ambient;
  float3 _diffuse;
  float3 _specular;
  float _specular_coeff;
  std::string _texture;
  int _texture_id;
};

#endif
