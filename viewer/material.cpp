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

#include <string>
#include "./material.h"
#include "./texture.h"

using namespace std;

Material::Material()
    : _ambient(make_float3(.3, .3, .3)),
      _diffuse(make_float3(.7, .7, .7)),
      _specular(make_float3(.1, .1, .1)), _specular_coeff(10),
      _texture_id(-1) {
}

Material::Material(const string& name)
    : _name(name),
      _ambient(make_float3(.1, .1, .1)),
      _diffuse(make_float3(.3, .3, .3)),
      _specular(make_float3(.1, .1, .1)), _specular_coeff(10),
      _texture_id(-1) {
}

Material Material::FromDiffuseAmbient(const float3& diff_amb, float coeff) {
  Material m;
  m.set_diffuse(diff_amb);
  m.set_ambient(diff_amb);
  // m.set_specular(diff_amb);
  m.set_specular(make_float3(1, 1, 1));
  m.set_specular_coeff(coeff);
  return m;
}

void Material::set_ambient(const float3& ambient) {
  _ambient = ambient;
}
void Material::set_diffuse(const float3& diffuse) {
  _diffuse = diffuse;
}
void Material::set_specular(const float3& specular) {
  _specular = specular;
}
void Material::set_specular_coeff(const float& coeff) {
  _specular_coeff = coeff;
}
void Material::set_texture(const string& texture) {
  _texture = texture;
}

void Material::LoadTexture(int texture_id) {
  ::LoadTexture(texture(), texture_id);
  _texture_id = texture_id;
}
