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

#ifndef __GL_3D_H__
#define __GL_3D_H__

#include <string>

#include "../bb.h"
#include "./common.h"
#include "./material.h"
#include "./scene.h"

class GL3D : public Scene {
 public:
  enum Justify { kLeftJustify, kRightJustify,
                 kTopJustify, kBottomJustify,
                 kCenterJustify };

 public:
  GL3D(const int win_width, const int win_height) {//, const float2& world_min, const float2& world_max) {
    Reshape(win_width, win_height);//, world_min, world_max);
  }
  virtual ~GL3D() {}

  virtual void Reshape(const int win_width, const int win_height);
                       // const float2& world_min, const float2& world_max);

  virtual int ProcessArgs(int argc, char** argv) { return 0; }
  virtual void Init() {}
  virtual void Mouse(int button, int state, int x, int y) {}
  virtual void MouseMotion(int x, int y) {}
  virtual void PassiveMouseMotion(int x, int y) {}
  virtual void Keyboard(unsigned char key, int x, int y) {}
  virtual void Special(unsigned char key, int x, int y) {}
  virtual void Display() {}

  static void InitColors();
  static float3 RandomColor();
  // Returned color may not be close to avoid
  static float3 RandomColor(int seed, const float3& avoid);

  float3 Win2Obj(const float3& p) const;
  float3 Obj2Win(const float3& obj) const;

  // p is in window coordinates.
  void BitmapString(const string& s, const int2& p,
                          void* font = GLUT_BITMAP_8_BY_13) const;

  // p is in object coordinates
  void BitmapString(const string& s, const float3& p,
                    void* font = GLUT_BITMAP_8_BY_13) const;

// obj is in world coordinates.
// xoff, yoff are in window coordinates
  void BitmapString(const string& s, const float3& obj,
                  int xoff, int yoff,
                    void* font = GLUT_BITMAP_8_BY_13) const;

  void BitmapString(int value, const float3& obj,
                  int xoff = 1, int yoff = 1,
                    void* font = GLUT_BITMAP_8_BY_13) const;

  void BitmapString(const string& s, const float3& obj,
                    Justify hjustify, Justify vjustify) const;

  void BitmapString(int value, const float3& obj,
                    Justify hjustify, Justify vjustify) const;

  void BitmapString(const string& s) const;

// size = edge length
  void glCube(const float3& obj, int size) const;
  void SetMaterial(const Material& m) const;
  void SetDiffuseAmbient(const float r, const float g, const float b) const;

 protected:
  int window_width;
  int window_height;
  float window_aspect;
  static std::vector<float3> colors;
};

#endif
