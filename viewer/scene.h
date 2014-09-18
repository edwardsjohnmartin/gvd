#ifndef __GL_SCENE_H__
#define __GL_SCENE_H__

#include <string>

#include "../opencl/vec.h"

class Scene {
 private:
  typedef BoundingBox<float2> BB2;

 public:
  Scene() {
  }
  // virtual ~Scene() {}

  virtual void Reshape(const int win_width, const int win_height) {}
  virtual void Reshape(const int win_width, const int win_height,
                       const BB2& viewport) {}
  // virtual void Reshape(const int2& w,
  //                      const float2& world_min, const float2& world_max) {}

  virtual void Init() {}
  virtual void Mouse(int button, int state, int x, int y) {}
  virtual void MouseMotion(int x, int y) {}
  virtual void PassiveMouseMotion(int x, int y) {}
  virtual void Keyboard(unsigned char key, int x, int y) {}
  virtual void Special(unsigned char key, int x, int y) {}
  virtual void Display() {}

 private:
  // int window_width;
  // int window_height;
  // GLfloat fwindow_width;
  // GLfloat fwindow_height;
  // float window_aspect;
  // BB win_obj;
};

#endif
