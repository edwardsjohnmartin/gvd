#include "./gl2d.h"

#include <sstream>
#include <string>

using namespace std;

void GL2D::Init() const {
}

void GL2D::Reshape(const int win_width, const int win_height,
                   const BB& viewport) {
  window_width = win_width;
  window_height = win_height;
  fwindow_width = static_cast<GLfloat>(window_width);
  fwindow_height = static_cast<GLfloat>(window_height);
  window_aspect = fwindow_width / fwindow_height;
  const float2 win_min_obj = make_float2(-window_aspect, -1);
  const float2 win_max_obj = make_float2(window_aspect, 1);
  win_obj = BB(win_min_obj, win_max_obj);

  // const float2 world_min(-window_aspect
  //                       + 2*window_aspect*(viewport.min()[0]/window_width),
  //                       -1 + 2*viewport.min()[1]/window_height);
  // const float2 world_max(-window_aspect
  //                       + 2*window_aspect*(viewport.max()[0]/window_width),
  //                       -1 + 2*viewport.max()[1]/window_height);
  // view_obj = BB(world_min, world_max);
  Init();
}

void GL2D::Reshape(const int win_width, const int win_height) {
  Reshape(win_width, win_height, BB(make_float2(0, 0),
                                    make_float2(win_width, win_height)));
  // window_width = win_width;
  // window_height = win_height;
  // fwindow_width = static_cast<GLfloat>(window_width);
  // fwindow_height = static_cast<GLfloat>(window_height);
  // window_aspect = fwindow_width / fwindow_height;
  // const float2 world_min(-window_aspect, -1);
  // const float2 world_max(window_aspect, 1);
  // win_obj = BB(world_min, world_max);
  // Init();
}

float2 GL2D::Win2Obj(const float2& w) const {
  // const float2& mi = win_obj.min();
  // const float2& ma = win_obj.max();
  // const GLfloat obj_width = (ma[0]-mi[0]);
  // const GLfloat obj_height = (ma[1]-mi[1]);
  // float2 winobj(
  //     (w[0] / fwindow_width) * (obj_width) + mi[0],
  //     ((window_height-w[1]) / fwindow_height) * obj_height + mi[1]);

  // const GLfloat view_width = view_obj.size()[0];
  // const GLfloat view_height = view_obj.size()[1];
  // return float2(view_width*(winobj[0]-mi[0])/obj_width + view_obj.min()[0],
  //              view_height*(winobj[1]-mi[1])/obj_height + view_obj.min()[1]);
  const float2& mi = view_obj.min();
  const float2& ma = view_obj.max();

  const GLfloat obj_width = win_obj.size().s[0];
  const GLfloat obj_height = win_obj.size().s[1];
  return make_float2(
      (w.s[0] / fwindow_width) * (obj_width) + win_obj.min().s[0],
      ((window_height-w.s[1]) / fwindow_height) * obj_height + win_obj.min().s[1]);
}

float GL2D::Win2Obj(const float w) const {
  const GLfloat obj_width = win_obj.size().s[0];
  return (w / fwindow_width) * (obj_width);
}

float2 GL2D::Obj2Win(const float2& v) const {
  const GLfloat obj_width = win_obj.size().s[0];
  const GLfloat obj_height = win_obj.size().s[1];
  const float2 offset = v - win_obj.min();
  return make_float2(
      fwindow_width * offset.s[0]/obj_width,
      fwindow_height * (1.0 - offset.s[1]/obj_height));
}

float2 GL2D::Obj2Win(const float x, const float y) const {
  return Obj2Win(make_float2(x, y));
}

void GL2D::BitmapString(const string& s, const float2& p,
                        void* font) const {
  // Draw white behind
  glColor3f(1.0f, 1.0f, 1.0f);
  for (int i = -1; i <= 1; ++i) {
    for (int j = -1; j <= 1; ++j) {
      glRasterPos2fv(Win2Obj(Obj2Win(p) + make_float2(i, j)).s);
      for (int i = 0; i < s.size(); ++i) {
        glutBitmapCharacter(font, s[i]);
      }
    }
  }

  // Draw black in front
  glColor3f(0.0f, 0.0f, 0.0f);
  glRasterPos2fv(Win2Obj(Obj2Win(p) + make_float2(0, 0)).s);
  for (int i = 0; i < s.size(); ++i) {
    glutBitmapCharacter(font, s[i]);
  }
}

void GL2D::BitmapString(const string& s, const float2& obj,
                  int xoff, int yoff,
                  void* font) const {
  glDisable(GL_TEXTURE_2D);
  glColor3f(0.0f, 0.0f, 0.0f);
  float2 w = Obj2Win(obj);
  BitmapString(s, Win2Obj(make_float2(w.s[0]+xoff, w.s[1]-yoff)), font);
}

void GL2D::BitmapString(int value, const float2& obj,
                  int xoff, int yoff,
                  void* font) const {
  stringstream ss;
  ss << value;
  BitmapString(ss.str(), obj, xoff, yoff, font);
}

void GL2D::BitmapString(const string& s, const float2& obj,
                  Justify hjustify, Justify vjustify) const {
  glDisable(GL_TEXTURE_2D);
  glColor3f(0.0f, 0.0f, 0.0f);
  float2 w = Obj2Win(obj);
  int xoff = 1, yoff = 1;
  if (hjustify == kRightJustify) {
    xoff = - 8 * s.size() - 1;
  } else if (hjustify == kCenterJustify) {
    xoff = - 4 * s.size() - 1;
  }
  if (vjustify == kTopJustify) {
    yoff = - 13 - 1;
  } else if (vjustify == kCenterJustify) {
    // yoff = - 8 - 1;
    yoff = - 6;
  }
  float2 p = Win2Obj(make_float2(w.s[0]+xoff, w.s[1]-yoff));
  BitmapString(s, p, GLUT_BITMAP_8_BY_13);
}

void GL2D::BitmapString(int value, const float2& obj,
                  Justify hjustify, Justify vjustify) const {
  stringstream ss;
  ss << value;
  BitmapString(ss.str(), obj, hjustify, vjustify);
}

void GL2D::BitmapString(const string& s) const {
  glDisable(GL_TEXTURE_2D);
  glColor3f(0.0f, 0.0f, 0.0f);
  float2 p = Win2Obj(make_float2(4, window_height-5));
  BitmapString(s, p, GLUT_BITMAP_8_BY_13);
}

float3 GL2D::RandomColor() const {
  float3 c = make_float3(0);
  for (int i = 0; i < 3; ++i) {
    c.s[i] = random() / static_cast<float>(RAND_MAX);
  }
  return c;
}

// Returned color may not be close to avoid
float3 GL2D::RandomColor(const float3& avoid) const {
  static const float DIST_THRESH = .7;
  float3 c = RandomColor();
  // while ((avoid - c).norm2() < DIST_THRESH) {
  while (length2(avoid - c) < DIST_THRESH) {
    c = RandomColor();
  }
  return c;
}

// Returned color may not be close to avoid
float3 GL2D::RandomColor(int seed, const float3& avoid) const {
  srandom(seed+1);
  return RandomColor(avoid);
}

void GL2D::SetColor(int i, const float3& avoid) const {
  // srandom(i+1);
  // float3 c = RandomColor(avoid);
  float3 c = RandomColor(i, avoid);
  glColor3fv(c.s);
}

