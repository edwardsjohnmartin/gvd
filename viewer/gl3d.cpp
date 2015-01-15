#include "./gl3d.h"

#include <fstream>
#include <sstream>

#include "./common.h"

vector<float3> GL3D::colors;

void GL3D::InitColors() {
  if (colors.empty()) {
    ifstream in("colors.config");
    if (in) {
      while (!in.eof()) {
        float r, g, b;
        in >> r >> g >> b;
        colors.push_back(make_float3(r, g, b));
      }
      in.close();
    }
  }
}

float3 GL3D::RandomColor() {
  float3 c = make_float3(0);
  for (int i = 0; i < 3; ++i) {
    c.s[i] = random() / static_cast<float>(RAND_MAX);
  }
  return c;
}

// Returned color may not be close to avoid
float3 GL3D::RandomColor(int seed, const float3& avoid) {
  InitColors();
  if (colors.size() > seed) {
    return colors[seed];
  }
  static const float DIST_THRESH = .7;
  srandom(seed+1);
  float3 c = RandomColor();
  // while ((avoid - c).norm2() < DIST_THRESH) {
  while (length2(avoid - c) < DIST_THRESH) {
    c = RandomColor();
  }
  return c;
}

// float3 GL3D::SetColor(int i) {
//   srandom(i+1);
//   float3 c = RandomColor();
//   glColor3fv(c);
//   return c;
// }

// float3 GL3D::SetColor(int i, const float3& avoid) {
//   srandom(i+1);
//   float3 c = RandomColor(avoid);
//   glColor3fv(c);
//   return c;
// }

float3 GL3D::Win2Obj(const float3& p) const {
  GLdouble modelview[16];
  GLdouble projection[16];
  GLint viewport[4];
  GLdouble objX, objY, objZ;
  glGetDoublev(GL_MODELVIEW_MATRIX, modelview);
  glGetDoublev(GL_PROJECTION_MATRIX, projection);
  glGetIntegerv(GL_VIEWPORT, viewport);
  gluUnProject(p.s[0], p.s[1], p.s[2],
               modelview, projection, viewport,
               &objX, &objY, &objZ);
  return make_float3(objX, objY, objZ);
}

void GL3D::Reshape(const int win_width, const int win_height) {
  // window_width = w[0];
  // window_height = w[1];
  window_width = win_width;
  window_height = win_height;
  window_aspect = window_width / static_cast<float>(window_height);
}

float3 GL3D::Obj2Win(const float3& obj) const {
  GLdouble modelview[16];
  GLdouble projection[16];
  GLint viewport[4];
  GLdouble winX, winY, winZ;
  glGetDoublev(GL_MODELVIEW_MATRIX, modelview);
  glGetDoublev(GL_PROJECTION_MATRIX, projection);
  glGetIntegerv(GL_VIEWPORT, viewport);
  gluProject(obj.s[0], obj.s[1], obj.s[2],
             modelview, projection, viewport,
             &winX, &winY, &winZ);
  return make_float3(winX, winY, winZ);
}

enum Justify { kLeftJustify, kRightJustify,
               kTopJustify, kBottomJustify,
               kCenterJustify };

// p is in window coordinates
void GL3D::BitmapString(const string& s, const int2& p,
                        void* font) const {
  glDisable(GL_TEXTURE_2D);
  glDisable(GL_LIGHTING);

  // // Draw white behind
  // glColor3f(1.0f, 1.0f, 1.0f);
  // for (int i = -1; i <= 1; ++i) {
  //   for (int j = -1; j <= 1; ++j) {
  //     glRasterPos2fv(Win2Obj(Obj2Win(p) + float2(i, j)));
  //     for (int i = 0; i < s.size(); ++i) {
  //       glutBitmapCharacter(font, s[i]);
  //     }
  //   }
  // }

  glMatrixMode(GL_MODELVIEW);
  glPushMatrix();
  glLoadIdentity();

  glMatrixMode(GL_PROJECTION);
  glPushMatrix();
  glLoadIdentity();
  gluOrtho2D(0, window_width, 0, window_height);

  // Draw white behind
  glColor3f(1.0f, 1.0f, 1.0f);
  for (int i = -1; i <= 1; ++i) {
    for (int j = -1; j <= 1; ++j) {
      glRasterPos2iv((p + make_int2(i, j)).s);
      for (int i = 0; i < s.size(); ++i) {
        glutBitmapCharacter(font, s[i]);
      }
    }
  }

  // for (int i = -1; i <= 1; ++i) {
  //   for (int j = -1; j <= 1; ++j) {
  //     const float3 pos = Win2Obj(Obj2Win(p) + make_float3(i, j, -.01));
  //     glRasterPos3fv(pos.s);
  //     for (int i = 0; i < s.size(); ++i) {
  //       glutBitmapCharacter(font, s[i]);
  //     }
  //   }
  // }

  // Draw black in front
  glColor3f(0.0f, 0.0f, 0.0f);
  glRasterPos2iv(p.s);
  for (int i = 0; i < s.size(); ++i) {
    glutBitmapCharacter(font, s[i]);
  }

  glPopMatrix();
  glMatrixMode(GL_MODELVIEW);
  glPopMatrix();

  glEnable(GL_LIGHTING);
}

// p is in object coordinates
void GL3D::BitmapString(const string& s, const float3& p, void* font) const {
  glDisable(GL_TEXTURE_2D);
  glDisable(GL_LIGHTING);

  // Draw white behind
  glColor3f(1.0f, 1.0f, 1.0f);
  for (int i = -1; i <= 1; ++i) {
    for (int j = -1; j <= 1; ++j) {
      const float3 pos = Win2Obj(Obj2Win(p) + make_float3(i, j, -.01));
      glRasterPos3fv(pos.s);
      for (int i = 0; i < s.size(); ++i) {
        glutBitmapCharacter(font, s[i]);
      }
    }
  }

  // Draw black in front
  glColor3f(0.0f, 0.0f, 0.0f);
  const float3 pos = Win2Obj(Obj2Win(p) + make_float3(0, 0, -.02));
  glRasterPos3fv(pos.s);
  for (int i = 0; i < s.size(); ++i) {
    glutBitmapCharacter(font, s[i]);
  }

  glEnable(GL_LIGHTING);
}

// obj is in world coordinates.
// xoff, yoff are in window coordinates
void GL3D::BitmapString(const string& s, const float3& obj,
                        int xoff, int yoff,
                        void* font) const {
  float3 w = Obj2Win(obj);
  float3 p = Win2Obj(w + make_float3(xoff, -yoff, 0));
  BitmapString(s, p, font);
}

void GL3D::BitmapString(int value, const float3& obj,
                        int xoff, int yoff,
                        void* font) const {
  stringstream ss;
  ss << value;
  BitmapString(ss.str(), obj, xoff, yoff, font);
}

void GL3D::BitmapString(const string& s, const float3& obj,
                  Justify hjustify, Justify vjustify) const {
  float3 w = Obj2Win(obj);
  int xoff = 1, yoff = 1;
  if (hjustify == kRightJustify) {
    xoff = - 8 * s.size() - 1;
  } else if (hjustify == kCenterJustify) {
    xoff = - 4 * s.size() - 1;
  }
  if (vjustify == kTopJustify) {
    yoff = - 13 - 1;
  } else if (vjustify == kCenterJustify) {
    yoff = - 8 - 1;
  }
  float3 p = Win2Obj(w + make_float3(xoff, 1, 0));
  BitmapString(s, p, GLUT_BITMAP_8_BY_13);
}

void GL3D::BitmapString(int value, const float3& obj,
                  Justify hjustify, Justify vjustify) const {
  stringstream ss;
  ss << value;
  BitmapString(ss.str(), obj, hjustify, vjustify);
}

void GL3D::BitmapString(const string& s) const {
  glDisable(GL_TEXTURE_2D);
  glColor3f(0.0f, 0.0f, 0.0f);
  float3 p = Win2Obj(make_float3(4, window_height-5, 0));
  glRasterPos2f(p.s[0], p.s[1]);
  for (int i = 0; i < s.size(); ++i) {
    glutBitmapCharacter(GLUT_BITMAP_HELVETICA_18, s[i]);
  }
}

// size = edge length
void GL3D::glCube(const float3& obj, int size) const {
  const float3 win = Obj2Win(obj);
  const float3 obj_offset =
      Win2Obj(make_float3(win.s[0]+size, win.s[1], win.s[2]));
  const double off = length(obj_offset-obj);
  glTranslatef(obj.s[0], obj.s[1], obj.s[2]);
  glutSolidCube(off);
  glTranslatef(-obj.s[0], -obj.s[1], -obj.s[2]);
}

void GL3D::SetMaterial(const Material& m) const {
  const GLfloat mat_emission[] = { 0, 0, 0, 1 };
  const GLfloat mat_specular[] = { m.specular().s[0], m.specular().s[1],
                                   m.specular().s[2], 1 };
  const GLfloat mat_diffuse[] = { m.diffuse().s[0], m.diffuse().s[1],
                                  m.diffuse().s[2], 1 };
  const GLfloat mat_ambient[] = { m.ambient().s[0], m.ambient().s[1],
                                  m.ambient().s[2], 1 };
  const GLfloat alpha = m.specular_coeff();
  glMaterialfv(GL_FRONT_AND_BACK, GL_EMISSION, mat_emission);
  glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, mat_specular);
  glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, mat_diffuse);
  glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, mat_ambient);
  glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, alpha);
  const int texture = m.texture_id();
  if (texture != -1) {
    glEnable(GL_TEXTURE_2D);
    glBindTexture(GL_TEXTURE_2D, texture);
  } else {
    glDisable(GL_TEXTURE_2D);
  }
}

void GL3D::SetDiffuseAmbient(
    const float r, const float g, const float b) const {
  float3 c = make_float3(r, g, b);
  Material m;
  m.set_diffuse(c);
  m.set_ambient(c);
  SetMaterial(m);
}
