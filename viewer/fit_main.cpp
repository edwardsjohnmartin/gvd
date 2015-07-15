#include <math.h>
#include <assert.h>
#include <iostream>
#include <fstream>
#include <utility>
#include <vector>
#include <stdexcept>
#include <algorithm>
#include <string>
#include <sstream>

#include "./common.h"
#include "../opencl/vec.h"
#include "../opencl/segment.h"
#include "../opencl/geom.h"

using namespace std;

int window_width = 600;
int window_height = 600;
GLfloat obj_left = 0;
GLfloat obj_right = 1;
GLfloat obj_bottom = 0;
GLfloat obj_top = 1;

// float_seg a(make_float2(0.2, 0.2), make_float2(0.8, 0.9));
// float_seg b(make_float2(0.1, 0.1), make_float2(0.9, 0.4));

// float_seg a(make_float2(0.5, 0.6), make_float2(0.3, 0.8));
// float_seg b(make_float2(0.2, 0.5), make_float2(0.8, 0.5));

// float_seg a(make_float2(0.2, 0.2), make_float2(0.8, 0.8));
// float_seg b(make_float2(0.2, 0.5), make_float2(0.8, 0.5));

// float_seg a(make_float2(0.2, 0.2), make_float2(0.8, 0.8));
// float_seg b(make_float2(0.5, 0.2), make_float2(0.5, 0.8));

// float_seg a(make_float2(0.5, 0.6), make_float2(0.3, 0.8));
// float_seg b(make_float2(0.2, 0.3), make_float2(0.8, 0.5));

// float_seg a(make_float2(0.2, 0.2), make_float2(0.8, 0.8));
// float_seg b(make_float2(0.3, 0.3), make_float2(0.7, 0.7));

float_seg a(make_float2(17136, 12798), make_float2(17133, 12800));
float_seg b(make_float2(17136, 12799), make_float2(17135, 12800));

float2 Win2Obj(const int x, const int y) {
  static GLfloat obj_width = obj_right - obj_left;
  static GLfloat obj_height = obj_top - obj_bottom;
  static GLfloat fwindow_width = static_cast<GLfloat>(window_width);
  static GLfloat fwindow_height = static_cast<GLfloat>(window_height);
  return make_float2(
      (x / fwindow_width) * (obj_width) + obj_left,
      ((window_height-y) / fwindow_height) * (obj_height) + obj_bottom);
}

int2 Obj2Win(const float x, const float y) {
  static GLfloat obj_width = obj_right - obj_left;
  static GLfloat obj_height = obj_top - obj_bottom;
  static GLfloat fwindow_width = static_cast<GLfloat>(window_width);
  static GLfloat fwindow_height = static_cast<GLfloat>(window_height);
  return make_int2(
      static_cast<int>((fwindow_width * (x - obj_left)) / obj_width),
      static_cast<int>(fwindow_height * (1.0 - (y - obj_bottom) / obj_height)));
}

int2 Obj2Win(const float2& v) {
  return Obj2Win(v.s[0], v.s[1]);
}

enum Justify { kLeftJustify, kRightJustify,
               kTopJustify, kBottomJustify,
               kCenterJustify };

// buf is in window coordinates
// void BitmapString(const string& s, float objx, float objy, int buf = 1,
//                   void* font = GLUT_BITMAP_HELVETICA_12) {
void BitmapString(const string& s, float objx, float objy,
                  int xoff = 1, int yoff = 1,
                  void* font = GLUT_BITMAP_8_BY_13) {
  glDisable(GL_TEXTURE_2D);
  glColor3f(0.0f, 0.0f, 0.0f);
  int2 w = Obj2Win(objx, objy);
  float2 p = Win2Obj(w.s[0]+xoff, w.s[1]-yoff);
  glRasterPos2f(p.s[0], p.s[1]);
  for (int i = 0; i < s.size(); ++i) {
    glutBitmapCharacter(font, s[i]);
  }
}

void BitmapString(int value, float objx, float objy,
                  int xoff = 1, int yoff = 1,
                  void* font = GLUT_BITMAP_8_BY_13) {
  stringstream ss;
  ss << value;
  BitmapString(ss.str(), objx, objy, xoff, yoff, font);
}

void BitmapString(const string& s, float objx, float objy,
                  Justify hjustify, Justify vjustify) {
  glDisable(GL_TEXTURE_2D);
  glColor3f(0.0f, 0.0f, 0.0f);
  int2 w = Obj2Win(objx, objy);
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
  float2 p = Win2Obj(w.s[0]+xoff, w.s[1]-yoff);
  glRasterPos2f(p.s[0], p.s[1]);
  for (int i = 0; i < s.size(); ++i) {
    glutBitmapCharacter(GLUT_BITMAP_8_BY_13, s[i]);
  }
}

void BitmapString(int value, float objx, float objy,
                  Justify hjustify, Justify vjustify) {
  stringstream ss;
  ss << value;
  BitmapString(ss.str(), objx, objy, hjustify, vjustify);
}

void BitmapString(const string& s) {
  glDisable(GL_TEXTURE_2D);
  glColor3f(0.0f, 0.0f, 0.0f);
  float2 p = Win2Obj(4, window_height-5);
  glRasterPos2f(p.s[0], p.s[1]);
  for (int i = 0; i < s.size(); ++i) {
    glutBitmapCharacter(GLUT_BITMAP_HELVETICA_18, s[i]);
  }
}

void Init() {
  glClearColor(1.0, 1.0, 1.0, 1.0);

  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  gluOrtho2D(obj_left, obj_right, obj_bottom, obj_top);
  glMatrixMode(GL_MODELVIEW);
}

enum kCircleType { DOTTED, SOLID, FILLED };

void glCircle(const float2& c, const double r, const kCircleType& t) {
  static const double M_2PI = 2 * M_PI;
  const double inc = M_2PI / 64;
  if (t == DOTTED)
    glBegin(GL_LINES);
  else if (t == SOLID)
    glBegin(GL_LINE_STRIP);
  else {
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
    glBegin(GL_POLYGON);
  }
  for (double d = 0; d < M_2PI; d += inc) {
    glVertex2f(c.s[0] + r*cos(d), c.s[1] + r*sin(d));
  }
  glVertex2f(c.s[0] + r, c.s[1]);
  glEnd();
}

// ra = inner radius
// rb = outer radius
void glAnnulus(const float2& c, const double ra, const double rb) {
  static const double M_2PI = 2 * M_PI;
  const double inc = M_2PI / 32;

  glBegin(GL_LINES);
  for (double d = 0; d < M_2PI; d += inc) {
    glVertex2f(c.s[0] + ra*cos(d), c.s[1] + ra*sin(d));
  }
  glVertex2f(c.s[0] + ra, c.s[1]);
  glEnd();

  glBegin(GL_LINE_STRIP);
  for (double d = 0; d < M_2PI; d += inc) {
    glVertex2f(c.s[0] + rb*cos(d), c.s[1] + rb*sin(d));
  }
  glVertex2f(c.s[0] + rb, c.s[1]);
  glEnd();
}

void glQuad(const float2& v, float size) {
  glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
  glBegin(GL_QUADS);
  glVertex2fv(v.s);
  glVertex2fv((v + make_float2(size, 0)).s);
  glVertex2fv((v + make_float2(size, size)).s);
  glVertex2fv((v + make_float2(0, size)).s);
  glEnd();
}

//      \__     2    |    1     __/
//         \__       |       __/
//            \__    |    __/
//       3       \__ | __/       0
//  ________________\|/________________
//                __/|\__     
//       4     __/   |   \__     7
//          __/      |      \__
//       __/    5    |    6    \__
//     _/            |            \_
int angle_octant(const float2& v) {
  const float theta = atan2(v.y, v.x);
  return ((theta / M_PI) * 4);
}

void DrawSeparator(float_seg A, float_seg B, const float& min_d) {
  vector<floatn> origins;
  vector<float> lengths;
  if (!multi_intersection(A, B)) {
    FitBoxes(A, B, min_d, &origins, &lengths);
  }

  glColor3f(0, 1, 0);
  for (int i = 0; i < origins.size(); ++i) {
    const floatn& o = origins[i];
    const float& d = lengths[i];
    glQuad(o, d/2);
    glQuad(o + make_float2(d/2, 0), d/2);
    glQuad(o + make_float2(0, d/2), d/2);
    glQuad(o + make_float2(d/2, d/2), d/2);
  }

  glColor3f(1, 0, 0);
  glPointSize(4);
  glBegin(GL_POINTS);
  for (int i = 0; i < origins.size(); ++i) {
    const floatn& o = origins[i];
    const float& d = lengths[i];
    glVertex2fv((o+make_uni_floatn(d/2)).s);
  }
  glEnd();
}

void Display() {
  glClear(GL_COLOR_BUFFER_BIT);

  // draw quadtree
  glColor3f(0, 0, 0);
  // glQuad(a.a(), 0.4);

  glBegin(GL_LINES);
  glVertex2fv(a.a().s);
  glVertex2fv(a.b().s);
  glVertex2fv(b.a().s);
  glVertex2fv(b.b().s);
  glEnd();

  DrawSeparator(a, b, 0.01);

  glFlush();
  glutSwapBuffers();
}

bool mod_a;
void Mouse(int button, int state, int x, int y) {
  if (button == GLUT_LEFT_BUTTON) {
    if (state == GLUT_DOWN) {
      mod_a = !(glutGetModifiers() & GLUT_ACTIVE_SHIFT);
      if (mod_a) {
        a.a() = Win2Obj(x, y);
      } else {
        b.a() = Win2Obj(x, y);
      }
      glutPostRedisplay();
    }
  }
}

void MouseMotion(int x, int y) {
  if (mod_a) {
    a.b() = Win2Obj(x, y);
  } else {
    b.b() = Win2Obj(x, y);
  }
  glutPostRedisplay();
}

void Keyboard(unsigned char key, int x, int y) {
  int i;
  float da;
  switch (key) {
    case 'a':
      break;
    case 'q':
      exit(EXIT_SUCCESS);
      break;
  }
  glutPostRedisplay();
}

int main(int argc, char** argv) {
  bool i1, i2;
  const floatn i = line_intersection(
      make_floatn(17136, 12798),
      make_floatn(17133, 12800),
      make_floatn(17136, 12799),
      make_floatn(17135, 12800),
      &i1, &i2);
  cout << i << " " << i1 << " " << i2 << endl;
  for (int j = 17000; j > 0; j-=1) {
    cout << endl;
    const floatn i = line_intersection(
        // make_floatn(17136, 12798),
        // make_floatn(17133, 12800),
        // make_floatn(17136, 12799),
        // make_floatn(17135, 12800),
        make_floatn(j+136, 12798),
        make_floatn(j+133, 12800),
        make_floatn(j+136, 12799),
        make_floatn(j+135, 12800),
        &i1, &i2);
    if (!i2) {
      cout << "j = " << j << ": " << i << " " << i1 << " " << i2 << endl;
      break;
    }
  }
  exit(0);

  cout << endl;
  cout << "Key commands:" << endl;
  cout << "  c - clear" << endl;
  cout << endl;

  glutInit(&argc, argv);
  glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB);
  glutInitWindowSize(window_width, window_height);
  glutInitWindowPosition(0, 0);
  glutCreateWindow("fit");
  glutDisplayFunc(Display);
  glutKeyboardFunc(Keyboard);
  glutMouseFunc(Mouse);
  glutMotionFunc(MouseMotion);
  Init();

  glutMainLoop();
}
