#include <math.h>
#include <assert.h>
#include <iostream>
#include <fstream>
#include <utility>
#include <vector>
#include <stdexcept>
#include <algorithm>
#include <string>

#include "./medial2.h"
#include "../shared_ptr.h"
#include "./pngio.h"

using namespace std;

int kWindowWidth = 600;
int kWindowHeight = 450;
int window_width = kWindowWidth;
int window_height = kWindowHeight;
oct::shared_ptr<GL2D> scene;
int cur_idx = 0;

void MyDisplay() {
  glClear(GL_COLOR_BUFFER_BIT);

  // Not sure why this doesn't work
  const float r = window_width / static_cast<float>(window_height);
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  gluOrtho2D(-r, r, -1, 1);

  // const float window_aspect = window_width/static_cast<float>(window_height);
  // glMatrixMode(GL_PROJECTION);
  // glLoadIdentity();
  // // from -100 to 100 in z
  // const GLfloat near = 1;
  // const GLfloat far = 201;
  // const GLfloat theta = 40.0;
  // const GLfloat z_eye = 101;
  // gluPerspective(theta, window_aspect, near, far);

  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();

  scene->Display();
  glFlush();
  glutSwapBuffers();
}

void Init() {
  glClearColor(1.0, 1.0, 1.0, 1.0);
  scene->Init();

  // const float r = window_width / static_cast<float>(window_height);
  // glMatrixMode(GL_PROJECTION);
  // glLoadIdentity();
  // gluOrtho2D(-r, r, -1, 1);
  // glMatrixMode(GL_MODELVIEW);
  // glLoadIdentity();
}

void Mouse(int button, int state, int x, int y) {
  scene->Mouse(button, state, x, y);
}

void MouseMotion(int x, int y) {
  scene->MouseMotion(x, y);
}

void PassiveMouseMotion(int x, int y) {
  scene->PassiveMouseMotion(x, y);
}

void Advance(const int inc) {
  cur_idx += inc;
  glutPostRedisplay();
}

void Keyboard(unsigned char key, int x, int y) {
  switch (key) {
    case 'Y': {
      GVDViewer2* gvd_viewer = (GVDViewer2*)scene.get();
      gvd_viewer->SetShowStatistics(false);
      // 30 fps is iMovie max
      const int fps = 30;
      const int secs = 30;
      const int n = fps*secs;
      const float max_zoom = 80000;
      const int n_size = 4;
      const float2 target = make_float2(0.0687421, 0.024298);
      gvd_viewer->Zoom(target, 1);
      for (int i = 0; i < n; ++i) {
        // Update zoom level
        // When i = 0, then zoom == 1
        // When i = n-1, then zoom == max_zoom.
        const float f = (i/(float)n);
        const float zoom = pow(max_zoom, f);
        gvd_viewer->Zoom(target, zoom);
        MyDisplay();

        stringstream ss;
        ss << "shot-" << setfill('0') << setw(n_size) << (i+1);
        string fn = ss.str() + ".png";
        writePngImage(fn.c_str(), window_width, window_height);

        // stringstream ss;
        // ss << "sshot-";
        // stringstream tmp2;
        // tmp2 << (i+1);
        // const string istr = tmp2.str();
        // for (int j = istr.size(); j < nstr.size(); ++j) {
        //   ss << "0";
        // }
        // ss << istr << ".png";

        // writePngImage(ss.str().c_str(), window_width, window_height);
      }
      break;
      }
    case ' ':
      Advance(1);
      break;
    case 12:
      if (glutGetModifiers() & GLUT_ACTIVE_CTRL) {
        glutFullScreen();
        glutPostRedisplay();
      }
      break;
    case 27:
      glutReshapeWindow(kWindowWidth, kWindowHeight);
      glutPostRedisplay();
      break;
    case 'q':
      exit(EXIT_SUCCESS);
      break;
    default:
      scene->Keyboard(key, x, y);
      break;
  }
}

void Special(int key, int x, int y) {
  switch (key) {
    case GLUT_KEY_DOWN:
    case GLUT_KEY_RIGHT:
      Advance(1);
      break;
    case GLUT_KEY_UP:
    case GLUT_KEY_LEFT:
      Advance(-1);
      break;
    case GLUT_KEY_F5:
      glutFullScreen();
      glutPostRedisplay();
      break;
    default:
      scene->Special(key, x, y);
      break;
  }
}

// Window size changed
void Reshape(int width, int height) {
  window_width = width;
  window_height = height;
  glViewport(0, 0, window_width, window_height);
  const float r = window_width / static_cast<float>(window_height);
  // scene->Reshape(int2(window_width, window_height),
  //               float2(-r, -1), float2(r, 1));
  scene->Reshape(window_width, window_height);
  Init();
}

int main(int argc, char** argv) {
  ifstream sizein("size.config");
  if (sizein) {
    sizein >> window_width >> window_height;
    sizein.close();
  }

  const float r = window_width / static_cast<float>(window_height);
  // scene.reset(new GVDViewer2(int2(window_width, window_height),
  //                         float2(-r, -1), float2(r, 1)));
  scene.reset(new GVDViewer2(window_width, window_height)),

  glutInit(&argc, argv);
  glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB);
  glutInitWindowSize(window_width, window_height);
  glutInitWindowPosition(0, 0);
  glutCreateWindow("octree");
  glutDisplayFunc(MyDisplay);
  glutKeyboardFunc(Keyboard);
  glutSpecialFunc(Special);
  glutMouseFunc(Mouse);
  glutMotionFunc(MouseMotion);
  glutPassiveMotionFunc(PassiveMouseMotion);
  glutReshapeFunc(Reshape);

  scene->ProcessArgs(argc, argv);

  Init();
  glutMainLoop();
}
