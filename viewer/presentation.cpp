#include <math.h>
#include <assert.h>
#include <iostream>
#include <fstream>
#include <utility>
#include <vector>
#include <stdexcept>
#include <algorithm>
#include <string>

#include "./common.h"
#include "../vec.h"
#include "./medial2.h"
#include "./medial3.h"
#include "./texture.h"
#include "../shared_ptr.h"
#include "../vector2.h"

using namespace std;

typedef BoundingBox<GLfloat, 2> BB;

int kWindowWidth = 600;
int kWindowHeight = 450;
int window_width = kWindowWidth;
int window_height = kWindowHeight;
GLuint* texture_ids;
typedef oct::shared_ptr<Scene> ScenePtr;
vector<ScenePtr> scenes;
int cur_idx = 0;

void DisplaySlide(const GLfloat x, const GLfloat y) {
  // Render the slide texture
  glDisable(GL_LIGHTING);
  glEnable(GL_TEXTURE_2D);
  const int texture = texture_ids[cur_idx];
  glBindTexture(GL_TEXTURE_2D, texture);

  const float2 tmin(0, 0);
  const float2 tmax(1, 1);
  const float2 vmin(-x, -y);
  const float2 vmax(x, y);

  glColor3f(1, 1, 1);
  glPolygonMode(GL_FRONT, GL_FILL);
  glBegin(GL_QUADS);
  glTexCoord2d(tmin[0], tmin[1]);
  glVertex2d(vmin[0], vmin[1]);
  glTexCoord2d(tmax[0], tmin[1]);
  glVertex2d(vmax[0], vmin[1]);
  glTexCoord2d(tmax[0], tmax[1]);
  glVertex2d(vmax[0], vmax[1]);
  glTexCoord2d(tmin[0], tmax[1]);
  glVertex2d(vmin[0], vmax[1]);
  glEnd();

  glDisable(GL_TEXTURE_2D);
}

void MyDisplay() {
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

  const float window_aspect = window_width/static_cast<float>(window_height);
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  // from -100 to 100 in z
  const GLfloat near = 1;
  const GLfloat far = 201;
  const GLfloat theta = 40.0;
  const GLfloat z_eye = 101;
  gluPerspective(theta, window_aspect, near, far);

  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
  gluLookAt(0, 0, z_eye,
            0, 0, 0,
            0, 1, 0);

  glTranslatef(0, 0, -(far - z_eye));
  const float y = far * tan((theta/2) * M_PI / 180);
  const float x = y * window_aspect;
  DisplaySlide(x, y);
  glLoadIdentity();

  scenes[cur_idx]->Display();

  glFlush();
  glutSwapBuffers();
}

void Init() {
  glClearColor(1.0, 1.0, 1.0, 1.0);
  scenes[cur_idx]->Init();
}

void Mouse(int button, int state, int x, int y) {
  scenes[cur_idx]->Mouse(button, state, x, y);
}

void MouseMotion(int x, int y) {
  scenes[cur_idx]->MouseMotion(x, y);
}

void PassiveMouseMotion(int x, int y) {
  scenes[cur_idx]->PassiveMouseMotion(x, y);
}

void Advance(const int inc) {
  if (cur_idx + inc >= 0 && cur_idx + inc < scenes.size()) {
    cur_idx += inc;
    scenes[cur_idx]->Init();
    glutPostRedisplay();
  }
}

void Keyboard(unsigned char key, int x, int y) {
  switch (key) {
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
      scenes[cur_idx]->Keyboard(key, x, y);
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
      scenes[cur_idx]->Special(key, x, y);
      break;
  }
}

// Window size changed
void Reshape(int width, int height) {
  cout << "Reshaping to " << width << " " << height << endl;
  window_width = width;
  window_height = height;
  glViewport(0, 0, window_width, window_height);
  // const float r = window_width / static_cast<float>(window_height);
  for (int i = 0; i < scenes.size(); ++i) {
    // scenes[i]->Reshape(int2(window_width, window_height));
    // scenes[i]->Reshape(window_width, window_height);
    scenes[i]->Reshape(window_width, window_height);
                       // BoundingBox<GLfloat, 2>(float2(0, 0), float2(500, 500)));
    // scenes[i]->Reshape(int2(window_width, window_height),
    //                    float2(-r, -1), float2(r, 1));
  }
  Init();
}

int main(int argc, char** argv) {
  glutInit(&argc, argv);
  glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
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

  const int n = argc-1;
  texture_ids = new GLuint[n];
  glGenTextures(n, texture_ids);
  for (int i = 0; i < n; ++i) {
    int w, h;
    cout << "Reading " << argv[i+1] << endl;
    png_texture_load(argv[i+1], &w, &h, texture_ids[i]);

    glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
  }

  for (int i = 0; i < n; ++i) {
    scenes.push_back(ScenePtr(new Scene()));
  }

  {
    Medial2* scene = new Medial2(window_width, window_height);
    scene->ReadMesh("../data2/test1.dat");
    scene->ReadMesh("../data2/test2.dat");
    scene->SetShowMedialAxis(false);
    scene->SetShowOctree(true);
    // scene->SetViewObj(BB(float2(-.9, -1), float2(.9, .8)));
    scene->SetViewObj(BB(float2(-1, -1), float2(1, 1)));
    scenes[0].reset(scene);
  // } {
  //   Medial2* scene = new Medial2(window_width, window_height);
  //   scene->ReadMesh("../data2/test1.dat");
  //   scene->ReadMesh("../data2/test2.dat");
  //   scene->SetShowMedialAxis(false);
  //   scene->SetShowOctree(true);
  //   scene->SetShowVertexLabel(true);
  //   scene->SetShowClosestPointLine(true);
  //   scenes[29].reset(scene);
  // } {
  //   Medial2* scene = new Medial2(window_width, window_height);
  //   scene->ReadMesh("../data2/test1.dat");
  //   scene->ReadMesh("../data2/test2.dat");
  //   scene->SetShowVertexLabel(true);
  //   scenes[30].reset(scene);
  // } {
  //   Medial3* scene = new Medial3(window_width, window_height);
  //   scene->ReadMesh("../data3/cube.obj");
  //   scene->ReadMesh("../data3/cube2.obj");
  //   scene->ReadMesh("../data3/tet.obj");
  //   oct::OctreeOptions o = oct::OctreeOptions::For3D();
  //   scene->SetShowMedial(false);
  //   scene->GenerateSurface(o);
  //   scenes[31].reset(scene);
  // // } {
  // //   Medial3* scene = new Medial3(window_width, window_height);
  // //   scene->ReadMesh("../data3/d000_tiles.obj");
  // //   scene->ReadMesh("../data3/d001_tiles.obj");
  // //   scene->ReadMesh("../data3/d002_tiles.obj");
  // //   scene->ReadMesh("../data3/d003_tiles.obj");
  // //   scene->ReadMesh("../data3/d004_tiles.obj");
  // //   scene->ReadMesh("../data3/d005_tiles.obj");
  // //   scene->ReadMesh("../data3/d006_tiles.obj");
  // //   scene->ReadMesh("../data3/d007_tiles.obj");
  // //   scene->ReadMesh("../data3/d008_tiles.obj");
  // //   scene->ReadMesh("../data3/d009_tiles.obj");
  // //   scene->ReadMesh("../data3/medial_8_d_0_9.obj", true);
  // //   // oct::OctreeOptions o = oct::OctreeOptions::For3D();
  // //   // scene->GenerateSurface(o);
  // //   scenes[32].reset(scene);
  // } {
  //   Medial3* scene = new Medial3(window_width, window_height);
  //   scene->ReadMesh("../data3/d001.obj");
  //   scene->ReadMesh("../data3/d009.obj");
  //   scene->ReadMesh("../data3/a001.obj");
  //   scene->ReadMesh("../data3/a017_tiles.obj");
  //   scene->ReadMesh("../data3/a020_tiles.obj");
  //   scene->ReadMesh("../data3/medial_8_2d_3a.obj", true);
  //   // oct::OctreeOptions o = oct::OctreeOptions::For3D();
  //   // scene->GenerateSurface(o);
  //   scene->SetShowMedial(false);
  //   scenes[32].reset(scene);
  // } {
  //   Medial3* scene = new Medial3(window_width, window_height);
  //   scene->ReadMesh("../data3/d001.obj");
  //   scene->ReadMesh("../data3/d009.obj");
  //   scene->ReadMesh("../data3/a001.obj");
  //   scene->ReadMesh("../data3/a017_tiles.obj");
  //   scene->ReadMesh("../data3/a020_tiles.obj");
  //   scene->ReadMesh("../data3/medial_10_2d_3a.obj", true);
  //   // oct::OctreeOptions o = oct::OctreeOptions::For3D();
  //   // scene->GenerateSurface(o);
  //   scenes[33].reset(scene);
  }

  Init();
  glutMainLoop();
}
