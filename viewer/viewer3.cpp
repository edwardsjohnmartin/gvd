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

#include "./gvdviewer3.h"
#include "../opencl.h"
#include "../orientation.h"
#include "./pngio.h"
#include "../vector3.h"
#include "./io.h"
#define _USE_MATH_DEFINES
#include <math.h>

using namespace std;

void ExplodeKnives();
void ExplodeRiceDwarf();
void Rotate360();
void Screenshot();
void SimpleScreenshot();

int window_width = 800, window_height = 600;
float window_aspect = window_width / static_cast<float>(window_height);

GVDViewer3 scene(window_width, window_height);
// GLfloat near = 1;
GLfloat near_ = -1;
GLfloat far_ = -1;

void MyDisplay() {
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

  const float window_aspect = window_width/static_cast<float>(window_height);
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  // from -100 to 100 in z
  // const GLfloat far = 201;
  const GLfloat theta = 40.0;
  if (near_ == -1) {
    BoundingBox3f bb = scene.bbox_full();
    const float zoom = 1;
    const float3 size = bb.size();
    const float3 eye = make_float3(
        0, 0, size.s[2]/2 + ((size.s[1]/2)/tan(20*M_PI/180.0))*1.1) * zoom;
    near_ = 1;
    far_ = length(eye) + length(size)*2;
  }
  gluPerspective(theta, window_aspect, near_, far_);

  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();

  scene.Display();
  glFlush();
  glutSwapBuffers();
}

void Init() {
  glClearColor(1.0, 1.0, 1.0, 1.0);
  scene.Init();
}

void Mouse(int button, int state, int x, int y) {
  scene.Mouse(button, state, x, y);
}

void MouseMotion(int x, int y) {
  scene.MouseMotion(x, y);
}

void Keyboard(unsigned char key, int x, int y) {
  //const float zoom_factor = 1.2;
  switch (key) {
    case 'c':
      near_ *= 0.9;
      glutPostRedisplay();
      break;
    case 'C':
      near_ *= 1.1;
      glutPostRedisplay();
      break;
    case 'u':
      far_ *= 0.9;
      glutPostRedisplay();
      break;
    case 'U':
      far_ *= 1.1;
      glutPostRedisplay();
      break;
    case 's':
      SimpleScreenshot();
      break;
    case 'S':
      Screenshot();
      break;
    case 'y': 
      Rotate360();
      break;
    case 'Y': 
      ExplodeRiceDwarf();
      break;
    case 'q':
    case 27:  // esc
      exit(0);
      break;
    default:
      scene.Keyboard(key, x, y);
      break;
  }
}

void Special(int key, int x, int y) {
  scene.Special(key, x, y);
}

void SimpleScreenshot() {
  // scene.WriteGvdMesh();
  scene.SetShowStatistics(false);
  /* function undefined. */
  //writePngImage("screenshot.png", window_width, window_height);
}

// To rename:
//    for f in *.png; do mv $f new-name${f/screenshot}; done
void Screenshot() {
  const int mm = scene.GetGVDMode();

  scene.WriteGvdMesh();

  scene.SetShowStatistics(false);
  /* function undefined. */
  //writePngImage("screenshot.png", window_width, window_height);

  // objects
  scene.SetShowMesh(true);
  scene.SetShowGVD(false);
  MyDisplay();
  /* function undefined. */
  //writePngImage("screenshot-o.png", window_width, window_height);

  // objects and GVD
  scene.SetShowMesh(true);
  scene.SetShowGVD(true);
  scene.SetGVDColor(false);
  scene.SetGVDMode(3);
  MyDisplay();
  /* function undefined. */
  //writePngImage("screenshot-om.png", window_width, window_height);

  // objects and inverted GVD
  scene.SetShowMesh(true);
  scene.SetShowGVD(true);
  scene.SetGVDColor(false);
  scene.SetGVDMode(3);
  scene.InvertGVD();
  MyDisplay();
  scene.InvertGVD();
  /* function undefined. */
  //writePngImage("screenshot-omi.png", window_width, window_height);

  // objects and colored GVD
  scene.SetShowMesh(true);
  scene.SetShowGVD(true);
  scene.SetGVDColor(true);
  scene.SetGVDMode(2);
  MyDisplay();
  /* function undefined. */
  //writePngImage("screenshot-omc.png", window_width, window_height);

  // GVD
  scene.SetShowMesh(false);
  scene.SetShowGVD(true);
  scene.SetGVDColor(false);
  scene.SetGVDMode(3);
  MyDisplay();
  /* function undefined. */
  //writePngImage("screenshot-m.png", window_width, window_height);

  // inverted GVD
  scene.SetShowMesh(false);
  scene.SetShowGVD(true);
  scene.SetGVDColor(false);
  scene.SetGVDMode(3);
  scene.InvertGVD();
  MyDisplay();
  scene.InvertGVD();
  /* function undefined. */
  //writePngImage("screenshot-mi.png", window_width, window_height);

  // colored GVD
  scene.SetShowMesh(false);
  scene.SetShowGVD(true);
  scene.SetGVDColor(true);
  scene.SetGVDMode(2);
  MyDisplay();
  /* function undefined. */
  //writePngImage("screenshot-mc.png", window_width, window_height);

  // Restore
  scene.SetShowMesh(true);
  scene.SetShowGVD(true);
  scene.SetGVDColor(false);
  scene.SetGVDMode(mm);
  MyDisplay();
}

// Writes images to shot-xxxx.png
void Rotate360() {
  scene.RotVec() = make_float3(0, 1, 0);
  scene.SetShowStatistics(false);
  const int fps = 30;
  const int secs_per_360 = 10;
  const int n = fps*secs_per_360;
  scene.RotAngle() = 0;
  const int n_size = 4;
  for (int i = 0; i < n; ++i) {
    scene.RotAngle() = (i/(double)n) * (2*M_PI);;
    MyDisplay();
    const int num = i+1;
    stringstream ss;
    ss << "shot-" << setfill('0') << setw(n_size) << num;
    string fn = ss.str() + ".png";
    /* function undefined. */
    //writePngImage(fn.c_str(), window_width, window_height);
  }
}

// Writes images to shot-xxxx.png
void ExplodeRiceDwarf() {
  scene.SetShowStatistics(false);
  const int n = 150;
  const int n_size = 4;
  // explode dataset
  double end = 0.35;
  // rice dwarf virus
  // double end = 500;
  double offset = ((atan((-n/2.0)/16))+M_PI/2)/M_PI;
  cout << offset << endl;
  for (int i = 0; i < n; ++i) {
    // const float e = ((atan((i-n/2.0)/16)+M_PI/2) / M_PI) * end - 33;
    const float e = ((atan((i-n/2.0)/16)+M_PI/2) / M_PI - offset) * end;
    scene.SetExplode(e);
    cout << scene.GetExplode() << endl;
    MyDisplay();
    const int num = i+1;
    stringstream ss;
    ss << "shot-" << setfill('0') << setw(n_size) << num;
    string fn = ss.str() + ".png";
    /* function undefined. */
    //writePngImage(fn.c_str(), window_width, window_height);
  }
}

// Use explode-knives found in the scripts directory.  First copy it
// to the release directory.
void ExplodeKnives() {
  scene.SetShowStatistics(false);
  // 8,000,000 pulls the knives completely out
  const int max_dist = 8000000;
  // const int max_dist = 800000;
  // const int max = 20000 * 400;

  const int num_frames = 200;
  // const int num_frames = 20;
  const double f = 0.25; // <-- factor to flatten out the velocity curve
  const double max_atan = atan(16);
  // compute how far we'll travel in native units
  double sum = 0;
  for (int frame = 0; frame < num_frames; ++frame) {
    // t ranges linearly from -16 to 16
    const double t = -16 + 32*frame/(double)(num_frames-1);
    // velocity ranges from 0 to 1 to 0
    const double v = pow(1 - fabs(atan(t))/max_atan, f);
    sum += v;
  }
  const double max_dist_native = sum;

  // Execute the explosion
  for (int frame = 0; frame < num_frames; ++frame) {
    const double t = -16 + 32*frame/(double)(num_frames-1);
    const double v = pow(1 - fabs(atan(t))/max_atan, f);
    const int inc = (int)(max_dist * (v/max_dist_native));

    // redirect output
    std::ostream cnull(0);
    streambuf *old_cout = cout.rdbuf();
    cout.rdbuf(cnull.rdbuf());

    // explode
    for (int cur = 0; cur < inc;) {
      cur += scene.Explode(inc - cur);
    }

    // restore output
    cout.rdbuf (old_cout);

    // One step of explosion done.  Save.
    stringstream ss;
    ss << "slice-" << setfill('0') << setw(5) << frame;
    const string base = ss.str();
    {
      // objects
      scene.SetShowMesh(true);
      scene.SetShowGVD(false);
      MyDisplay();
      const string fn = "slice/o/" + base + ".png";
      /* function undefined. */
      //writePngImage(fn.c_str(), window_width, window_height);
    } {
      // objects and GVD
      scene.SetShowMesh(true);
      scene.SetShowGVD(true);
      scene.SetGVDColor(false);
      MyDisplay();
      const string fn = "slice/om/" + base + ".png";
      /* function undefined. */
      //writePngImage(fn.c_str(), window_width, window_height);
    } {
      // objects and inverted GVD
      scene.SetShowMesh(true);
      scene.SetShowGVD(true);
      scene.SetGVDColor(false);
      scene.InvertGVD();
      MyDisplay();
      scene.InvertGVD();
      const string fn = "slice/omi/" + base + ".png";
      /* function undefined. */
      //writePngImage(fn.c_str(), window_width, window_height);
    } {
      // objects and colored GVD
      scene.SetShowMesh(true);
      scene.SetShowGVD(true);
      scene.SetGVDColor(true);
      MyDisplay();
      const string fn = "slice/omc/" + base + ".png";
      /* function undefined. */
      //writePngImage(fn.c_str(), window_width, window_height);
    } {
      // GVD
      scene.SetShowMesh(false);
      scene.SetShowGVD(true);
      scene.SetGVDColor(false);
      MyDisplay();
      const string fn = "slice/m/" + base + ".png";
      /* function undefined. */
      //writePngImage(fn.c_str(), window_width, window_height);
    } {
      // inverted GVD
      scene.SetShowMesh(false);
      scene.SetShowGVD(true);
      scene.SetGVDColor(false);
      scene.InvertGVD();
      MyDisplay();
      scene.InvertGVD();
      const string fn = "slice/mi/" + base + ".png";
      /* function undefined. */
      //writePngImage(fn.c_str(), window_width, window_height);
    } {
      // colored GVD
      scene.SetShowMesh(false);
      scene.SetShowGVD(true);
      scene.SetGVDColor(true);
      MyDisplay();
      const string fn = "slice/mc/" + base + ".png";
      /* function undefined. */
      //writePngImage(fn.c_str(), window_width, window_height);
    } {
      // Write object meshes
      const std::vector<Mesh>& meshes = scene.Meshes();
      for (int label = 1; label < meshes.size(); ++label) {
        stringstream ss2;
        ss2 << "slice/object-meshes/" << base << "-" << label << ".obj";
        const string fn = ss2.str();
        ofstream out(fn.c_str());
        WriteObj(out, meshes[label]);
      }
    } {
      // Write gvd meshes
      const std::vector<Mesh>& meshes = scene.GVDMeshes();
      for (int label = 1; label < meshes.size(); ++label) {
        stringstream ss2;
        ss2 << "slice/gvd-meshes/" << base << "-" << label << ".obj";
        const string fn = ss2.str();
        ofstream out(fn.c_str());
        WriteObj(out, meshes[label]);
      }
    }
    cout << "**** frame = " << frame
          << " inc = " << inc
          << " v = " << v
          << endl;
  }
}

int main(int argc, char *argv[]) {
  using namespace oct;

  // Initialize GLUT
  glutInit(&argc, argv);
  glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA | GLUT_DEPTH);
  glutInitWindowSize(window_width, window_height);

  glutInitWindowPosition(100, 100);
  glutCreateWindow("Object viewer");
  glutMouseFunc(Mouse);
  glutMotionFunc(MouseMotion);
  glutKeyboardFunc(Keyboard);
  glutSpecialFunc(Special);
  glutDisplayFunc(MyDisplay);
  Init();

  const int ret = scene.ProcessArgs(argc, argv);
  if (ret != 0) return ret;

  glutMainLoop();

  return 0;
}
