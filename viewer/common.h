#ifndef __COMMON_GL_INCLUDES_H__
#define __COMMON_GL_INCLUDES_H__

#ifdef _WIN32
#include <windows.h>
#endif

// #define GL_GLEXT_PROTOTYPES
#ifdef __MAC__
  #include <OpenGL/gl.h>
  #include <GLUT/glut.h>
  #include <OpenGL/glext.h>
#else
  #define GL_GLEXT_PROTOTYPES
#   define ANT_UNIX
#   include <X11/cursorfont.h>
#   define GLX_GLXEXT_LEGACY
#   include <GL/glx.h>
#   include <X11/Xatom.h>
#   include <unistd.h>
#   include <malloc.h>
#   undef _WIN32
#   undef WIN32
#   undef _WIN64
#   undef WIN64
#   undef _WINDOWS
#   undef ANT_WINDOWS
#   undef ANT_OSX
#	    include <GL/gl.h>  // must be included after windows.h
  #define GL_GLEXT_PROTOTYPES
  #include <GL/gl.h>
  #include <GL/glu.h>
  #include <GL/glut.h>
  #include <GL/glext.h>
  // #define GLX_GLXEXT_LEGACY
  // #include <GL/glx.h>
  // #define GL_GLEXT_PROTOTYPES
  // #include <GL/gl.h>
  // #include <GL/glu.h>
  // #include <GL/glut.h>
  // #include <GL/glext.h>
  // // // #include <GL/glew.h>
  // // #include <GL/gl.h>
  // // #include <GL/glut.h>
  // // #include <GL/glext.h>
#endif

#include <algorithm>
#include <iostream>
#include <iomanip>

// void Init();
// void Display();
// void Keyboard(unsigned char key, int x, int y);
// void DisableLighting();
// void EnableLighting();
// void DrawAxis();
// void DrawAxes();
void PrintMatrix(GLint matrix);
void PrintMatrix();
void LoadMatrix(GLfloat* m);
void MultMatrix(GLfloat* m);

inline void PrintMatrix(GLfloat* m) {
  using namespace std;
  cout.precision(2);
  int w = 6;
  for (int i = 0; i < 4; ++i) {
    cout << setprecision(2) << setw(w) << m[i] << " "
        << setprecision(2) << setw(w) << m[i+4] << " "
        << setprecision(2) << setw(w) << m[i+8] << " "
        << setprecision(2) << setw(w) << m[i+12] << " "
        << endl;
  }
  cout << endl;
}

inline void PrintMatrix(GLint matrix) {
  GLfloat m[16];
  glGetFloatv(matrix, m);
  PrintMatrix(m);
}

inline void PrintMatrix() {
  PrintMatrix(GL_MODELVIEW_MATRIX);
}

inline void LoadMatrix(GLfloat* m) {
  // transpose to column-major
  for (int i = 0; i < 4; ++i) {
    for (int j = i; j < 4; ++j) {
      std::swap(m[i*4+j], m[j*4+i]);
    }
  }
  glLoadMatrixf(m);
}

inline void MultMatrix(GLfloat* m) {
  // transpose to column-major
  for (int i = 0; i < 4; ++i) {
    for (int j = i; j < 4; ++j) {
      std::swap(m[i*4+j], m[j*4+i]);
    }
  }
  glMultMatrixf(m);
}

bool IsExtensionSupported(const char* szTargetExtension);

#endif
