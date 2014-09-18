#ifndef __PNGIO_H__
#define __PNGIO_H__

// #ifdef __MAC__

// #include <iostream>

// static void writePngImage(const char *filepath, int width, int height) {
//   std::cout << "writePngImage not currently implemented on Mac OS X"
//             << std::endl;
// }

// #else

#include <png.h>
#include <iostream>
#include <stdlib.h>
#include <stdint.h>
#include "./common.h"

static void readPixelValues(GLubyte* pixels, const int width, const int height) {
  glReadPixels(0, 0, width, height, GL_RGB, GL_UNSIGNED_BYTE, pixels);
  return;
}

static void writePngImage(const char *filepath, int width, int height)
{

  
  FILE * fp;
  png_structp png_ptr = NULL;
  png_infop info_ptr = NULL;
  // size_t x;
  size_t y;
  png_byte ** row_pointers = NULL;
  
  int pixel_size = 3;
  int depth = 8;

  GLubyte *pixels = (GLubyte*) malloc(pixel_size* width * height);
  glutSwapBuffers();
  readPixelValues(pixels, width, height);
  glutSwapBuffers();
  
  fp = fopen (filepath, "wb");
  if (! fp) {
    cout << "Failed to Open file for writing" << endl;
    return;
  }

  png_ptr = png_create_write_struct (PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
  if (png_ptr == NULL) {
    cout << "Png_create_write_struct failed" << endl;
    return;
  }

  info_ptr = png_create_info_struct (png_ptr);
  if (info_ptr == NULL) {
    cout << "Png_create_info_struct failed" << endl;
    return;
  }

  if (setjmp (png_jmpbuf (png_ptr))) {
    cout << "Error during png init_io" << endl;
    return;
  }

  png_set_IHDR (png_ptr,
                info_ptr,
                width,
                height,
                depth,
                PNG_COLOR_TYPE_RGB,
                PNG_INTERLACE_NONE,
                PNG_COMPRESSION_TYPE_DEFAULT,
                PNG_FILTER_TYPE_DEFAULT);


  row_pointers = (png_byte**) png_malloc (png_ptr, height * sizeof (png_byte *));
  for (int y = 0; y < height; ++y) {
    png_byte *row = 
        (png_byte*) png_malloc (png_ptr, sizeof (GLubyte) * width * pixel_size);
    row_pointers[height-y-1] = row;
    for (int x = 0; x < width; ++x) {
      *row++ = pixels[pixel_size*(y*width + x)];
      *row++ = pixels[pixel_size*(y*width + x)+1];
      *row++ = pixels[pixel_size*(y*width + x)+2];
    }
  }

  png_init_io (png_ptr, fp);
  png_set_rows (png_ptr, info_ptr, row_pointers);
  png_write_png (png_ptr, info_ptr, PNG_TRANSFORM_IDENTITY, NULL);

  for (y = 0; y < height; y++) {
    png_free (png_ptr, row_pointers[y]);
  }
  png_free (png_ptr, row_pointers);

  fclose (fp);
  // cout << "Screenshot written" << endl;
  
  return;
}

// #endif

#endif
