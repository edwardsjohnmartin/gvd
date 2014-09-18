#include "./texture.h"

#ifdef __MAC__

void LoadTexture(const std::string& filename, int texture_id) {
  std::cout << "LoadTexture not currently implemented on Mac OS X"
            << std::endl;
}

#else

#include <jpeglib.h>
#include <jerror.h>
#include <png.h>

#include <string>
#include <iostream>
#include <cstring>
#include <cstdio>

#include "./common.h"

using namespace std;

namespace {
typedef unsigned char BYTE;
typedef GLuint UINT;
}  // end empty namespace
void LoadTexture(const string& filename, int texture_id) {
  
  jpeg_decompress_struct info;
  unsigned char* data;
  const string fn = filename;
  if (!fn.empty() && LoadJPEG(fn.c_str(), &data, &info, true)) {
    InitTexture(texture_id, data, info.output_width,
                info.output_height, info.num_components, true);
  } else {
    cout << "Failed to load " << filename << endl;
  }

}


// Load a texture into OpenGL
GLuint InitTexture(GLuint texture, unsigned char* data,
                   int width, int height, int channels, bool wrap) {
  typedef unsigned char BYTE;

  // select our current texture
  glBindTexture(GL_TEXTURE_2D, texture);

  if (channels == 3) {
    glTexImage2D(
        GL_TEXTURE_2D, 0, GL_RGB8, width, height, 0, GL_RGB,
        GL_UNSIGNED_BYTE, data);
    // glGenerateMipmap(GL_TEXTURE_2D);  //Generate mipmaps now!!!
    gluBuild2DMipmaps(
        GL_TEXTURE_2D, GL_RGB, width, height,
        GL_RGB, GL_UNSIGNED_BYTE, data);
  } else if (channels == 4) {
    glTexImage2D(
        GL_TEXTURE_2D, 0, GL_RGBA8, width, height, 0, GL_RGBA,
        GL_UNSIGNED_BYTE, data);
    // glGenerateMipmap(GL_TEXTURE_2D);  //Generate mipmaps now!!!
    gluBuild2DMipmaps(
        GL_TEXTURE_2D, GL_RGBA8, width, height,
        GL_RGBA, GL_UNSIGNED_BYTE, data);
  }

  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
  // glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
  // glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
  // glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER,
  //                 GL_LINEAR_MIPMAP_NEAREST);
  // glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER,
  //                 GL_LINEAR_MIPMAP_NEAREST);

  // Mix texture with color for shading
  glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);

  // If wrap is true, the texture wraps over at the edges (repeat)
  //       ... false, the texture ends at the edges (clamp)
  if (wrap) {
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
  } else {
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP);
  }

  free(data);

  return texture;
}

bool LoadJPEG(const char* FileName, BYTE** out_data,
              jpeg_decompress_struct* out_info,
              bool Fast) {
  FILE* file = fopen(FileName, "rb");
  jpeg_decompress_struct& info = *out_info;
  // error handler
  struct jpeg_error_mgr err;

  // Tell the jpeg decompression handler to send the errors to err
  info.err = jpeg_std_error(&err);
  // Set info defaults
  jpeg_create_decompress(&info);

  if (!file) {
    fprintf(stderr, "Error reading JPEG file %s\n", FileName);
    return false;
  }

  // Tell the jpeg lib the file we're reading and read it
  jpeg_stdio_src(&info, file);
  jpeg_read_header(&info, TRUE);

  if (Fast) {
    info.do_fancy_upsampling = FALSE;
  }

  jpeg_start_decompress(&info);

  const size_t x = info.output_width;
  const size_t y = info.output_height;
  const int channels = info.num_components;
  const size_t size = x * y * channels;

  // Turn the uncompressed data into something ogl can read
  BYTE* data = new BYTE[size];

  BYTE* p1 = data;
  BYTE** p2 = &p1;
  int numlines = 0;

  while (info.output_scanline < info.output_height) {
    numlines = jpeg_read_scanlines(&info, p2, 1);
    *p2 += numlines * channels * info.output_width;
  }

  jpeg_finish_decompress(&info);
  fclose(file);

  // If scan lines are not a multiple of 4, then change channels to 4.
  // We only support converting from 3 right now.
  if ((x * channels) % 4 != 0 && channels == 3) {
    BYTE* new_data = new BYTE[x * y * 4];
    memset(new_data, x*y*4, 1);
    for (int i = 0; i < y; ++i) {
      for (int j = 0; j < x; ++j) {
        memcpy(new_data+(i*x+j)*4, data+(i*x+j)*3, 3);
      }
    }
    delete data;
    data = new_data;

    info.num_components = 4;
  }

  *out_data = data;

  return true;
}

GLuint png_texture_load(const char * file_name, int * width, int * height,
                        const GLuint texture) {
  
  png_byte header[8];

  FILE *fp = fopen(file_name, "rb");
  if (fp == 0) {
    perror(file_name);
    return 0;
  }

  // read the header
  fread(header, 1, 8, fp);

  if (png_sig_cmp(header, 0, 8)) {
    fprintf(stderr, "error: %s is not a PNG.\n", file_name);
    fclose(fp);
    return 0;
  }

  png_structp png_ptr = png_create_read_struct(
      PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
  if (!png_ptr) {
    fprintf(stderr, "error: png_create_read_struct returned 0.\n");
    fclose(fp);
    return 0;
  }

  // create png info struct
  png_infop info_ptr = png_create_info_struct(png_ptr);
  if (!info_ptr) {
    fprintf(stderr, "error: png_create_info_struct returned 0.\n");
    png_destroy_read_struct(&png_ptr, (png_infopp)NULL, (png_infopp)NULL);
    fclose(fp);
    return 0;
  }

  // create png info struct
  png_infop end_info = png_create_info_struct(png_ptr);
  if (!end_info) {
    fprintf(stderr, "error: png_create_info_struct returned 0.\n");
    png_destroy_read_struct(&png_ptr, &info_ptr, (png_infopp) NULL);
    fclose(fp);
    return 0;
  }

  // the code in this if statement gets called if libpng encounters an error
  if (setjmp(png_jmpbuf(png_ptr))) {
    fprintf(stderr, "error from libpng\n");
    png_destroy_read_struct(&png_ptr, &info_ptr, &end_info);
    fclose(fp);
    return 0;
  }

  // init png reading
  png_init_io(png_ptr, fp);

  // let libpng know you already read the first 8 bytes
  png_set_sig_bytes(png_ptr, 8);

  // read all the info up to the image data
  png_read_info(png_ptr, info_ptr);

  // variables to pass to get info
  int bit_depth, color_type;
  png_uint_32 temp_width, temp_height;

  // get info about png
  png_get_IHDR(png_ptr, info_ptr, &temp_width,
               &temp_height, &bit_depth, &color_type,
               NULL, NULL, NULL);

  if (width) *width = temp_width;
  if (height) *height = temp_height;

  // Update the png info struct.
  png_read_update_info(png_ptr, info_ptr);

  // Row size in bytes.
  int rowbytes = png_get_rowbytes(png_ptr, info_ptr);

  // glTexImage2d requires rows to be 4-byte aligned
  rowbytes += 3 - ((rowbytes-1) % 4);

  // Allocate the image_data as a big block, to be given to opengl
  // png_byte * image_data;
  // image_data =
  //    (png_byte*) malloc(rowbytes * temp_height * sizeof(png_byte)+15);
  const size_t img_num_bytes = rowbytes * temp_height * sizeof(png_byte)+15;
  png_byte* image_data = static_cast<png_byte*>(malloc(img_num_bytes));
  if (image_data == NULL) {
    fprintf(stderr, "error: could not allocate memory for PNG image data\n");
    png_destroy_read_struct(&png_ptr, &info_ptr, &end_info);
    fclose(fp);
    return 0;
  }

  // row_pointers is for pointing to image_data for reading the png with libpng
  // png_bytep * row_pointers = (png_bytep*)malloc(
  //      temp_height * sizeof(png_bytep));
  const int row_num_bytes = temp_height * sizeof(png_bytep);
  png_bytep * row_pointers = static_cast<png_bytep*>(malloc(row_num_bytes));
  if (row_pointers == NULL) {
    fprintf(stderr, "error: could not allocate memory for PNG row pointers\n");
    png_destroy_read_struct(&png_ptr, &info_ptr, &end_info);
    free(image_data);
    fclose(fp);
    return 0;
  }

  // set the individual row_pointers to point at the correct offsets
  // of image_data
  int i;
  for (i = 0; i < temp_height; i++) {
    row_pointers[temp_height - 1 - i] = image_data + i * rowbytes;
  }

  // read the png into image_data through row_pointers
  png_read_image(png_ptr, row_pointers);

  // Generate the OpenGL texture object
  // GLuint texture;
  // glGenTextures(1, &texture);
  glBindTexture(GL_TEXTURE_2D, texture);

  GLenum type = GL_UNSIGNED_BYTE;
  //if (info_ptr->bit_depth == 16) {
  if (png_get_bit_depth(png_ptr, info_ptr) == 16) {
    type = GL_UNSIGNED_SHORT;
  }
  GLint internal_format = GL_RGB8;
  GLenum format = GL_RGB;
  //if (info_ptr->color_type == PNG_COLOR_TYPE_RGBA) {
  if (png_get_color_type(png_ptr, info_ptr) == PNG_COLOR_TYPE_RGBA) {
    format = GL_RGBA;
  }
  glTexImage2D(GL_TEXTURE_2D, 0, internal_format, temp_width, temp_height,
               0, format, type, image_data);

  // gluBuild2DMipmaps(GL_TEXTURE_2D, GL_RGB, temp_width, temp_height,
  //                   GL_RGB, GL_UNSIGNED_BYTE, image_data);

  // glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MIN_LOD, 0);
  // glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MAX_LOD, 0);
  glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
  glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
  // glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
  // glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
  // glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER,
  //                 GL_LINEAR_MIPMAP_LINEAR);
  // glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER,
  //                 GL_LINEAR_MIPMAP_LINEAR);

  // clean up
  png_destroy_read_struct(&png_ptr, &info_ptr, &end_info);
  free(image_data);
  free(row_pointers);
  fclose(fp);

  return texture;
}

#endif
