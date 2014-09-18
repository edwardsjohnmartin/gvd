#ifndef __TEXTURE_H__
#define __TEXTURE_H__

#ifdef __MAC__

#include <iostream>

void LoadTexture(const std::string& filename, int texture_id);

#else

#include <cstring>
#include <cstdio>
#include <string>

#include <jpeglib.h>

#include "./common.h"
// Loads a texture into OpenGL.  Currently supports only JPEG.
void LoadTexture(const std::string& filename, int texture_id);


GLuint InitTexture(GLuint texture, unsigned char* data,
                   int width, int height, int channels, bool wrap);

bool LoadJPEG(const char* FileName, unsigned char** data,
              jpeg_decompress_struct* out_info,
              bool Fast = true);

GLuint png_texture_load(const char * file_name, int * width, int * height,
                        const GLuint texture);

#endif
#endif
