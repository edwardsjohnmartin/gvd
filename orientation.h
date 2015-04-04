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

#ifndef __ORIENTATION_H__
#define __ORIENTATION_H__

#include <vector>

#include "./opencl/triangle.h"
#include "./opencl/vec.h"

void orient_lean(std::vector<Triangle>& triangles);
std::vector<Triangle> orient_new(const std::vector<Triangle>& triangles);
std::vector<Triangle> orient(const std::vector<Triangle>& triangles);
std::vector<Triangle> split(const std::vector<Triangle>& triangles,
                            std::vector<float3>& vertices);

#endif
