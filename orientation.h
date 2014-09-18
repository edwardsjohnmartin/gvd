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
