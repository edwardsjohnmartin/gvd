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

#ifndef __TRI_BOX_H__
#define __TRI_BOX_H__

#include "./defs.h"
#include "./triangle.h"
#include "./vec.h"

// int triBoxOverlap(const float boxcenter[3],
//                   const float boxhalfsize[3],
//                   const float triverts[3][3]);
//                   // const float** triverts);

bool TriBoxOverlap(const int3 boxcenter,
                   const int boxhalfsize,
                   const __GLOBAL__ int3* vertices,
                   const Triangle tri);

#endif
