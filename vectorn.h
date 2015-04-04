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

#ifdef OCT2D
#include "./vector2.h"
namespace oct {
typedef LabeledGeometry2 LabeledGeometry;
}

#else
#include "./vector3.h"
namespace oct {
typedef LabeledGeometry3 LabeledGeometry;
}
#endif

