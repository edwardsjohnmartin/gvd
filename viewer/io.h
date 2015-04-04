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

#ifndef __IO_H__
#define __IO_H__

#include <string>
#include <iostream>
#include <fstream>

#include "./mesh.h"

bool ParseObj(const std::string& fn, Mesh& mesh);
void WriteObj(std::ostream& out, const Mesh& mesh);
void WritePly(std::ostream& out, const Mesh& mesh);

// Writes a density function file usable by gx_ccvt.
// The density at each point is computed using adjacent
// triangle areas.  Currently supports only triangles.
void WriteDensity(std::ostream& out, const Mesh& mesh);

#endif
