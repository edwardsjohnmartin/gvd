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
