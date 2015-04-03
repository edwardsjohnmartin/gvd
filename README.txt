gvd-viewer2 / gvd-viewer3
A 2- and 3-Dimensional Generalized Voronoi Diagram approximator
Version 1.0

Copyright 2015 by John Martin Edwards
Please send bugs and comments to edwardsjohnmartin@gmail.com

There is no warranty.

gvd-viewerx constructs a quadtree/octree that resolves between objects,
propagates a distance transform over the octree vertices, then builds a GVD
approximation using the labels and distances on the vertices. For more
information, please see:

http://sci.utah.edu/~jedwards/research/gvd/index.html

------------------------------------------------------------------------------

If you use gvd-viewerx for a publication, please cite

John Edwards, "Approximating the Generalized Voronoi Diagram of closely spaced objects," in Computer Graphics Forum (Eurographics 2015), 2015.

------------------------------------------------------------------------------

These programs may be freely redistributed under the condition that the
copyright notices (including the copy of this notice in the code comments
and the copyright notice printed when the `-h' switch is selected) are
not removed, and no compensation is received.  Private, research, and
institutional use is free.  You may distribute modified versions of this
code UNDER THE CONDITION THAT THIS CODE AND ANY MODIFICATIONS MADE TO IT
IN THE SAME FILE REMAIN UNDER COPYRIGHT OF THE ORIGINAL AUTHOR, BOTH
SOURCE AND OBJECT CODE ARE MADE FREELY AVAILABLE WITHOUT CHARGE, AND
CLEAR NOTICE IS GIVEN OF THE MODIFICATIONS.  Distribution of this code as
part of a commercial system is permissible ONLY BY DIRECT ARRANGEMENT
WITH THE AUTHOR.  (If you are not directly supplying this code to a
customer, and you are instead telling them how they can obtain it for
free, then you are not required to make any arrangement with me.)

------------------------------------------------------------------------------

BUILD INSTRUCTIONS

  > mkdir build
  > cd build
  > cmake ..
  > make

EXECUTION EXAMPLES

To run the 2D GVD viewer in interactive mode, use
  > ./gvd-viewer2
In the blank window click and drag to create polygons.

To run the 2D GVD viewer on pre-determined polygons, use
  > ./gvd-viewer2 -l 8 ../data2/test[5-8].dat
This creates the GVD on a sample dataset. Each .dat file is a list of vertex
coordinates of the polygon.

To run the 3D GVD viewer on pre-determined polyhedra, use
  > ./gvd-viewer3 -l 8 ../data3/simple*.obj
This creates the GVD on a sample dataset. Currently the GVD viewer supports
polyhedra input only in .obj format.

To write the computed GVD, press 'w'.

OPENCL

gvd-viewerx has an OpenCL acceleration option. Before building, use ccmake
or a CMake GUI interface to set the build variable OPENCL_ACCEL to ON.
