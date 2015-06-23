cmake_minimum_required (VERSION 2.6)
project (GVD)

SET(SRCS
  ../OctreeUtils.cpp
  ../opencl.cpp
  ../octree.cpp
  ../octree_gpu.cpp
  ../gpu.cpp
  ../ambiguous.cpp
  ../vertex.cpp
  ../vertex_network.cpp
  ../mvertex_network.cpp
  ../vertices_gpu_state.cpp
  ../options.cpp
  ../geometry_cpp.cpp
  ../opencl/bit.cpp
#  ../opencl/distance3.cpp
  ../opencl/geometry.cpp
  ../opencl/uvertex_network.cpp
  ../opencl/cl_ambiguous.cpp
  ../opencl/cl_octree.cpp
)

SET(GVD2_SRCS
  ../vector2.cpp
  ./gl2d.cpp
  ./gvdviewer2.cpp
  ../jrs-triangle.c
)

SET(OCT2_SRCS
  ../vector2.cpp
  ./gl2d.cpp
  ./octviewer2.cpp
  ../jrs-triangle.c
  ../karras.cpp
)

SET(GVD3_SRCS
  ../vector3.cpp
  ../tile3.cpp
  ./gl3d.cpp
  ./gvdviewer3.cpp
  ./texture.cpp
  ./mesh.cpp
  ./material.cpp
  ./io.cpp
  ../opencl/tribox.cpp
  ../opencl/distance3.cpp
 ../orientation.cpp
)

SET(SEGMENT2_SRCS
  ../vector2.cpp
  ./gl2d.cpp
  ./segment2.cpp
  ../jrs-triangle.c
)

ADD_DEFINITIONS(-DANSI_DECLARATORS -DTRILIBRARY)
ADD_DEFINITIONS(-U__STRICT_ANSI__)

IF(APPLE)
  ADD_DEFINITIONS(-D__MAC__)
  #SET(CMAKE_CXX_FLAGS -Wno-deprecated-declarations)
  #list(APPEND CMAKE_CXX_FLAGS -Wno-deprecated-declarations)
ENDIF(APPLE)

FIND_PACKAGE(OpenGL)
FIND_PACKAGE(GLUT)
FIND_PACKAGE(JPEG)
FIND_PACKAGE(PNG)

if (JPEG_FOUND)
include_directories(${JPEG_INCLUDE_DIR})
endif (JPEG_FOUND)

if (PNG_FOUND)
include_directories(${PNG_INCLUDE_DIR})
endif (PNG_FOUND)

option(BUILD_GVD2 "Build 2D GVD" ON)
if(BUILD_GVD2)
  ADD_EXECUTABLE(gvd-viewer2 ${SRCS} ${GVD2_SRCS} viewer2.cpp)
  set_target_properties (gvd-viewer2 PROPERTIES COMPILE_DEFINITIONS "OCT2D")
  TARGET_LINK_LIBRARIES(gvd-viewer2 ${GLUT_LIBRARIES} ${OPENGL_LIBRARIES} ${JPEG_LIBRARIES} ${PNG_LIBRARIES} ${OPENCL_LIBRARY})
endif(BUILD_GVD2)


option(BUILD_GVD3 "Build 3D GVD" ON)
if(BUILD_GVD3)
  ADD_EXECUTABLE(gvd-viewer3 ${SRCS} ${GVD3_SRCS} viewer3.cpp)
  set_target_properties (gvd-viewer3 PROPERTIES COMPILE_DEFINITIONS "OCT3D")
  TARGET_LINK_LIBRARIES(gvd-viewer3 ${GLUT_LIBRARIES} ${OPENGL_LIBRARIES} ${JPEG_LIBRARIES} ${PNG_LIBRARIES} ${OPENCL_LIBRARY})
endif(BUILD_GVD3)

option(BUILD_OCT2 "Build new parallel octree" ON)
if(BUILD_OCT2)
  ADD_EXECUTABLE(oct-viewer2 ${SRCS} ${OCT2_SRCS} oviewer2.cpp)
  set_target_properties (oct-viewer2 PROPERTIES COMPILE_DEFINITIONS "OCT2D")
  TARGET_LINK_LIBRARIES(oct-viewer2 ${GLUT_LIBRARIES} ${OPENGL_LIBRARIES} ${JPEG_LIBRARIES} ${PNG_LIBRARIES} ${OPENCL_LIBRARY})
endif(BUILD_OCT2)

#ADD_EXECUTABLE(segment2 ${SRCS} ${SEGMENT2_SRCS})
#set_target_properties (segment2 PROPERTIES COMPILE_DEFINITIONS "OCT2D")
#TARGET_LINK_LIBRARIES(segment2 ${GLUT_LIBRARIES} ${OPENGL_LIBRARIES} ${JPEG_LIBRARIES} ${PNG_LIBRARIES} ${OPENCL_LIBRARY})

## Copy opencl files
#file(GLOB open_cl_files
#  ${CMAKE_SOURCE_DIR}/opencl/*.cpp
#  ${CMAKE_SOURCE_DIR}/opencl/*.h
#  ${CMAKE_SOURCE_DIR}/opencl/*.cl
#)
#
#message(STATUS "$<TARGET_FILE_DIR:viewer2>")
#
#add_custom_target(copy)
#foreach(open_cl_file ${open_cl_files})
#  add_custom_command(TARGET viewer2 PRE_BUILD
#                     COMMAND ${CMAKE_COMMAND} -E
#                         copy ${open_cl_file} $<TARGET_FILE_DIR:viewer2>)
##                         copy ${open_cl_file} ${viewer2dir}/opencl)
#endforeach()
##add_dependencies(viewer2 copy)
