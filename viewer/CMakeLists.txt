cmake_minimum_required (VERSION 2.6)
project (Viewer)

SET(SRCS
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

SET(VIEWER2_SRCS
  ../vector2.cpp
  ./gl2d.cpp
  ./medial2.cpp
  ../jrs-triangle.c
)

SET(VIEWER3_SRCS
  ../vector3.cpp
  ../tile3.cpp
  ./gl3d.cpp
  ./medial3.cpp
  ./texture.cpp
  ./mesh.cpp
  ./material.cpp
  ./io.cpp
  ../opencl/tribox.cpp
  ../opencl/distance3.cpp
 ../orientation.cpp
)

ADD_DEFINITIONS(-DANSI_DECLARATORS -DTRILIBRARY)

IF(APPLE)
  ADD_DEFINITIONS(-D__MAC__)
  SET(CMAKE_CXX_FLAGS -Wno-deprecated-declarations)
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

ADD_EXECUTABLE(viewer2 ${SRCS} ${VIEWER2_SRCS} viewer2.cpp)
set_target_properties (viewer2 PROPERTIES COMPILE_DEFINITIONS "OCT2D")
TARGET_LINK_LIBRARIES(viewer2 ${GLUT_LIBRARIES} ${OPENGL_LIBRARIES} ${JPEG_LIBRARIES} ${PNG_LIBRARIES} ${OPENCL_LIBRARY})

ADD_EXECUTABLE(viewer3 ${SRCS} ${VIEWER3_SRCS} viewer3.cpp)
set_target_properties (viewer3 PROPERTIES COMPILE_DEFINITIONS "OCT3D")
TARGET_LINK_LIBRARIES(viewer3 ${GLUT_LIBRARIES} ${OPENGL_LIBRARIES} ${JPEG_LIBRARIES} ${PNG_LIBRARIES} ${OPENCL_LIBRARY})

#ADD_EXECUTABLE(presentation ${VIEWER2_SRCS} ${VIEWER3_SRCS} presentation.cpp)
#TARGET_LINK_LIBRARIES(presentation ${GLUT_LIBRARIES} ${OPENGL_LIBRARIES} ${JPEG_LIBRARIES} ${PNG_LIBRARIES})

#add_custom_target(cpplint
#    COMMAND ../cpplint.py --filter=-readability/todo,-readability/streams,-runtime/threadsafe_fn,-legal,-build/namespaces,-build/header_guard,-runtime/sizeof,-build/include_order,-runtime/references,-runtime/explicit ${CMAKE_SOURCE_DIR}/*.cpp ${CMAKE_SOURCE_DIR}/*.h ${CMAKE_SOURCE_DIR}/../*.cpp ${CMAKE_SOURCE_DIR}/../*.cpp)
