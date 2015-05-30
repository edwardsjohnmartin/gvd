cmake_minimum_required (VERSION 2.6)

SET(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

#------------------------------------------------------------
# Use C++ 2011
#------------------------------------------------------------
include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-std=c++11" COMPILER_SUPPORTS_CXX11)
CHECK_CXX_COMPILER_FLAG("-std=c++0x" COMPILER_SUPPORTS_CXX0X)
if(COMPILER_SUPPORTS_CXX11)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
elseif(COMPILER_SUPPORTS_CXX0X)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")
else()
        message(STATUS "The compiler ${CMAKE_CXX_COMPILER} has no C++11 support. Please use a different C++ compiler.")
endif()

#------------------------------------------------------------
# Use OpenCL depending on setting
#------------------------------------------------------------
OPTION(OPENCL_ACCEL "Use OpenCL acceleration" OFF)
if(OPENCL_ACCEL)
  FIND_PACKAGE(OpenCL)
  if (OPENCL_FOUND)
    INCLUDE_DIRECTORIES(${OPENCL_INCLUDE_DIR})
    ADD_DEFINITIONS(-D__OPEN_CL_SUPPORT__)
  else (OPENCL_FOUND)
    set(OPENCL_LIBRARY "")
  endif (OPENCL_FOUND)
endif(OPENCL_ACCEL)

set(EXECUTABLE_OUTPUT_PATH "${CMAKE_BINARY_DIR}")

add_subdirectory (viewer)
