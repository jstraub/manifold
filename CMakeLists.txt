# Copyright (c) 2015, Julian Straub <jstraub@csail.mit.edu>
# Licensed under the MIT license. See the license file LICENSE.

project(manifold)
cmake_minimum_required(VERSION 2.8.1)
set(POD_NAME manifold)
include(cmake/pods.cmake)

if(COMMAND cmake_policy)
  cmake_policy(SET CMP0003 NEW)
endif(COMMAND cmake_policy)

set(CMAKE_COLOR_MAKEFILE ON)
IF(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release) # Release; debug; ...
  set(CMAKE_BUILD_TYPE Debug) # Release; debug; ...
ENDIF()

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)
set(CMAKE_EXPORT_COMPILE_COMMANDS 1 )

if (CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  set (CMAKE_INSTALL_PREFIX "${CMAKE_BINARY_DIR}/install" CACHE PATH "default install path" FORCE )
endif()

find_package (Eigen3 REQUIRED)

# add flags depending on build type
set(CMAKE_CXX_FLAGS_RELEASE "-g -Wall -Wno-deprecated -O2 -finline-functions -DNDEBUG -std=c++11 ")
set(CMAKE_CXX_FLAGS_DEBUG "-g -Wno-deprecated -std=c++11")

include_directories(
  include
  ${EIGEN3_INCLUDE_DIR}
  )
link_directories(
)

add_executable(testSO3GD ./test/SO3GD.cpp)
add_executable(testSE3GD ./test/SE3GD.cpp)
add_executable(testSO3 ./test/SO3.cpp)
add_executable(testSE3 ./test/SE3.cpp)
add_executable(testSO3_incSurfNormAlign ./test/SO3_incSurfNormAlign.cpp)
add_executable(testS ./test/S.cpp)
add_executable(testSO3_gmm ./test/SO3_gmm.cpp)

pods_install_executables(testSO3 testS)
#pods_install_libraries(optRot)
pods_install_headers(
  ./include/manifold/S.h
  ./include/manifold/S_impl.hpp
  ./include/manifold/SO3.h
  ./include/manifold/SE3.h
  ./include/manifold/manifold.h
  ./include/manifold/gradientDescent.h
  ./include/manifold/gradientDescentSO3.h
  ./include/manifold/gradientDescentSE3.h
  ./include/manifold/SO3_impl.hpp
  ./include/manifold/SE3_impl.hpp
  DESTINATION manifold)
