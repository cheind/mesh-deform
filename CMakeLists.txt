cmake_minimum_required(VERSION 2.8)

project(mesh-deform)

if (WIN32)
		add_definitions("-D_SCL_SECURE_NO_WARNINGS -D_CRT_SECURE_NO_WARNINGS")
else()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
    set(CMAKE_XCODE_ATTRIBUTE_CLANG_CXX_LIBRARY "libc++")
endif()

set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
set(DEFORM_LINK_TARGETS)

set(DEFORM_EIGEN_DIR "../eigen" CACHE PATH "Where is the include directory of Eigen located")
find_package(OpenMesh REQUIRED)

include_directories("inc" ${OPENMESH_INCLUDE_DIR} ${DEFORM_EIGEN_DIR})
list(APPEND DEFORM_LINK_TARGETS ${OPENMESH_LIBRARIES})

add_library(deform
    inc/deform/mesh.h
    inc/deform/arap.h
    src/arap.cpp
)

add_executable(deform_tests
    tests/catch.hpp
    tests/test_weights.cpp
)
target_link_libraries(deform_tests deform ${DEFORM_LINK_TARGETS})

add_executable(deform_example
    examples/deform_example.cpp
)
target_link_libraries(deform_example deform ${DEFORM_LINK_TARGETS})
