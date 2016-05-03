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

set(EIGEN_INCLUDE_DIR "../eigen" CACHE PATH "Where is the include directory of Eigen located")
find_package(OpenMesh REQUIRED)
find_package(OpenSceneGraph QUIET COMPONENTS osgViewer osgFX osgGA)

include_directories("inc" ${OPENMESH_INCLUDE_DIR} ${EIGEN_INCLUDE_DIR})
list(APPEND DEFORM_LINK_TARGETS ${OPENMESH_LIBRARIES})

add_library(deform
    inc/deform/mesh.h
    inc/deform/cotan_matrix.h
    inc/deform/arap.h
    src/arap.cpp
    src/cotan_matrix.cpp
)

add_executable(deform_tests
    tests/catch.hpp
    tests/test_cotan.cpp
)
target_link_libraries(deform_tests deform ${DEFORM_LINK_TARGETS})

set(DEFORM_ETC_DIR "${CMAKE_CURRENT_SOURCE_DIR}/etc")
configure_file(examples/example_config.h.in example_config.h)
include_directories(${CMAKE_CURRENT_BINARY_DIR})
include_directories("examples")

if (OPENSCENEGRAPH_FOUND)
		set(OSG_COMMON_FILES examples/osg_viewer.h examples/osg_viewer.cpp)

    include_directories(${OPENSCENEGRAPH_INCLUDE_DIRS})
    add_executable(deform_sphere
        examples/deform_sphere.cpp
				${OSG_COMMON_FILES}
    )
    target_link_libraries(deform_sphere deform ${DEFORM_LINK_TARGETS} ${OPENSCENEGRAPH_LIBRARIES})
endif()
