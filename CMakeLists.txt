cmake_minimum_required(VERSION 2.8)

project(mesh-deform)

if (WIN32)
		add_definitions("-D_SCL_SECURE_NO_WARNINGS -D_CRT_SECURE_NO_WARNINGS -D_ITERATOR_DEBUG_LEVEL=0")
else()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
    set(CMAKE_XCODE_ATTRIBUTE_CLANG_CXX_LIBRARY "libc++")
endif()

set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
set(DEFORM_LINK_TARGETS)

find_package(Eigen3 REQUIRED)
find_package(OpenMesh QUIET)
find_package(OpenSceneGraph QUIET COMPONENTS osgViewer osgFX osgGA)
find_package(Sophus QUIET)

include_directories("inc" ${OPENMESH_INCLUDE_DIR} ${EIGEN3_INCLUDE_DIR} ${SOPHUS_INCLUDE_DIR})
list(APPEND DEFORM_LINK_TARGETS ${OPENMESH_LIBRARIES})

add_library(deform
    inc/deform/arap.h
    inc/deform/openmesh_adapter.h
    inc/deform/trajectory.h
    inc/deform/deformation_util.h
    src/unused.cpp
)
set_target_properties(deform PROPERTIES LINKER_LANGUAGE CXX)

add_executable(deform_tests
    tests/catch.hpp
    tests/accessor.h
    tests/test_cotan.cpp
    tests/test_trajectory.cpp
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

    add_executable(deform_bar
        examples/deform_bar.cpp
        ${OSG_COMMON_FILES}
    )
    target_link_libraries(deform_bar deform ${DEFORM_LINK_TARGETS} ${OPENSCENEGRAPH_LIBRARIES})

    if (SOPHUS_FOUND)
        add_executable(deform_trajectory
            examples/deform_trajectory.cpp
            ${OSG_COMMON_FILES}
        )
        target_link_libraries(deform_trajectory deform ${DEFORM_LINK_TARGETS} ${OPENSCENEGRAPH_LIBRARIES})
    endif ()

endif()
