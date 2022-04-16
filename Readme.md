## mesh-deform
**mesh-deform** is cross-platform C++11 library capable of deforming meshes under rigidity constraints, leading to an interactive visually appealing (not necessarily physically correct) surface modeling approach.

The implementation is based on *Sorkine, Olga, and Marc Alexa. "As-rigid-as-possible surface modeling." Symposium on Geometry processing. Vol. 4. 2007.*


### Features

 - Mesh independent implementation running interactive frame rates
 - Keyframe-based 6DoF trajectory interpolation
 - OpenMesh surface adapter
 - Visualizer

### Compilation

**mesh-deform** is a header only project requiring the following dependencies:
 - Eigen: http://eigen.tuxfamily.org/

Optionally:
 - Sophus: https://github.com/stevenlovegrove/Sophus/
 - OpenMesh: https://www.openmesh.org/
 - OSG: http://www.openscenegraph.org/

### Getting started

You may want to browse the source of one of the examples to get a feeling for how this library is meant to be used.

The basic procedure is: Given a source mesh, a set of handles that are said to be fixed in space, another set of handles dragged to an offset position the library iteratively computes the positions of the remaining vertices so that every vertex is moved as rigid as possible given its direct neighborhood.

### Videos
Click to view.

[![Sphere Example](https://img.youtube.com/vi/h7YqRzAhKVw/maxresdefault.jpg)](https://www.youtube.com/watch?v=h7YqRzAhKVw)
[![Bar Example](https://img.youtube.com/vi/0BI0iYPgbyo/maxresdefault.jpg )](https://www.youtube.com/watch?v=0BI0iYPgbyo)
