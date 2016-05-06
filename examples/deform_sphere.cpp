/**
 This file is part of mesh-deform.
 
 Copyright(C) 2016 Christoph Heindl
 All rights reserved.
 
 This software may be modified and distributed under the terms
 of the BSD license.See the LICENSE file for details.

 This examples shows the deformation of a sphere like object. In this particular
 case a single vertex is pinned to its initial location and a harmonic motion
 is applied to the vertex on the opposite side of the sphere. All other vertices
 are unconstrained.
 */

#define _USE_MATH_DEFINES

// Include as-rigid-as-possible deformation algorithm
#include <deform/arap.h>

// Include an adapter that provides mesh conversion from OpenMesh types.
#include <deform/openmesh_adapter.h>

// For OpenMesh default triangular meshes an OSG 3D viewer is provided.
#include "osg_viewer.h"
#include "example_config.h"

#include <string>
#include <iostream>

int main(int argc, char **argv) {

    // Define the mesh type. Here we use OpenMesh defaults.
    typedef OpenMesh::TriMesh_ArrayKernelT<> Mesh;

    
    // Load the triangulated sphere.
    std::string pathToSource = std::string(DEFORM_ETC_DIR) + std::string("/sphere.obj");
    
    Mesh mesh;
    if (!OpenMesh::IO::read_mesh(mesh, pathToSource)) {
        std::cerr << "Failed to read mesh" << std::endl;
        return -1;
    }
    
  
    // Mesh-deform is independent of the types used to represent triangular surfaces.
    // We need to provide an adapter so that Mesh-deform is able to handle OpenMesh.
    // For OpenMesh a ready-to-use adapter is provided.
    
    typedef deform::OpenMeshAdapter<> Adapter;
    Adapter ma(mesh);
    
    // Instance the deformation algorithm.
    deform::AsRigidAsPossibleDeformation<Adapter> arap(ma);
    
    // Set anchors. Semantically anchors are vertices that are pinned to a specific location.
    // Here we just pin one vertex to its initial location.
    
    Mesh::VertexHandle va = mesh.vertex_handle(37);
    arap.setConstraint(va.idx(), deform::convert::toEigen(mesh.point(va)));
    
    // Motion related variables
    
    Mesh::Point p = mesh.point(mesh.vertex_handle(32));
    double pi = -M_PI_2;
    double add = 0.01;
    
    // Create an OSG viewer to visualize incremental deformation.
    deform::example::OSGViewer viewer(argc, argv, mesh, [&p, &pi, &add, &arap](Mesh &mesh, double time) {
        
        // Update motion
        
        pi += add;
        if (std::abs(pi) > M_PI_2) {
            add *= -1.0;
            pi += add;
        }
        
        // Set handles. Semantically handles are vertices that are dragged around by some motion. Note,
        // although semantically they have a different meaning conceptually its the same as pinning a vertex
        // to a location.
        
        Mesh::VertexHandle vh = mesh.vertex_handle(32);
        arap.setConstraint(vh.idx(), deform::convert::toEigen(p + Mesh::Point(0, 0, (float)(std::sin(pi)))));
        
        // Run the algorithm for 5 iterations. Usually enough when the target constrain locations are not
        // far from the initial locations. Mesh is automatically updated when method completes.
        
        arap.deform(5);
        
        return true;
    });
    viewer.run();
    
    return 0;

}
