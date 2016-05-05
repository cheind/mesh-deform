/**
 This file is part of mesh-deform.
 
 Copyright(C) 2016 Christoph Heindl
 All rights reserved.
 
 This software may be modified and distributed under the terms
 of the BSD license.See the LICENSE file for details.
 */


#include <deform/arap.h>
#include <deform/openmesh_adapter.h>
#include <string>
#include <iostream>

#include "osg_viewer.h"
#include "example_config.h"

int main(int argc, char **argv) {
    
    std::string pathToSource = std::string(DEFORM_ETC_DIR) + std::string("/sphere.obj");

    typedef OpenMesh::TriMesh_ArrayKernelT<> Mesh;
    
    Mesh mesh;
    if (!OpenMesh::IO::read_mesh(mesh, pathToSource)) {
        std::cerr << "Failed to read mesh" << std::endl;
        return -1;
    }
    
  
    Mesh::Point p = mesh.point(mesh.vertex_handle(32));
    double pi = -M_PI_2;
    double add = 0.01;
    
    deform::OpenMeshAdapter<> ma(mesh);
    deform::AsRigidAsPossibleDeformation< deform::OpenMeshAdapter<> > arap(ma);
    
    // Set anchors
    Mesh::VertexHandle va = mesh.vertex_handle(37);
    arap.setConstraint(va.idx(), deform::convert::toEigen(mesh.point(va)));
    
    deform::example::OSGViewer viewer(argc, argv, mesh, [&p, &pi, &add, &arap](Mesh &mesh, double time) {
        
        pi += add;
        if (std::abs(pi) > M_PI_2) {
            add *= -1.0;
            pi += add;
        }
        
        Mesh::VertexHandle vh = mesh.vertex_handle(32);
        arap.setConstraint(vh.idx(), deform::convert::toEigen(p + Mesh::Point(0, 0, (float)(std::sin(pi)))));
        arap.deform(3);
        
        return true;
    });
    viewer.run();
    
    return 0;

}
