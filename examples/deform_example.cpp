/**
 This file is part of mesh-deform.
 
 Copyright(C) 2016 Christoph Heindl
 All rights reserved.
 
 This software may be modified and distributed under the terms
 of the BSD license.See the LICENSE file for details.
 */


#include <deform/arap.h>
#include <deform/cotan_matrix.h>
#include <iostream>

int main(int argc, char **argv) {

    if (argc != 2) {
        std::cerr << argv[0] << " plane.ply anchors.txt handles.txt" << std::endl;
        return -1;
    }

    deform::Mesh mesh;
    if (!OpenMesh::IO::read_mesh(mesh, argv[1])) {
        std::cerr << "Failed to read mesh" << std::endl;
        return -1;
    }
    
    const int nids = 1;

    //const int anchorIds[] = { 0,1,2,3,40,41,60,61,80,81,100,101,120,121,140,141,160,161,180,181,200,201,220,221,240,241,260,261,280,281,300,301,320,321,340,341,360,361,380,381 };
    const int anchorIds[] = { 37 };
    const int handleIds[] = { 36,37,38,39,58,59,78,79,98,99,118,119,138,139,158,159,178,179,198,199,218,219,238,239,258,259,278,279,298,299,318,319,338,339,358,359,378,379,398,399 };

    deform::AsRigidAsPossibleDeformation arap(&mesh);

    for (int i = 0; i < nids; ++i) {
        deform::Mesh::VertexHandle vh = mesh.vertex_handle(anchorIds[i]);
        arap.setConstraint(vh, mesh.point(vh));
    }

    for (int i = 0; i < nids; ++i) {
        //deform::Mesh::VertexHandle vh = mesh.vertex_handle(handleIds[i]);
        deform::Mesh::VertexHandle vh = mesh.vertex_handle(32);
        arap.setConstraint(vh, mesh.point(vh) + deform::Mesh::Point(0, 0, 1.0));
    }

    arap.deform(10);   

    std::cout << "ok" << std::endl;

    if (!OpenMesh::IO::write_mesh(mesh, "deformed.ply")) {
        std::cerr << "Failed to write mesh" << std::endl;
        return -1;
    }
    

    return 0;

}
