/**
 This file is part of mesh-deform.
 
 Copyright(C) 2016 Christoph Heindl
 All rights reserved.
 
 This software may be modified and distributed under the terms
 of the BSD license.See the LICENSE file for details.
 */


#include <deform/arap.h>
#include <iostream>

int main(int argc, char **argv) {

    if (argc != 4) {
        std::cerr << argv[0] << " plane.ply anchors.txt handles.txt" << std::endl;
        return -1;
    }

    deform::Mesh mesh;
    if (!OpenMesh::IO::read_mesh(mesh, argv[1])) {
        std::cerr << "Failed to read mesh" << std::endl;
        return -1;
    }

    deform::AsRigidAsPossibleDeform arap(&mesh);
    
}
