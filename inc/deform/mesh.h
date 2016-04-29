/**
 This file is part of mesh-deform.
 
 Copyright(C) 2016 Christoph Heindl
 All rights reserved.
 
 This software may be modified and distributed under the terms
 of the BSD license.See the LICENSE file for details.
 */

#ifndef DEFORM_MESH_H
#define DEFORM_MESH_H

#ifndef _USE_MATH_DEFINES
    #define _USE_MATH_DEFINES
#endif

#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>

namespace deform {

    typedef OpenMesh::TriMesh_ArrayKernelT<> Mesh;
       
}

#endif