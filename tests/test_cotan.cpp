/**
 This file is part of mesh-deform.
 
 Copyright(C) 2016 Christoph Heindl
 All rights reserved.
 
 This software may be modified and distributed under the terms
 of the BSD license.See the LICENSE file for details.
*/

#define CATCH_CONFIG_MAIN
#include "catch.hpp"

#include <deform/arap.h>
#include <deform/openmesh_adapter.h>
#include <iostream>

#include "accessor.h"

TEST_CASE("cotan_weights")
{
    typedef OpenMesh::TriMesh_ArrayKernelT<> Mesh;
    
    Mesh m;
    
    m.add_vertex(Mesh::Point(0.f,0.f,0.f));
    m.add_vertex(Mesh::Point(1.f,0.f,0.f));
    m.add_vertex(Mesh::Point(1.f,1.f,0.f));
    m.add_vertex(Mesh::Point(0.f,1.f,0.f));
    
    m.add_face(Mesh::VertexHandle(0), Mesh::VertexHandle(1), Mesh::VertexHandle(2));
    m.add_face(Mesh::VertexHandle(0), Mesh::VertexHandle(2), Mesh::VertexHandle(3));
    
    typedef deform::AsRigidAsPossibleDeformation< deform::OpenMeshAdapter<> > ARAP;
    
    deform::OpenMeshAdapter<> adapter(m);
    
    ARAP arap(adapter);
    arap.deform(0);
    
    Eigen::MatrixXf sp = deform::PrivateAccessor<ARAP>::cotanWeights(arap);
    
    REQUIRE(sp.rows() == 4);
    REQUIRE(sp.cols() == 4);
    
    Eigen::Matrix4f expected;
    expected <<
    0.f, 0.5f, 0.f, 0.5f,
    0.5f, 0.f, 0.5f, 0.f,
    0.f, 0.5f, 0.f, 0.5f,
    0.5f, 0.f, 0.5f, 0.f;
    
    REQUIRE(sp.isApprox(expected, 1e-4f));
}
    