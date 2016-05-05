/**
 This file is part of mesh-deform.
 
 Copyright(C) 2016 Christoph Heindl
 All rights reserved.
 
 This software may be modified and distributed under the terms
 of the BSD license.See the LICENSE file for details.
*/

#define CATCH_CONFIG_MAIN
#include "catch.hpp"

#include <deform/mesh.h>
#include <deform/cotan_matrix.h>
#include <deform/arap.h>

TEST_CASE("test_cotanmatrix")
{
    deform::Mesh m;
    
    
    
    m.add_vertex(deform::Mesh::Point(0.f,0.f,0.f));
    m.add_vertex(deform::Mesh::Point(1.f,0.f,0.f));
    m.add_vertex(deform::Mesh::Point(1.f,1.f,0.f));
    m.add_vertex(deform::Mesh::Point(0.f,1.f,0.f));
    
    m.add_face(deform::Mesh::VertexHandle(0), deform::Mesh::VertexHandle(1), deform::Mesh::VertexHandle(2));
    m.add_face(deform::Mesh::VertexHandle(0), deform::Mesh::VertexHandle(2), deform::Mesh::VertexHandle(3));
    
    
    OpenMesh::EPropHandleT<double> w;
    m.add_property(w);
    deform::cotanWeights(m, w);
    
    Eigen::SparseMatrix<double> C;
    deform::cotanMatrix(m, w, C);
    
    Eigen::Matrix4f is = C;
    Eigen::Matrix4f expected;
    expected << 1.f, -0.5f, 0.f, -0.5f,
                -0.5f, 1.f, -0.5f, 0.f,
                0.f, -0.5f, 1.f, -0.5f,
                -0.5f, 0.f, -0.5f, 1.f;
    
    REQUIRE(is.isApprox(expected, 1e-4f));
    
    
    deform::OpenMeshAdapter ma(m);
    deform::AsRigidAsPossibleDeformation2<deform::OpenMeshAdapter> arap(ma);
    arap.setConstraint(0, Eigen::Vector3f(0,0,0));
    arap.deform(1);
}
    