/**
 This file is part of mesh-deform.
 
 Copyright(C) 2016 Christoph Heindl
 All rights reserved.
 
 This software may be modified and distributed under the terms
 of the BSD license.See the LICENSE file for details.
 */

#ifndef DEFORM_COTAN_MATRIX_H
#define DEFORM_COTAN_MATRIX_H

#include <deform/mesh.h>
#include <Eigen/Sparse>

namespace deform {
    
    void cotanWeights(Mesh &m, OpenMesh::EPropHandleT<double> &w);
    void cotanMatrix(Mesh &m, const OpenMesh::EPropHandleT<double> &w, Eigen::SparseMatrix<double> &c);
    void cotanMatrixWithConstraints(Mesh &m, const OpenMesh::EPropHandleT<double> &w, const OpenMesh::VPropHandleT<bool> &vc, Eigen::SparseMatrix<double> &c);
    
}

#endif