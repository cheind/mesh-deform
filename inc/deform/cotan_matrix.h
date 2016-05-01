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
    
    void cotanWeights(Mesh &m, OpenMesh::HPropHandleT<float> &w);
    void cotanMatrix(Mesh &m, const OpenMesh::HPropHandleT<float> &w, Eigen::SparseMatrix<float> &c);
    
}

#endif