/**
 This file is part of mesh-deform.
 
 Copyright(C) 2016 Christoph Heindl
 All rights reserved.
 
 This software may be modified and distributed under the terms
 of the BSD license.See the LICENSE file for details.
 */

#ifndef DEFORM_ARAP_H
#define DEFORM_ARAP_H

#include <deform/mesh.h>
#include <memory>

namespace deform {
    
    class AsRigidAsPossibleDeform {
    public:
        AsRigidAsPossibleDeform(Mesh * mesh);
        ~AsRigidAsPossibleDeform();
        AsRigidAsPossibleDeform(const AsRigidAsPossibleDeform &other) = delete;
        AsRigidAsPossibleDeform &operator=(const AsRigidAsPossibleDeform &other) = delete;

        void setConstraint(Mesh::VertexHandle v, const Mesh::Point &pos);
        void deform(size_t nIterations);

    private:

        void copyPositions();
        void restorePositions();
        void attachMeshProperties();
        void releaseMeshProperties();
        void computeEdgeWeights();
        void estimateRotations();
        void estimatePositions();
        void computeL();
        void computeB();


        struct data;
        std::unique_ptr<data> _data;
    };

}

#endif