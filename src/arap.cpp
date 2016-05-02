/**
 This file is part of mesh-deform.
 
 Copyright(C) 2016 Christoph Heindl
 All rights reserved.
 
 This software may be modified and distributed under the terms
 of the BSD license.See the LICENSE file for details.
 */

#include <deform/arap.h>
#include <deform/cotan_matrix.h>

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/StdVector>
#include <iostream>

namespace deform {

    inline Eigen::Vector3d toEigen(const Mesh::Point &p) {
        return Eigen::Vector3d(p[0], p[1], p[2]);
    }

    inline Mesh::Point toMesh(Eigen::Ref<const Eigen::Vector3d> p) {
        return Mesh::Point(Mesh::Scalar(p.x()), Mesh::Scalar(p.y()), Mesh::Scalar(p.z()));
    }

    typedef OpenMesh::VPropHandleT<bool> VertexBoolProperty;
    typedef OpenMesh::VPropHandleT<Eigen::Vector3d> VertexPointProperty;
    typedef OpenMesh::EPropHandleT<double> EdgeWeightProperty;

    struct AsRigidAsPossibleDeformation::data {
        Mesh *mesh;

        VertexBoolProperty isConstrained;
        VertexPointProperty constraintLocations; // target
        EdgeWeightProperty cotanWeights;

        Eigen::Matrix3Xd points, pointsPrime;

        std::vector<Eigen::Matrix3d, Eigen::aligned_allocator<Eigen::Matrix3d> > rotations;
        Eigen::SparseMatrix<double> Lt;
        Eigen::Matrix<double, Eigen::Dynamic, 3> b;
        
        //Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
        Eigen::SimplicialLDLT<Eigen::SparseMatrix<double> > solver;
        //Eigen::SimplicialLLT< Eigen::SparseMatrix<float> > solver;
        bool dirty;

        data()
            :mesh(0), dirty(true)
        {}
    };

    AsRigidAsPossibleDeformation::AsRigidAsPossibleDeformation(Mesh *mesh)
        :_data(new data())
    {
        if (!mesh) {
            throw std::runtime_error("Mesh cannot be null.");
        }

        _data->mesh = mesh;
        attachMeshProperties();        
        
        const size_t nverts = mesh->n_vertices();
        _data->points.resize(3, nverts);
        for (size_t i = 0; i < nverts; ++i) {
            _data->points.col(i) = toEigen(mesh->point(Mesh::VertexHandle(static_cast<int>(i))));
        }
        _data->pointsPrime = _data->points;
    }

   
    AsRigidAsPossibleDeformation::~AsRigidAsPossibleDeformation()
    {
        releaseMeshProperties();
    }
   
    void AsRigidAsPossibleDeformation::setConstraint(Mesh::VertexHandle v, const Mesh::Point & pos)
    {
        bool isc = _data->mesh->property(_data->isConstrained, v);
        _data->dirty |= !isc;
        _data->mesh->property(_data->isConstrained, v) = true;
        _data->mesh->property(_data->constraintLocations, v) = toEigen(pos);
        _data->pointsPrime.col(v.idx()) = toEigen(pos);
    }

    void AsRigidAsPossibleDeformation::deform(size_t nIterations)
    {
        if (_data->dirty) {
            const size_t nv = _data->mesh->n_vertices();
            
            // Initial rotations are identity
            _data->rotations.clear();
            _data->rotations.resize(nv, Eigen::Matrix3d::Identity());

            // Init right hand side of linear system
            _data->b.resize(nv, 3);
            _data->b.setZero();

            // Pre-factor matrix A (L)        
            computeL();

            _data->dirty = false;
        }

        for (size_t i = 0; i < nIterations; ++i) {
            // Estimate rotations
            estimateRotations();

            // Estimate positions
            estimatePositions();
        }

        restorePositions();
    }

    void AsRigidAsPossibleDeformation::restorePositions()
    {
        Mesh *m = _data->mesh;

        Eigen::Matrix3Xd &p = _data->pointsPrime;

        for (size_t i = 0; i < m->n_vertices(); ++i) {
            m->point(Mesh::VertexHandle(static_cast<int>(i))) = toMesh(p.col(i));
        }
    }

    void AsRigidAsPossibleDeformation::attachMeshProperties()
    {
        Mesh *m = _data->mesh;

        VertexBoolProperty &vc = _data->isConstrained;
        m->add_property(vc);
        std::fill(m->property(vc).data_vector().begin(), m->property(vc).data_vector().end(), false);

        EdgeWeightProperty &ew = _data->cotanWeights;
        m->add_property(ew);

        VertexPointProperty &vcp = _data->constraintLocations;
        m->add_property(vcp);
    }

    void AsRigidAsPossibleDeformation::releaseMeshProperties()
    {
        Mesh *m = _data->mesh;
        if (!m)
            return;

        m->remove_property(_data->isConstrained);
        m->remove_property(_data->cotanWeights);
        m->remove_property(_data->constraintLocations);
    }

    void AsRigidAsPossibleDeformation::estimateRotations()
    {
        Mesh &m = *_data->mesh;
        for (auto viter = m.vertices_begin(); viter != m.vertices_end(); ++viter) {

            int i = viter->idx();

            const auto &vi = _data->points.col(i);
            const auto &vip = _data->pointsPrime.col(i);

            Eigen::Matrix3d cov = Eigen::Matrix3d::Zero();

            for (auto vh = m.voh_begin(*viter); vh != m.voh_end(*viter); ++vh) {
                Mesh::VertexHandle vhj = m.to_vertex_handle(*vh);
                int j = vhj.idx();
                
                const double w = m.property(_data->cotanWeights, m.edge_handle(*vh));

                const auto &vj = _data->points.col(j);
                const auto &vjp = _data->pointsPrime.col(j);

                cov += w * (vi - vj) * ((vip - vjp).transpose());
            }

            Eigen::JacobiSVD<Eigen::Matrix3d> svd(cov, Eigen::ComputeFullU | Eigen::ComputeFullV);
            const Eigen::Matrix3d &v = svd.matrixV();
            const Eigen::Matrix3d ut = svd.matrixU().transpose();
            
            Eigen::Matrix3d id = Eigen::Matrix3d::Identity(3, 3);
            id(2, 2) = (v * ut).determinant();  // accounts for required flip
            _data->rotations[i] = (v * id * ut);
        }
    }

    void AsRigidAsPossibleDeformation::estimatePositions()
    {
        computeB();

        // Solve Lp' = b
        for (size_t i = 0; i < 3; ++i) {
            Eigen::VectorXd Ltb = _data->Lt * _data->b.col(i);
            _data->pointsPrime.row(i) = _data->solver.solve(Ltb).transpose();            
        }
    }

    void AsRigidAsPossibleDeformation::computeL()
    {
        Mesh &m = *_data->mesh;
        
        Eigen::SparseMatrix<double> L;
        cotanWeights(m, _data->cotanWeights);
        cotanMatrixWithConstraints(m, _data->cotanWeights, _data->isConstrained, L); 
        
        _data->Lt = L.transpose();
        _data->solver.compute(_data->Lt * L);        
    }

    void AsRigidAsPossibleDeformation::computeB()
    {
        Mesh &m = *_data->mesh;

        _data->b.setZero();

        for (auto viter = m.vertices_begin(); viter != m.vertices_end(); ++viter) {
            int i = viter->idx();

            if (m.property(_data->isConstrained, *viter)) {
                // If this vertex is constraint we enforce lhs = rhs
                _data->b.row(i) = m.property(_data->constraintLocations, *viter);
            }
            else {
                // else w/2 * (Ri+Rj) * (pi - pj)
                Eigen::Vector3d sum = Eigen::Vector3d::Zero();

                for (auto vh = m.voh_begin(*viter); vh != m.voh_end(*viter); ++vh) {
                    
                    Mesh::VertexHandle vj = m.to_vertex_handle(*vh);
                    int j = vj.idx();

                    Eigen::Matrix3d r = (_data->rotations[i] + _data->rotations[j]);
                    Eigen::Vector3d v = _data->points.col(i) - _data->points.col(j);
                    Eigen::Vector3d t = r * v;
                    const double w = m.property(_data->cotanWeights, m.edge_handle(*vh));
                    sum += w * t * 0.5;
                }
                _data->b.row(i) = sum.transpose();
            }
        }
    }
}