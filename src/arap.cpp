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
    typedef OpenMesh::VPropHandleT<int> VertexIntProperty;
    typedef OpenMesh::VPropHandleT<Eigen::Vector3d> VertexPointProperty;
    typedef OpenMesh::EPropHandleT<double> EdgeWeightProperty;

    struct AsRigidAsPossibleDeformation::data {
        Mesh *mesh;

        VertexBoolProperty isConstrained;
        VertexIntProperty freeIdx;
        VertexPointProperty constraintLocations; // target
        EdgeWeightProperty cotanWeights;

        Eigen::Matrix3Xd points, pointsPrime;

        std::vector<Eigen::Matrix3d, Eigen::aligned_allocator<Eigen::Matrix3d> > rotations;
        Eigen::SparseMatrix<double> L;
        Eigen::Matrix3Xd b;
        int nFreeVariables;
        
        Eigen::SimplicialLDLT<Eigen::SparseMatrix<double> > solver;

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

        std::fill(mesh->property(_data->isConstrained).data_vector().begin(), mesh->property(_data->isConstrained).data_vector().end(), false);
        
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

    void AsRigidAsPossibleDeformation::attachMeshProperties()
    {
        Mesh *m = _data->mesh;
        
        m->add_property(_data->isConstrained);
        m->add_property(_data->cotanWeights);
        m->add_property(_data->constraintLocations);
        m->add_property(_data->freeIdx);
    }

    void AsRigidAsPossibleDeformation::releaseMeshProperties()
    {
        Mesh *m = _data->mesh;
        if (!m)
            return;

        m->remove_property(_data->isConstrained);
        m->remove_property(_data->cotanWeights);
        m->remove_property(_data->constraintLocations);
        m->remove_property(_data->freeIdx);
    }

    void AsRigidAsPossibleDeformation::restorePositions()
    {
        Mesh *m = _data->mesh;

        Eigen::Matrix3Xd &p = _data->pointsPrime;

        for (size_t i = 0; i < m->n_vertices(); ++i) {
            m->point(Mesh::VertexHandle(static_cast<int>(i))) = toMesh(p.col(i));
        }
    }

    void AsRigidAsPossibleDeformation::deform(size_t nIterations)
    {
        if (_data->dirty) {
            const size_t nv = _data->mesh->n_vertices();

            // Initial rotations for all variables including constrained are identity
            _data->rotations.clear();
            _data->rotations.resize(nv, Eigen::Matrix3d::Identity());

            // Compute number of free variables and mapping of indices.
            this->setupFreeVariableMap();

            // Pre-factor matrix A (L)        
            prepareLinearSystem();

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

    void AsRigidAsPossibleDeformation::prepareLinearSystem()
    {
        Mesh &m = *_data->mesh;
        Eigen::SparseMatrix<double> &L = _data->L;
        Eigen::Matrix3Xd &b = _data->b;

        // Compute edge weights.
        cotanWeights(m, _data->cotanWeights);

        // Setup L and constant part of b.
        L.resize(_data->nFreeVariables, _data->nFreeVariables);
        L.reserve(Eigen::VectorXi::Constant(_data->nFreeVariables, 1, 7)); // Average numbers of neighbors on closed 2d manifold
        L.setZero();

        b.resize(3, _data->nFreeVariables);
        b.setZero();

        typedef Eigen::Triplet<double> T;
        std::vector<T> triplets;

        for (auto vi = m.vertices_begin(); vi != m.vertices_end(); ++vi) {
            
            // Skip constrained vertex rows and cols
            if (m.property(_data->isConstrained, *vi))
                continue;

            int iidx = m.property(_data->freeIdx, *vi);

            // If not constrained setup cotan matrix.
            for (auto h = m.voh_begin(*vi); h != m.voh_end(*vi); ++h) {
                Mesh::VertexHandle vj = m.to_vertex_handle(*h);

                int jidx = m.property(_data->freeIdx, vj);
                
                Mesh::EdgeHandle eh = m.edge_handle(*h);
                const double weight = m.property(_data->cotanWeights, eh);

                if (m.property(_data->isConstrained, vj)) {
                    b.col(iidx) += weight * m.property(_data->constraintLocations, vj);
                } else {
                    triplets.push_back(T(iidx, jidx, -weight));
                }

                triplets.push_back(T(iidx, iidx, weight));
            }
        }

        L.setFromTriplets(triplets.begin(), triplets.end());
        _data->solver.compute(L);
        
    }

    void AsRigidAsPossibleDeformation::setupFreeVariableMap()
    {
        Mesh &m = *_data->mesh;                
        VertexIntProperty &map = _data->freeIdx;
        VertexBoolProperty &constrained = _data->isConstrained;

        int freeIdx = 0;
        for (auto v = m.vertices_begin(); v != m.vertices_end(); ++v) {
            int id = m.property(constrained, *v) ? -1 : freeIdx++;
            m.property(map, *v) = id;
        }

        _data->nFreeVariables = freeIdx;
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
        Mesh &m = *_data->mesh;

        Eigen::Matrix3Xd rhs = _data->b;

        for (auto vi = m.vertices_begin(); vi != m.vertices_end(); ++vi) {
            if (m.property(_data->isConstrained, *vi))
                continue;

            int iidx = m.property(_data->freeIdx, *vi);

            // else w/2 * (Ri+Rj) * (pi - pj)
            for (auto vh = m.voh_begin(*vi); vh != m.voh_end(*vi); ++vh) {

                Mesh::VertexHandle vj = m.to_vertex_handle(*vh);
                int j = vj.idx();

                Eigen::Matrix3d r = (_data->rotations[vi->idx()] + _data->rotations[vj.idx()]);
                Eigen::Vector3d v = _data->points.col(vi->idx()) - _data->points.col(vj.idx());
                Eigen::Vector3d t = r * v;
                const double w = m.property(_data->cotanWeights, m.edge_handle(*vh));
                rhs.col(iidx) += w * t * 0.5;
            }
        }
        
        // Solve Lp' = b

        for (int d = 0; d < 3; ++d) {
            Eigen::VectorXd pprime = _data->solver.solve(rhs.row(d).transpose());

            int idx = 0;
            for (auto vi = m.vertices_begin(); vi != m.vertices_end(); ++vi) {
                if (!m.property(_data->isConstrained, *vi)) {
                    _data->pointsPrime(d, vi->idx()) = pprime(idx++);
                }                    
            }
        }
    }


    void AsRigidAsPossibleDeformation::computeB()
    {
        /*
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
        */
    }
}