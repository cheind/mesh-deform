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

    inline Eigen::Vector3f toEigen(const Mesh::Point &p) {
        return Eigen::Vector3f(p[0], p[1], p[2]);
    }

    inline Mesh::Point toMesh(Eigen::Ref<const Eigen::Vector3f> p) {
        return Mesh::Point(p.x(), p.y(), p.z());
    }

    typedef Mesh::Scalar Scalar;
    typedef OpenMesh::VPropHandleT<bool> VertexBoolProperty;
    typedef OpenMesh::VPropHandleT<Mesh::Point> VertexPointProperty;

    struct AsRigidAsPossibleDeformation::data {
        Mesh *mesh;
        VertexBoolProperty isConstrained;
        VertexPointProperty meshPoints; // target

        Eigen::Matrix3Xf points;
        std::vector<Eigen::Matrix3f, Eigen::aligned_allocator<Eigen::Matrix3f> > rotations;
        Eigen::Matrix<float, Eigen::Dynamic, 3> b;
        Eigen::SparseMatrix<float> L;
        Eigen::SimplicialLDLT< Eigen::SparseMatrix<float> > solver;
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
        copyPositions();
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
        _data->mesh->property(_data->meshPoints, v) = pos;
    }

    void AsRigidAsPossibleDeformation::deform(size_t nIterations)
    {
        if (_data->dirty) {
            const size_t nv = _data->mesh->n_vertices();
            
            // Initial rotations are identity
            _data->rotations.clear();
            _data->rotations.resize(nv, Eigen::Matrix3f::Identity());

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

    void AsRigidAsPossibleDeformation::copyPositions()
    {
        Mesh *m = _data->mesh;

        Eigen::Matrix3Xf &points = _data->points;
        points.resize(3, m->n_vertices());

        for (size_t i = 0; i < m->n_vertices(); ++i) {
            points.col(i) = toEigen(m->point(Mesh::VertexHandle(static_cast<int>(i))));
        }
    }

    void AsRigidAsPossibleDeformation::restorePositions()
    {
        Mesh *m = _data->mesh;

        Eigen::Matrix3Xf &points = _data->points;

        for (size_t i = 0; i < m->n_vertices(); ++i) {
            m->point(Mesh::VertexHandle(static_cast<int>(i))) = toMesh(points.col(i));
        }
    }

    void AsRigidAsPossibleDeformation::attachMeshProperties()
    {
        Mesh *m = _data->mesh;

        VertexBoolProperty &vc = _data->isConstrained;
        m->add_property(vc);
        std::fill(m->property(vc).data_vector().begin(), m->property(vc).data_vector().end(), false);

        _data->meshPoints = m->points_pph();
    }

    void AsRigidAsPossibleDeformation::releaseMeshProperties()
    {
        Mesh *m = _data->mesh;
        if (!m)
            return;

        m->remove_property(_data->isConstrained);
    }

    void AsRigidAsPossibleDeformation::estimateRotations()
    {
        Mesh &m = *_data->mesh;
        for (auto viter = m.vertices_begin(); viter != m.vertices_end(); ++viter) {
            const Eigen::Vector3f vi = toEigen(m.point(*viter));
            const Eigen::Vector3f vip = _data->points.col(viter->idx());

            Eigen::Matrix3f cov = Eigen::Matrix3f::Zero();

            for (auto vh = m.voh_begin(*viter); vh != m.voh_end(*viter); ++vh) {
                
                Mesh::VertexHandle vhj = m.to_vertex_handle(*vh);
                const float w = -_data->L.coeff(viter->idx(), vhj.idx());

                const Eigen::Vector3f vj = toEigen(m.point(vhj));
                const Eigen::Vector3f &vjp = _data->points.col(vhj.idx());

                cov += w * (vi - vj) * ((vip - vjp).transpose());
            }

            Eigen::JacobiSVD<Eigen::Matrix3f> svd(cov, Eigen::ComputeFullU | Eigen::ComputeFullV);
            const Eigen::Matrix3f &v = svd.matrixV();
            const Eigen::Matrix3f ut = svd.matrixU().transpose();
            
            Eigen::Matrix3f id = Eigen::Matrix3f::Identity(3, 3);
            id(2, 2) = (v * ut).determinant();
            _data->rotations[viter->idx()] = (v * id * ut); // accounts for required flip
        }
    }

    void AsRigidAsPossibleDeformation::estimatePositions()
    {
        computeB();
        
        // Solve Lp' = b
        for (size_t i = 0; i < 3; ++i) {
            _data->points.row(i) = _data->solver.solve(_data->b.col(i)).transpose();
        }        
    }

    void AsRigidAsPossibleDeformation::computeL()
    {
        Mesh &m = *_data->mesh;
        
        OpenMesh::HPropHandleT<float> hew;
        m.add_property(hew);
        cotanWeights(m, hew);
        cotanMatrix(m, hew, _data->L);
        m.remove_property(hew);
        
        _data->solver.compute(_data->L);
    }

    void AsRigidAsPossibleDeformation::computeB()
    {
        Mesh &m = *_data->mesh;

        _data->b.setZero();

        for (auto viter = m.vertices_begin(); viter != m.vertices_end(); ++viter) {
            int i = viter->idx();

            if (m.property(_data->isConstrained, *viter)) {
                // If this vertex is constraint we enforce lhs = rhs
                _data->b.row(i) = toEigen(m.point(*viter)).transpose();
            }
            else {
                // else w/2 * (Ri+Rj) * (pi - pj)
                Eigen::Vector3f sum = Eigen::Vector3f::Zero();

                for (auto vh = m.voh_begin(*viter); vh != m.voh_end(*viter); ++vh) {
                    
                    Mesh::VertexHandle vj = m.to_vertex_handle(*vh);
                    
                    const float w = -_data->L.coeff(i, vj.idx());

                    Eigen::Vector3f t = (_data->rotations[i] + _data->rotations[vj.idx()]) * (toEigen(m.point(*viter)) - toEigen(m.point(vj)));

                    sum += 0.5f * w * t;
                }
                _data->b.row(i) = sum.transpose();
            }
        }
    }
}