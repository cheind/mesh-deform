/**
 This file is part of mesh-deform.
 
 Copyright(C) 2016 Christoph Heindl
 All rights reserved.
 
 This software may be modified and distributed under the terms
 of the BSD license.See the LICENSE file for details.
 */

#include <deform/arap.h>

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/StdVector>

namespace deform {

    inline Eigen::Vector3f toEigen(const Mesh::Point &p) {
        return Eigen::Vector3f(p[0], p[1], p[2]);
    }

    inline Mesh::Point toMesh(Eigen::Ref<const Eigen::Vector3f> p) {
        return Mesh::Point(p.x(), p.y(), p.z());
    }

    typedef Mesh::Scalar Scalar;
    typedef OpenMesh::EPropHandleT<Mesh::Scalar> EdgeWeightProperty;
    typedef OpenMesh::VPropHandleT<bool> VertexBoolProperty;
    typedef OpenMesh::VPropHandleT<Mesh::Point> VertexPointProperty;

    struct AsRigidAsPossibleDeform::data {
        Mesh *mesh;
        EdgeWeightProperty weights;
        VertexBoolProperty isConstrained;
        VertexPointProperty meshPoints;

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

    AsRigidAsPossibleDeform::AsRigidAsPossibleDeform(Mesh *mesh)
        :_data(new data())
    {
        if (!mesh) {
            throw std::runtime_error("Mesh cannot be null.");
        }

        _data->mesh = mesh;
        attachMeshProperties();
        copyPositions();        
        computeEdgeWeights();
    }

   
    AsRigidAsPossibleDeform::~AsRigidAsPossibleDeform()
    {
        releaseMeshProperties();
    }
   
    void AsRigidAsPossibleDeform::setConstraint(Mesh::VertexHandle v, const Mesh::Point & pos)
    {
        bool isc = _data->mesh->property(_data->isConstrained, v);
        _data->dirty |= !isc;
        _data->mesh->property(_data->isConstrained, v) = true;
        _data->mesh->property(_data->)
        _data->points[1].col(v.idx()) = toEigen(pos);
    }

    void AsRigidAsPossibleDeform::deform(size_t nIterations)
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
    }

    void AsRigidAsPossibleDeform::copyPositions()
    {
        Mesh *m = _data->mesh;

        Eigen::Matrix3Xf &points = _data->points;
        points.resize(3, m->n_vertices());

        for (size_t i = 0; i < m->n_vertices(); ++i) {
            points.col(i) = toEigen(m->point(Mesh::VertexHandle(static_cast<int>(i))));
        }
    }

    void AsRigidAsPossibleDeform::attachMeshProperties()
    {
        Mesh *m = _data->mesh;

        EdgeWeightProperty &w = _data->weights;
        m->add_property(w);
        std::fill(m->property(w).data_vector().begin(), m->property(w).data_vector().end(), Scalar(1.f));

        VertexBoolProperty &vc = _data->isConstrained;
        m->add_property(vc);
        std::fill(m->property(vc).data_vector().begin(), m->property(vc).data_vector().end(), false);

        _data->meshPoints = m->points_pph();
    }

    void AsRigidAsPossibleDeform::releaseMeshProperties()
    {
        Mesh *m = _data->mesh;
        if (!m)
            return;

        m->remove_property(_data->weights);
        m->remove_property(_data->isConstrained);
    }

    void AsRigidAsPossibleDeform::computeEdgeWeights()
    {
        Mesh &m = *_data->mesh;
        Eigen::Matrix3Xf &pos = _data->points;
        EdgeWeightProperty &w = _data->weights;

        // Compute cotangent weights for each edge.
        for (auto eiter = m.edges_begin(); eiter != m.edges_end(); ++eiter) {
            
            float cot_alpha = 0.f;
            float cot_beta = 0.f;

            Mesh::HalfedgeHandle he0 = m.halfedge_handle(*eiter, 0);
            if (he0.is_valid()) {
                // Face present
                Mesh::HalfedgeHandle he1 = m.next_halfedge_handle(he0);
                Mesh::HalfedgeHandle he2 = m.next_halfedge_handle(he1);
                
                Eigen::Vector3f v0 = pos.col(m.from_vertex_handle(he1).idx()) - pos.col(m.to_vertex_handle(he1).idx());
                Eigen::Vector3f v1 = pos.col(m.to_vertex_handle(he2).idx()) - pos.col(m.from_vertex_handle(he2).idx());

                // cot(alpha) = cos(alpha) / sin(alpha)
                // (dot(v0,v1) * cos(alpha)) / (|v0||v1| * sin(alpha))
                // (dot(v0,v1) / |cross(v0, v1)|
                cot_alpha = v0.dot(v1) / (1e-6f * v0.cross(v1).norm());
            }

            he0 = m.halfedge_handle(*eiter, 1);
            if (he0.is_valid()) {
                // Face present
                Mesh::HalfedgeHandle he1 = m.next_halfedge_handle(he0);
                Mesh::HalfedgeHandle he2 = m.next_halfedge_handle(he1);

                Eigen::Vector3f v0 = pos.col(m.from_vertex_handle(he1).idx()) - pos.col(m.to_vertex_handle(he1).idx());
                Eigen::Vector3f v1 = pos.col(m.to_vertex_handle(he2).idx()) - pos.col(m.from_vertex_handle(he2).idx());

                cot_beta = v0.dot(v1) / (1e-6f * v0.cross(v1).norm());
            }

            m.property(w, *eiter) = (cot_alpha + cot_beta) * 0.5f;
        }        

    }

    void AsRigidAsPossibleDeform::estimateRotations()
    {
        Mesh &m = *_data->mesh;
        for (auto viter = m.vertices_begin(); viter != m.vertices_end(); ++viter) {
            const Eigen::Vector3f &vi = _data->points[0].col(viter->idx());
            const Eigen::Vector3f &vip = _data->points[1].col(viter->idx());

            Eigen::Matrix3f cov = Eigen::Matrix3f::Zero();

            for (auto vv = m.vv_begin(*viter); vv; ++vv) {
                Mesh::HalfedgeHandle heh = vv.current_halfedge_handle();
                Mesh::EdgeHandle eh = m.edge_handle(heh);
                float w = m.property(_data->weights, eh);

                const Eigen::Vector3f &vj = _data->points[0].col(vv->idx());
                const Eigen::Vector3f &vjp = _data->points[1].col(vv->idx());

                cov += (vi - vj) * (vip - vjp).transpose() * w;
            }

            Eigen::JacobiSVD<Eigen::Matrix3f> svd(cov, Eigen::ComputeThinU | Eigen::ComputeThinV); 
            const Eigen::Matrix3f &v = svd.matrixV();
            const Eigen::Matrix3f ut = svd.matrixU().transpose();
            
            Eigen::Matrix3f id = Eigen::Matrix3f::Identity(3, 3);
            id(2, 2) = (v * ut).determinant();
            _data->rotations[viter->idx()] = (v * id * ut); // accounts for required flip
        }
    }

    void AsRigidAsPossibleDeform::estimatePositions()
    {
        computeB();
        
        // Solve Lp' = b
        for (size_t i = 0; i < 3; ++i) {
            _data->points[0].row(i) = _data->solver.solve(_data->b.col(i)).transpose();
        }        
    }

    void AsRigidAsPossibleDeform::computeL()
    {
        Mesh &m = *_data->mesh;
       
        typedef Eigen::Triplet<float> T;
        std::vector<T> triplets;

        for (auto viter = m.vertices_begin(); viter != m.vertices_end(); ++viter) {
            int i = viter->idx();
            
            if (m.property(_data->isConstrained, *viter)) {
                // If this vertex is constraint we enforce lhs = rhs
                triplets.push_back(T(i, i, 1.f));              
            } else {
                // Otherwise we need to loop over connected vertices.
                float sumw = 0.0f;
                for (auto vv = m.vv_begin(*viter); vv; ++vv) {
                    Mesh::HalfedgeHandle heh = vv.current_halfedge_handle();
                    Mesh::EdgeHandle eh = m.edge_handle(heh);
                    float w = m.property(_data->weights, eh);

                    sumw += w;
                    triplets.push_back(T(i, vv->idx(), -w));
                }
                triplets.push_back(T(i, i, sumw));
            }

        }

        _data->L.setZero();
        _data->L.setFromTriplets(triplets.begin(), triplets.end());

        _data->solver.analyzePattern(_data->L);
        _data->solver.compute(_data->L);
    }

    void AsRigidAsPossibleDeform::computeB()
    {
        Mesh &m = *_data->mesh;

        _data->b.setZero();

        for (auto viter = m.vertices_begin(); viter != m.vertices_end(); ++viter) {
            int i = viter->idx();

            if (m.property(_data->isConstrained, *viter)) {
                // If this vertex is constraint we enforce lhs = rhs
                _data->b.row(i) = _data->points[1].col(i).transpose();
            }
            else {
                // else w/2 * (Ri+Rj) * (pi - pj)
                Eigen::Vector3f sum = Eigen::Vector3f::Zero();

                for (auto vv = m.vv_begin(*viter); vv; ++vv) {
                    int j = vv->idx();

                    Mesh::HalfedgeHandle heh = vv.current_halfedge_handle();
                    Mesh::EdgeHandle eh = m.edge_handle(heh);
                    float w = m.property(_data->weights, eh);

                    sum += 0.5f * w * (_data->rotations[i] + _data->rotations[j]) * (_data->points[0].col(i) - _data->points[0].col(j));
                }
                _data->b.row(i) = sum.transpose();
            }
        }
    }
}