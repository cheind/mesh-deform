/**
 This file is part of mesh-deform.
 
 Copyright(C) 2016 Christoph Heindl
 All rights reserved.
 
 This software may be modified and distributed under the terms
 of the BSD license.See the LICENSE file for details.
 */

#include <deform/cotan_matrix.h>
#include <iostream>

namespace deform {
    
    void cotanWeights(Mesh &m, OpenMesh::HPropHandleT<float> &w) {
        // https://igl.ethz.ch/projects/bbw/a-cotangent-laplacian-for-images-as-surfaces-2012-jacobson-sorkine.pdf
        
        std::fill(m.property(w).data_vector().begin(), m.property(w).data_vector().end(), 0.f);
        
        // Compute edge lengths
        OpenMesh::EPropHandleT<float> edgelengths;
        m.add_property(edgelengths);
        
        for (auto e = m.edges_begin(); e != m.edges_end(); ++e) {
            Mesh::HalfedgeHandle heh0 =  m.halfedge_handle(*e, 0);
            m.property(edgelengths, *e) = (m.point(m.from_vertex_handle(heh0)) - m.point(m.to_vertex_handle(heh0))).length();
        }
        
        // Compute weights
        for (auto f = m.faces_begin(); f != m.faces_end(); ++f) {
            
            Mesh::HalfedgeHandle he0 = m.halfedge_handle(*f);
            Mesh::HalfedgeHandle he1 = m.next_halfedge_handle(he0);
            Mesh::HalfedgeHandle he2 = m.next_halfedge_handle(he1);
            
            const float l0 = m.property(edgelengths, m.edge_handle(he0));
            const float l1 = m.property(edgelengths, m.edge_handle(he1));
            const float l2 = m.property(edgelengths, m.edge_handle(he2));
            
            const float semip = 0.5f * (l0 + l1 + l2);
            const float area = std::sqrt(semip*(semip - l0)*(semip - l1) * (semip - l2));
            const float denom = 1.f / (4.f * area);
            
            
            const float cot0 = (-l0 * l0 + l1 * l1 + l2 * l2) * denom;
            const float cot1 = (l0 * l0 - l1 * l1 + l2 * l2) * denom;
            const float cot2 = (l0 * l0 + l1 * l1 - l2 * l2) * denom;
            

            m.property(w, he0) = cot0;
            m.property(w, he1) = cot1;
            m.property(w, he2) = cot2;
        }

        m.remove_property(edgelengths);
    }
    
    void cotanMatrix(Mesh &m, const OpenMesh::HPropHandleT<float> &w, Eigen::SparseMatrix<float> &c) {
        const size_t nverts = m.n_vertices();
        
        c.resize(nverts, nverts);
        c.reserve(Eigen::VectorXi::Constant(nverts, 1, 7)); // Average numbers of neighbors on closed 2d manifold
        
        typedef Eigen::Triplet<float> T;
        std::vector<T> triplets;
        
        for (auto v = m.vertices_begin(); v != m.vertices_end(); ++v) {
            for (auto h = m.voh_begin(*v); h != m.voh_end(*v); ++h) {
                Mesh::HalfedgeHandle opp = m.opposite_halfedge_handle(*h);
                
                const float weight = (m.property(w, *h) + m.property(w, opp)) / 2;
                
                triplets.push_back(T(v->idx(), m.to_vertex_handle(*h).idx(), -weight));
                triplets.push_back(T(v->idx(), v->idx(), weight));
            }
        }
        
        c.setZero();
        c.setFromTriplets(triplets.begin(), triplets.end());
    }
    
    
}