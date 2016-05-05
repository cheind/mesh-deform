/**
 This file is part of mesh-deform.
 
 Copyright(C) 2016 Christoph Heindl
 All rights reserved.
 
 This software may be modified and distributed under the terms
 of the BSD license.See the LICENSE file for details.
 */

#ifndef DEFORM_MESH_H
#define DEFORM_MESH_H

#ifndef _USE_MATH_DEFINES
    #define _USE_MATH_DEFINES
#endif

#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
#include <Eigen/Core>

namespace deform {
    
    
    typedef OpenMesh::TriMesh_ArrayKernelT<> Mesh;

    
    inline Eigen::Vector3f toEigenF(const Mesh::Point &p) {
        return Eigen::Vector3f(p[0], p[1], p[2]);
    }
    
    inline Mesh::Point toMeshF(Eigen::Ref<const Eigen::Vector3f> p) {
        return Mesh::Point(Mesh::Scalar(p.x()), Mesh::Scalar(p.y()), Mesh::Scalar(p.z()));
    }

    
    class OpenMeshAdapter {
    public:
        typedef float Scalar;
        
        OpenMeshAdapter(Mesh &m)
        :_m(m) {}
        
        Eigen::Matrix<Scalar, 3, Eigen::Dynamic> vertices() const {
            
            Eigen::Matrix<Scalar, 3, Eigen::Dynamic> points(3, numberOfVertices());
            for (int i = 0; i < numberOfVertices(); ++i) {
                points.col(i) = toEigenF(_m.point(Mesh::VertexHandle(i)));
            }
            
            return points;
        }
        
        void vertices(const Eigen::Matrix<Scalar, 3, Eigen::Dynamic> &p) {
            for (size_t i = 0; i < numberOfVertices(); ++i) {
                _m.point(Mesh::VertexHandle(static_cast<int>(i))) = toMeshF(p.col(i));
            }

        }
        
        Eigen::Matrix<int, 3, Eigen::Dynamic> faces() const {
            Eigen::Matrix<int, 3, Eigen::Dynamic> tris(3, numberOfFaces());
            
            
            for (auto f = _m.faces_begin(); f != _m.faces_end(); ++f) {
                auto v = _m.fv_begin(*f);
                tris(0, f->idx()) = v->idx(); ++v;
                tris(1, f->idx()) = v->idx(); ++v;
                tris(2, f->idx()) = v->idx();
                
            }
            
            return tris;
        }
        
        int numberOfFaces() const {
            return static_cast<int>(_m.n_faces());
        }
        
        int numberOfVertices() const {
            return static_cast<int>(_m.n_vertices());
        }
        
    private:
        Mesh &_m;

    };
    
}

#endif