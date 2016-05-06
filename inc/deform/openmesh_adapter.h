/**
 This file is part of mesh-deform.
 
 Copyright(C) 2016 Christoph Heindl
 All rights reserved.
 
 This software may be modified and distributed under the terms
 of the BSD license.See the LICENSE file for details.
 */

#ifndef DEFORM_OPENMESH_ADAPTER_H
#define DEFORM_OPENMESH_ADAPTER_H

#ifndef _USE_MATH_DEFINES
    #define _USE_MATH_DEFINES
#endif

#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
#include <Eigen/Core>

namespace deform {
    
    namespace convert {
        
        /** 
            Helper function to convert from OpenMesh to Eigen.
        */
        template<class Scalar>
        Eigen::Matrix<Scalar, 3, 1> toEigen(const OpenMesh::VectorT<Scalar, 3> &v) {
            return Eigen::Matrix<Scalar, 3, 1>(v[0], v[1], v[2]);
        }
        
        /**
            Helper function to convert from Eigen to OpenMesh.
         */
        template<class Scalar>
        OpenMesh::VectorT<Scalar, 3> toOpenMesh(const Eigen::Matrix<Scalar, 3, 1> &v) {
            return OpenMesh::VectorT<Scalar, 3>(v.x(), v.y(), v.z());
        }
        
    }
    
    /**
        This adapter connects OpenMesh based triangular surfaces to deform algorithms.
     
     
    */
    template<class Traits = OpenMesh::DefaultTraits>
    class OpenMeshAdapter {
    public:
        typedef OpenMesh::TriMesh_ArrayKernelT<Traits> Mesh;
        
        /** Floating point precision type of mesh. Required. */
        typedef typename Mesh::Scalar Scalar;
        
        typedef typename Eigen::Matrix<Scalar, 3, 1> VertexType;
        typedef typename Eigen::Matrix<int, 3, 1> FaceType;
        
        
        /**
            Initialize from mesh. 
            
            Mesh is not copied and needs to outlive this object.
        */
        OpenMeshAdapter(Mesh &mesh)
        :_mesh(mesh)
        {}
        
        
        /** 
            Retrieve a specific vertex position. Required.
        */
        VertexType vertexLocation(int idx) const {
            return convert::toEigen(_mesh.point(_mesh.vertex_handle(idx)));
        }
        
        /**
            Set a specific vertex position. Required.
        */
        void vertexLocation(int idx, const VertexType &v) {
            _mesh.set_point(_mesh.vertex_handle(idx), convert::toOpenMesh(v));
        }
        
        /**
            Retrieve triangular vertex indices of face. Required.
         
            By convention you should return the vertices in counterclockwise order.
        */
        FaceType face(int idx) const {
            auto v = _mesh.fv_begin(_mesh.face_handle(idx));
            
            FaceType f(3, 1);
            f(0) = v->idx(); ++v;
            f(1) = v->idx(); ++v;
            f(2) = v->idx();
            
            return f;
        }
        
        /** 
            Return the total number of faces. Required.
        */
        int numberOfFaces() const {
            return static_cast<int>(_mesh.n_faces());
        }
        
        /**
            Return the total number of vertices. Required.
        */
        int numberOfVertices() const {
            return static_cast<int>(_mesh.n_vertices());
        }

        
    private:
        Mesh &_mesh;
    };
    
    
    }

#endif