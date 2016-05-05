/**
 This file is part of mesh-deform.
 
 Copyright(C) 2016 Christoph Heindl
 All rights reserved.
 
 This software may be modified and distributed under the terms
 of the BSD license.See the LICENSE file for details.
 */

#ifndef DEFORM_ARAP_H
#define DEFORM_ARAP_H

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/SVD>
#include <vector>
#include <unordered_map>
#include <iostream>

namespace deform {
    
    template<class T> class PrivateAccessor;
    
    template<class Mesh>
    class AsRigidAsPossibleDeformation {
    public:
        typedef typename Mesh::Scalar Scalar;
        typedef int Index;
        typedef typename Eigen::Matrix<Scalar, 3, 1> Point;
        
        AsRigidAsPossibleDeformation(Mesh &mesh)
            : _mesh(mesh), _dirty(true)
        {
            initializeMeshTopology();
        }
        
        void setConstraint(Index vidx, Eigen::Ref<const Point> loc) {
            _constrainedLocations[vidx] = loc;
            _dirty = true;
        }
        
        bool deform(Index numberOfIterations) {
            if (_dirty) {
                
                initializeMeshGeometry();
                computeCotanWeights();
                initializeRotations();
                initializeFreeVariableMapping();
                initializeConstraints();
                
                if (_numberOfFreeVariables == _mesh.numberOfVertices())
                    return true;
                
                if (!setupLinearSystem())
                    return false;
                
                _dirty = false;
            }
            
            for (Index i = 0; i < numberOfIterations; ++i) {
                
                // Estimate rotations
                estimateRotations();
                
                // Estimate positions
                estimatePositions();
            }
            
            // Update vertex locations at mesh
            
            for (Index i = 0; i < _mesh.numberOfVertices(); ++i) {
                _mesh.vertexLocation(i, _pprime.col(i));
            }

            return true;
        }
        
        const Eigen::Matrix<Scalar, 3, Eigen::Dynamic> &updatedPositions() const {
            return _pprime;
        }
        
    private:
        template<class> friend class PrivateAccessor;
        
        void initializeMeshTopology() {
            
            _faces.resize(3, _mesh.numberOfFaces());
            for (Index f = 0; f < _mesh.numberOfFaces(); ++f) {
                _faces.col(f) = _mesh.face(f);
            }
        }
        
        void initializeMeshGeometry() {
            _p.resize(3, _mesh.numberOfVertices());
            for (Index v = 0; v < _mesh.numberOfVertices(); ++v) {
                _p.col(v) = _mesh.vertexLocation(v);
            }
            _pprime = _p;
        }
        
        void computeCotanWeights() {
            // Based on https://igl.ethz.ch/projects/bbw/a-cotangent-laplacian-for-images-as-surfaces-2012-jacobson-sorkine.pdf
            
            typedef Eigen::Triplet<Scalar> T;
            std::vector<T> triplets;
            triplets.reserve(_faces.cols() * 3);
            
            // For each face
            for (Eigen::DenseIndex f = 0; f < _faces.cols(); ++f) {
                auto vids = _faces.col(f);
                
                auto v0 = _p.col(vids(0));
                auto v1 = _p.col(vids(1));
                auto v2 = _p.col(vids(2));
                
                
                // TODO: How to handle near degenerate triangles?
                // https://igl.ethz.ch/projects/deformation-survey/deformation_survey.pdf
                Scalar l0 = std::min<Scalar>(Scalar(1e-8), (v1 - v0).squaredNorm());
                Scalar l1 = std::min<Scalar>(Scalar(1e-8), (v2 - v1).squaredNorm());
                Scalar l2 = std::min<Scalar>(Scalar(1e-8), (v0 - v2).squaredNorm());
                
                l0 = std::sqrt(l0);
                l1 = std::sqrt(l1);
                l2 = std::sqrt(l2);
                
                const Scalar semip = Scalar(0.5) * (l0 + l1 + l2);
                const Scalar area = std::sqrt(semip*(semip - l0)*(semip - l1) * (semip - l2));
                
                const Scalar denom = Scalar(1.0) / (Scalar(4.0) * area);
                
                const Scalar cot0 = (-l0 * l0 + l1 * l1 + l2 * l2) * denom;
                const Scalar cot1 = (l0 * l0 - l1 * l1 + l2 * l2) * denom;
                const Scalar cot2 = (l0 * l0 + l1 * l1 - l2 * l2) * denom;
                
                std::pair<Index, Index> e0 = undirectedEdge(vids(0), vids(1));
                std::pair<Index, Index> e1 = undirectedEdge(vids(1), vids(2));
                std::pair<Index, Index> e2 = undirectedEdge(vids(2), vids(0));
                
                // Note that duplicate triplets are summed.
                triplets.push_back(T(e0.first, e0.second, cot0 * Scalar(0.5)));
                triplets.push_back(T(e1.first, e1.second, cot1 * Scalar(0.5)));
                triplets.push_back(T(e2.first, e2.second, cot2 * Scalar(0.5)));
                
                // We also fill the lower triangular region for convenience
                triplets.push_back(T(e0.second, e0.first, cot0 * Scalar(0.5)));
                triplets.push_back(T(e1.second, e1.first, cot1 * Scalar(0.5)));
                triplets.push_back(T(e2.second, e2.first, cot2 * Scalar(0.5)));
            }
            
            _edgeWeights.resize(_mesh.numberOfVertices(), _mesh.numberOfVertices());
            _edgeWeights.reserve(Eigen::VectorXi::Constant(_mesh.numberOfVertices(),7));
            _edgeWeights.setZero();
            _edgeWeights.setFromTriplets(triplets.begin(), triplets.end());
        }
        
        void initializeRotations() {
            _rotations.clear();
            _rotations.resize(_mesh.numberOfVertices(), RotationMatrix::Identity());
        }
        
        void initializeFreeVariableMapping() {
            
            _freeIdxMap.resize(_mesh.numberOfVertices());
            
            auto end = _constrainedLocations.end();
            Index freeIdx = 0;
            for (Index i = 0; i < _mesh.numberOfVertices(); ++i) {
                Index idx = _constrainedLocations.find(i) != end ? -1 : freeIdx++;
                _freeIdxMap[i] = idx;
            }
            _numberOfFreeVariables = freeIdx;
        }
        
        void initializeConstraints() {
            for (auto i = _constrainedLocations.begin(); i != _constrainedLocations.end(); ++i) {
                _pprime.col(i->first) = i->second;
            }
        }
        
        bool setupLinearSystem() {
            
            _L.resize(_numberOfFreeVariables, _numberOfFreeVariables);
            _L.reserve(Eigen::VectorXi::Constant(_numberOfFreeVariables, 7));
            _L.setZero();
            
            _bFixed.resize(3, _numberOfFreeVariables);
            _bFixed.setZero();
            
            _b.resize(3, _numberOfFreeVariables);
            
            typedef Eigen::Triplet<Scalar> T;
            std::vector<T> triplets;
            triplets.reserve(_numberOfFreeVariables * 7);
            
            // We use the implicitly stored topology of the mesh in the edgeWeight matrix
            // to determine one-ring relations of vertices.
            
            for (Index vi = 0; vi < _edgeWeights.outerSize(); ++vi) {
                
                Index iidx = _freeIdxMap[vi];
                
                if (iidx == -1) {
                    // Constrained, skip
                    continue;
                }
                
                for (typename SparseMatrix::InnerIterator vj(_edgeWeights, vi); vj; ++vj) {
                    
                    Index jidx = _freeIdxMap[vj.col()];
                    
                    const Scalar wij = vj.value();
                    
                    if (jidx == -1) {
                        // Neighbor is constrained, move to right side.
                        _bFixed.col(iidx) += wij * _constrainedLocations[vj.col()];
                    } else {
                        triplets.push_back(T(iidx, jidx, -wij));
                    }
                    
                    triplets.push_back(T(iidx, iidx, wij));
                }
            }
            
            _L.setFromTriplets(triplets.begin(), triplets.end());
            _solver.compute(_L);
            
            return _solver.info() == Eigen::Success;
        }
        
        void estimateRotations() {
            
            // Use the structure of the edge sparse matrix to traverse mesh topology.
            
            for (Index vi = 0; vi < _edgeWeights.outerSize(); ++vi) {
                
                const auto &pi = _p.col(vi);
                const auto &ppi = _pprime.col(vi);
                
                typedef Eigen::Matrix<Scalar, 3, 3> S;
                S cov = Eigen::Matrix<Scalar, 3, 3>::Zero();
                
                for (typename SparseMatrix::InnerIterator vj(_edgeWeights, vi); vj; ++vj) {
                    
                    const Scalar wij = vj.value();
                    
                    const auto &pj = _p.col(vj.col());
                    const auto &ppj = _pprime.col(vj.col());
                    
                    cov += wij * (pi - pj) * ((ppi - ppj).transpose());
                }
                
                Eigen::JacobiSVD<S> svd(cov, Eigen::ComputeFullU | Eigen::ComputeFullV);
                const S &v = svd.matrixV();
                const S ut = svd.matrixU().transpose();
                
                S id = S::Identity(3, 3);
                id(2, 2) = (v * ut).determinant();  // accounts for required flip
                _rotations[vi] = (v * id * ut);
            }
        }
        
        void estimatePositions() {
            
            _b = _bFixed;
            
            for (Index vi = 0; vi < _edgeWeights.outerSize(); ++vi) {
                
                Index iidx = _freeIdxMap[vi];
                
                if (iidx == -1) {
                    // Constrained, skip
                    continue;
                }
                
                for (typename SparseMatrix::InnerIterator vj(_edgeWeights, vi); vj; ++vj) {
                    
                    const Scalar wij = vj.value();
                    
                    RotationMatrix r = _rotations[vi] + _rotations[vj.col()];
                    Point pt = r * (_p.col(vi) - _p.col(vj.col())) * wij * Scalar(0.5);
                    _b.col(iidx) += pt;
                }
            }
            
            // Solve Lp' = b;
            
            Eigen::Matrix<Scalar, Eigen::Dynamic, 1> u;
            for (int d = 0; d < 3; ++d) {
                u = _solver.solve(_b.row(d).transpose());
                
                // Update pprime
                Index idx = 0;
                for (Index vi = 0; vi < _mesh.numberOfVertices(); ++vi) {
                    if (_freeIdxMap[vi] != -1) {
                        _pprime(d, vi) = u(idx++);
                    }
                }
            }
        }
        
        std::pair<Index, Index> undirectedEdge(Index a, Index b) {
            if (a > b)
                return std::pair<Index, Index>(b, a);
            else
                return std::pair<Index, Index>(a, b);
        }
                              
        typedef std::unordered_map<Index, Point> MapIndexToPoint;
        typedef Eigen::Matrix<Scalar, 3, 3> RotationMatrix;
        typedef Eigen::SparseMatrix<Scalar, Eigen::RowMajor> SparseMatrix;
        
        Mesh &_mesh;
        Eigen::Matrix<Scalar, 3, Eigen::Dynamic> _p, _pprime;
        Eigen::Matrix<Index, 3, Eigen::Dynamic> _faces;
        SparseMatrix _edgeWeights;
        std::vector< Eigen::Matrix<Scalar, 3, 3> > _rotations;
        
        Index _numberOfFreeVariables;
        std::vector<Index> _freeIdxMap;
        
        MapIndexToPoint _constrainedLocations;
        bool _dirty;
        
        
        typedef Eigen::SparseMatrix<Scalar> LMatrixType;
        LMatrixType _L;
        
        Eigen::Matrix<Scalar, 3, Eigen::Dynamic> _bFixed, _b;
        Eigen::SimplicialLDLT<LMatrixType> _solver;
        
    };
}

#endif