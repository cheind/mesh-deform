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
    
    /**
        Forward declaration of a class that has access to internals of AsRigidAsPossibleDeformation.
        Used during tests.
    */
    template<class T> class PrivateAccessor;
    
    /**
        As-rigid-as-possible (ARAP) surface deformation.
     
        Based on
        Sorkine, Olga, and Marc Alexa. "As-rigid-as-possible surface modeling." 
        Symposium on Geometry processing. Vol. 4. 2007.
     
        This class is an implementation of the algorithm outlined in paper "As-rigid-as-possible surface modeling". It
        takes a triangular mesh and a sparse set of anchored target vertex positions. The algorithm then optimizes
        the remaining vertex positions in such a way that local rigidity is conserved and target positions are reached.
     
        The implementation is mesh surface representation agnostic. In order to use it with a custom mesh type, you need to specify 
        a mesh (adapter) that provides a minimalistic set of mesh accessors. See examples for such adapters.
     
        \note By default the precision used in floating point computations is the same as the precision of the vertex locations.
              In case you observe simulation artefacts after running the simulation for some time, you might want to instruct
              AsRigidAsPossibleDeformation to use higher precision floating point type using the PrecisionType
              template argument.
    */
    template<class MeshType, class PrecisionType = typename MeshType::Scalar>
    class AsRigidAsPossibleDeformation {
    public:
        /** Floating point precision used in calculations. */
        typedef PrecisionType Scalar;
        
        /** Mesh type we are working on. */
        typedef MeshType Mesh;
        
        /** Type for indices */
        typedef int Index;
        
        /**
            Construct from mesh. 
            
            \param mesh Mutable reference to mesh. Mesh needs to outlive this object.
        */
        AsRigidAsPossibleDeformation(Mesh &mesh)
            : _mesh(mesh), _dirty(true)
        {
            initializeMeshTopology();
        }
        
        /**
            Set positional constrains for a specific vertex.
         
            Constraining a vertex to a specific location provides boundary conditions
            for the optimization.
            
            \param vidx Index of vertex to constrain.
            \param loc Location to pin vertex at.
        */
        template<class S>
        void setConstraint(Index vidx, const Eigen::Matrix<S, 3, 1 > &loc) {
            _constrainedLocations[vidx] = loc.template cast<Scalar>();
            _dirty = true;
        }
        
        /**
            Run the deformation optimization for a given number of iterations.
         
            The initial guess for the optimization is provided by the current vertex positions. 
            Since the objective function is not convex, the constrained locations
            should not be too far apart from this initial guess. Usually 2-5 iterations suffice 
            when the initial guess is good.
         
            The mesh topology is not supposed to change during successive calls of deform.
         
            \param numberOfIterations Number of iterations to perform.
            \returns Indication of success.
        */
         
        bool deform(Index numberOfIterations) {
            if (_dirty) {
                
                // Do all the initialization that only needs to be done
                // once for the mesh.
                
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
            
            // Estimate positions
            estimatePositions();
            
            for (Index i = 0; i < numberOfIterations; ++i) {
                
                // Estimate rotations
                estimateRotations();
                
                // Estimate positions
                estimatePositions();
            }
            
            // Back propagate updated vertex positions.
            
            for (Index i = 0; i < _mesh.numberOfVertices(); ++i) {
                _mesh.vertexLocation(i, _pprime.col(i).template cast<typename Mesh::Scalar>());
            }

            return true;
        }
        
    private:
        template<class> friend class PrivateAccessor;
        
        
        /**
            Initialize the mesh topology.
            
            Done at construction of this object.
        */
        void initializeMeshTopology() {
            
            _faces.resize(3, _mesh.numberOfFaces());
            for (Index f = 0; f < _mesh.numberOfFaces(); ++f) {
                _faces.col(f) = _mesh.face(f);
            }
        }
        
        /**
            Initialize mesh geometry
         
            Copies vertices and initializes target positions.
        */
        void initializeMeshGeometry() {
            _p.resize(3, _mesh.numberOfVertices());
            for (Index v = 0; v < _mesh.numberOfVertices(); ++v) {
                _p.col(v) = _mesh.vertexLocation(v).template cast<Scalar>();
            }
            _pprime = _p;
        }
        
        /**
            Compute cotan edge weights.
         
            Formulae taken from
            Jacobson, Alec, and Olga Sorkine. "A cotangent Laplacian for images as surfaces." 
            ACM Trans. Graph 25.3 (2012): 646-653.
         
            The formulae are robustified by handling edge cases such as degenerate triangles. 
            Note that the sparse matrix created by this method also defines the edge connections
            of vertices. This fact is used elsewhere in this algorithm to iterate over neighboring
            vertices.
        */
        void computeCotanWeights() {
            
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
                Scalar l0 = std::max<Scalar>(Scalar(1e-8), (v1 - v0).squaredNorm());
                Scalar l1 = std::max<Scalar>(Scalar(1e-8), (v2 - v1).squaredNorm());
                Scalar l2 = std::max<Scalar>(Scalar(1e-8), (v0 - v2).squaredNorm());
                
                l0 = std::sqrt(l0);
                l1 = std::sqrt(l1);
                l2 = std::sqrt(l2);
                
                const Scalar semip = Scalar(0.5) * (l0 + l1 + l2);
                const Scalar area = std::max<Scalar>(Scalar(1e-8), std::sqrt(semip*(semip - l0)*(semip - l1) * (semip - l2)));
                
                const Scalar denom = Scalar(1.0) / (Scalar(4.0) * area);
                
                Scalar cot0 = (-l0 * l0 + l1 * l1 + l2 * l2) * denom;
                Scalar cot1 = (l0 * l0 - l1 * l1 + l2 * l2) * denom;
                Scalar cot2 = (l0 * l0 + l1 * l1 - l2 * l2) * denom;
                
                cot0 = std::max<Scalar>(Scalar(1e-8), cot0);
                cot1 = std::max<Scalar>(Scalar(1e-8), cot1);
                cot2 = std::max<Scalar>(Scalar(1e-8), cot2);
                
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
        
        /**
            Initializes the rotations per vertex.
         
            All rotations are set to identity.
        */
        void initializeRotations() {
            _rotations.clear();
            _rotations.resize(_mesh.numberOfVertices(), RotationMatrix::Identity());
        }
        
        /**
            Initialize the free variable mapping.
         
            In order to use Cholesky decomposition we need to ensure a symmetric matrix.
            For the sparse matrix L we will only generate a row per free (unconstrained)
            vertex. All constrained vertices are moved to the right hand side of the 
            linear system of equations Lp'=b
         
            This method builds a map from vertex index to free variable index.
        */
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
        
        /**
            Initialize geometric constraints.
        */
        void initializeConstraints() {
            for (auto i = _constrainedLocations.begin(); i != _constrainedLocations.end(); ++i) {
                _pprime.col(i->first) = i->second;
            }
        }
        
        /**
            Setup the linear system of equations.
         
            As argued in the paper, the matrix L can be prefactored for every call to deform. 
            This methods builds L only over the unconstrained vertices and moves all references
            to constrained vertex positions (i.e known variables) to the right hand side.
         
            \returns Indication of success.
        */
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
        
        /**
            Estimate best matching rotations for each cell.
         
            This method estimates the best matching rotation matrix between a cell made of
            original vertex positions and a cell consisting of optimized positions.
         
            Rotation estimation for known correspondences is a closed form step based on SVD
            decomposition of a weighted covariance matrix. 
         
            Based on the flip-flop optimization scheme presented in the paper, this step
            assumes the optimized positions to be known and solves only for the rotations.
        */
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
        
        /**
            Optimize vertex locations.
         
            Assumes that the rotations are known and solves for the unknown positions.
            Note that the right hand side is not initialized to zero but to a constant 
            vector stemming from erased constrained variables from the left hand side.
        */
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
        
        /**
            Sort an edge, turning making a pair of vertices into an undirected edge.
        */
        std::pair<Index, Index> undirectedEdge(Index a, Index b) {
            if (a > b)
                return std::pair<Index, Index>(b, a);
            else
                return std::pair<Index, Index>(a, b);
        }
        
        typedef Eigen::Matrix<Scalar, 3, 1> Point;
        typedef std::unordered_map<Index, Point> MapIndexToPoint;
        typedef Eigen::Matrix<Scalar, 3, 3> RotationMatrix;
        typedef Eigen::SparseMatrix<Scalar, Eigen::RowMajor> SparseMatrix;
        typedef Eigen::Matrix<Scalar, 3, Eigen::Dynamic> PointMatrix;
        typedef Eigen::Matrix<Index, 3, Eigen::Dynamic> FaceMatrix;
        typedef Eigen::SparseMatrix<Scalar> LSolverMatrix;
        
        Mesh &_mesh;
        PointMatrix _p, _pprime;
        FaceMatrix _faces;
        SparseMatrix _edgeWeights;
        std::vector<RotationMatrix> _rotations;
        
        Index _numberOfFreeVariables;
        std::vector<Index> _freeIdxMap;
        MapIndexToPoint _constrainedLocations;
        
        LSolverMatrix _L;
        PointMatrix _bFixed, _b;
        Eigen::SimplicialLDLT<LSolverMatrix> _solver;
        
        bool _dirty;
        
    };
}

#endif