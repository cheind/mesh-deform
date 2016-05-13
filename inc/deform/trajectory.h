/**
 This file is part of mesh-deform.
 
 Copyright(C) 2016 Christoph Heindl
 All rights reserved.
 
 This software may be modified and distributed under the terms
 of the BSD license.See the LICENSE file for details.
 */

#ifndef DEFORM_TRAJECTORY_H
#define DEFORM_TRAJECTORY_H

#include <Eigen/Dense>
#include <Eigen/Geometry>
#include <Eigen/StdVector>

#ifdef _WIN32
#pragma warning (push)
#pragma warning (disable : 4244)
#endif

#include <sophus/sophus.hpp>
#include <sophus/se3.hpp>

#ifdef _WIN32
#pragma warning (pop)
#endif

#include <iostream>

namespace deform {
    
    template<class PrecisionType>
    class TrajectorySE3 {
    public:

        EIGEN_MAKE_ALIGNED_OPERATOR_NEW;

        /** Floating point precision used in calculations. */
        typedef PrecisionType Scalar;

        /** Group type */
        typedef Sophus::SE3Group<Scalar> SE3Group;
               
        /** Transformation matrix type. */
        typedef Eigen::Transform<Scalar, 3, Eigen::Affine> Transform;
        
        /**
            Construct empty trajectory
        */
        TrajectorySE3(Scalar timeStep)
            :_timeStep(timeStep)
        {
            _C << 
                6, 0, 0, 0,
                5, 3, -3, 1,
                1, 3, 3, -2,
                0, 0, 0, 1;
            _C *= Scalar(1) / 6;
        }

        void addKeyframe(const Transform &transform) {
            _keyframes.push_back(transform);
            const size_t n = _keyframes.size();
            
            if (n > 1) {
                Transform delta = _keyframes[n - 2].inverse(Eigen::Isometry) * _keyframes[n - 1];
               
                SE3Group g(delta.matrix());
                _omega.push_back(g.log());

                if (n == 2) {
                    _omega[0] = g.inverse().log();
                }
            } else {
                _omega.push_back(SE3Group::Tangent::Zero()); // Will be fixed when second keyframe is given.
            }
        }

        Transform operator()(Scalar time) const {
    
            const Scalar t0(0);
            const Scalar s = (time - t0) / _timeStep;
            
            const int i = static_cast<int>(std::floor(s));
            eigen_assert(i >= 0);

            const Scalar u = s - Scalar(i);
            const Eigen::Matrix<Scalar, 4, 1> b = _C * Eigen::Matrix<Scalar, 4, 1>(Scalar(1), u, u*u, u*u*u);
            
            Transform trans = previousTransform(i);
            
            Transform prod = Transform::Identity();
            for (int j = 1; j < 4; ++j) {
                typename SE3Group::Tangent o = omega(i + j);
                prod = prod * SE3Group::exp(b(j) * o).affine3();
            }
            
            return trans * prod;
        }
        
    private:

        typename SE3Group::Tangent omega(int idx) const {
            if ((size_t)idx < _omega.size()) {
                return _omega[idx];
            } else {
                // Extrapolate
                return _omega.back();
            }
        }

        Transform previousTransform(int idx) const {
            if (idx == 0) {
                // Extrapolate
                return SE3Group::exp(_omega[0]).affine3() * _keyframes[0];
            } else {
                return _keyframes[idx - 1];
            }
        }

        Scalar _timeStep;
        std::vector<Transform, Eigen::aligned_allocator<Transform> > _keyframes;
        std::vector<typename SE3Group::Tangent, Eigen::aligned_allocator<typename SE3Group::Tangent> > _omega;
        Eigen::Matrix<Scalar, 4, 4> _C;
    };
}

#endif