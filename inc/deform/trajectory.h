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
#include <unsupported/Eigen/Splines>

#ifdef _WIN32
#pragma warning (push)
#pragma warning (disable : 4244)
#endif

#include <sophus/sophus.hpp>
#include <sophus/se3.hpp>

#ifdef _WIN32
#pragma warning (pop)
#endif

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
        TrajectorySE3()
            :_dirty(false)
        {}

        Transform addKeyPose(const Transform &transform) {
            _transforms.push_back(transform);
            _dirty = true;
            return transform;
        }

        Transform operator()(Scalar time) {
            if (_dirty) {
                Eigen::Matrix<Scalar, 6, Eigen::Dynamic> points(6, _transforms.size());
                for (size_t i = 0; i < _transforms.size(); ++i) {
                    points.col(i) = SE3Group(_transforms[i].matrix()).log();
                }
                _spline = Eigen::SplineFitting<SplineType>::Interpolate(points, 3);
                _dirty = false;
             }

            typename SE3Group::Tangent tangent = _spline(time);
            return SE3Group::exp(tangent).affine3();
        }
        
    private:
        typedef Eigen::Spline<Scalar, 6, 3> SplineType;
        typedef std::vector<Transform, Eigen::aligned_allocator<Transform> > TransformVector;
        
        bool _dirty;
        SplineType _spline;

        TransformVector _transforms;
    };
}

#endif