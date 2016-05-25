/**
 This file is part of mesh-deform.
 
 Copyright(C) 2016 Christoph Heindl
 All rights reserved.
 
 This software may be modified and distributed under the terms
 of the BSD license.See the LICENSE file for details.
 */

#ifndef DEFORM_DEFORMATION_PATH_H
#define DEFORM_DEFORMATION_PATH_H

#include <deform/trajectory.h>
#include <deform/arap.h>

namespace deform {
    
    template<class MeshType>
    class DeformationUtil {
    public:

        EIGEN_MAKE_ALIGNED_OPERATOR_NEW;

        typedef MeshType Mesh;

        /** Floating point precision used in calculations. */
        typedef typename MeshType::Scalar Scalar;
              
        /** Transformation matrix type. */
        typedef Eigen::Transform<Scalar, 3, Eigen::Affine> Transform;
        

        template<class HandleIterator>
        DeformationUtil(const Mesh &mesh, HandleIterator handlesBegin, HandleIterator handlesEnd, const Transform &origin = Transform::Identity())
            :_origin(origin)
        {
            _originInv = origin.inverse(Eigen::Isometry);
            _handles.insert(_handles.begin(), handlesBegin, handlesEnd);
            _points.resize(3, _handles.size());
            
            for (size_t i = 0; i < _handles.size(); ++i) {
                _points.col(i) = mesh.vertexLocation(_handles[i]);
            }

        }

        template<class ARAP>
        void updateConstraints(const Transform &t, ARAP &arap)
        {
            Transform tabs = _origin * t * _originInv;

            for (size_t i = 0; i < _handles.size(); ++i) {
                arap.setConstraint(_handles[i], tabs * _points.col(i));
            }
        
        }
        
    private:
        std::vector<int> _handles;
        Eigen::Matrix<Scalar, 3, Eigen::Dynamic> _points;
        Transform _origin;
        Transform _originInv;
    };
}

#endif