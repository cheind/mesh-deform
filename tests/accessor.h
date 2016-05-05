/**
 This file is part of mesh-deform.
 
 Copyright(C) 2016 Christoph Heindl
 All rights reserved.
 
 This software may be modified and distributed under the terms
 of the BSD license.See the LICENSE file for details.
*/

#ifndef DEFORM_ACCESSOR_H
#define DEFORM_ACCESSOR_H

namespace deform {
    
    template<class ARAP>
    class PrivateAccessor {
    public:
        static typename ARAP::SparseMatrix cotanWeights(const ARAP &a) {
            return a._edgeWeights;
        }
    };
}


#endif
