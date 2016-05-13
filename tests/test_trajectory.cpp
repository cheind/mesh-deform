/**
 This file is part of mesh-deform.
 
 Copyright(C) 2016 Christoph Heindl
 All rights reserved.
 
 This software may be modified and distributed under the terms
 of the BSD license.See the LICENSE file for details.
*/

#include "catch.hpp"

#include <deform/trajectory.h>
#include <iostream>


TEST_CASE("trajectory")
{
    typedef deform::TrajectorySE3<float> Trajectory;


    Trajectory path(1.f);

    // time 0
    Trajectory::Transform trans = Trajectory::Transform::Identity();
    path.addKeyframe(trans);

    // time 1
    trans *= Eigen::Translation3f(0, 0, 1);
    path.addKeyframe(trans);

    // time 2
    trans *= Eigen::Translation3f(0, 0, 1);
    path.addKeyframe(trans);

    // time 3
    trans *= Eigen::Translation3f(0, 0, 1);
    path.addKeyframe(trans);

    // time 4
    trans = trans * Eigen::Translation3f(0, 0, 1) * Eigen::AngleAxisf(M_PI / 4, Eigen::Vector3f::UnitX());
    std::cout << "Expected transform" << std::endl << trans.matrix() << std::endl;
    path.addKeyframe(trans);

    // time 5
    trans *= Eigen::Translation3f(0, 0, 1);
    path.addKeyframe(trans);

    // time 6
    trans *= Eigen::Translation3f(0, 0, 1);
    path.addKeyframe(trans);

    // time 7
    trans *= Eigen::Translation3f(0, 0, 1);
    path.addKeyframe(trans);

    // time 8
    trans *= Eigen::Translation3f(0, 0, 1);
    path.addKeyframe(trans);


    std::cout << "Transform is" << std::endl;
    std::cout << path(4.f).matrix() << std::endl;
}
    