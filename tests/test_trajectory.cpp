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


    Trajectory path;

    // time 0
    Trajectory::Transform pose0 = Trajectory::Transform::Identity();
    path.addKeyPose(pose0);

    Trajectory::Transform pose1 = Eigen::Translation3f(0, 0, 1) * pose0;
    path.addKeyPose(pose1);

    Trajectory::Transform pose2 = Eigen::Translation3f(0, 0, 1) * Eigen::AngleAxisf((float)M_PI / 4.f, Eigen::Vector3f::UnitX()) * pose1;
    path.addKeyPose(pose2);

    Trajectory::Transform pose3 = Eigen::Translation3f(0, 0, 1) * pose2;
    path.addKeyPose(pose3);

    REQUIRE(path(0.f).matrix().isApprox(pose0.matrix(), 1e-3f));
    REQUIRE(path(1.f).matrix().isApprox(pose3.matrix(), 1e-3f));
}
    