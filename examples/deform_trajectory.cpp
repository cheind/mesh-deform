/**
 This file is part of mesh-deform.
 
 Copyright(C) 2016 Christoph Heindl
 All rights reserved.
 
 This software may be modified and distributed under the terms
 of the BSD license.See the LICENSE file for details.

 This examples shows the deformation of a sphere like object. In this particular
 case a single vertex is pinned to its initial location and a harmonic motion
 is applied to the vertex on the opposite side of the sphere. All other vertices
 are unconstrained.
 */

#define _USE_MATH_DEFINES

#include <deform/arap.h>
#include <deform/openmesh_adapter.h>
#include <deform/trajectory.h>
#include <deform/deformation_util.h>

#include "osg_viewer.h"
#include "example_config.h"

#include <Eigen/Geometry>

#include <string>
#include <iostream>


const std::vector<int> handles = { 4,5,7,11,14,15,16,17,18,19,26,27,38,39,45,46,48,49,50,51,52,53,54,55,56,76,77,131,132,201,202,203,204,205,206,207,208,209,210,211,212,213,214,215,216,217,218,219,220,221,222,223,224,225,226,227,228,229,230,231,232,233,234,235,240,241,257,307,308,324,325,327,328,329,330,331,332,333,334,354,355 };
const std::vector<int> anchors = {0,1,2,3,79,80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,96,97,98,99,100,101,102,103,104,105,106,107,108,109,110,111,112,113,114,115,116,117,118,119,120,121,122,123,124,125,126,357,358,359,360,361,362,363,364,365,366,367,368,369,370,371,372,373,374,375,376,377,378,379,380,381,382,383,384,385};


typedef OpenMesh::TriMesh_ArrayKernelT<> Mesh;
typedef deform::OpenMeshAdapter<> Adapter;
typedef deform::AsRigidAsPossibleDeformation<Adapter, double> ARAP;
typedef deform::TrajectorySE3<float> Trajectory;
typedef deform::DeformationUtil<Adapter> DeformationUtil;


int main(int argc, char **argv) {
    
    std::string pathToSource = std::string(DEFORM_ETC_DIR) + std::string("/bar.obj");
    
    Mesh mesh;
    if (!OpenMesh::IO::read_mesh(mesh, pathToSource)) {
        std::cerr << "Failed to read mesh" << std::endl;
        return -1;
    }
    
    Adapter ma(mesh);
    ARAP arap(ma);
    
    // Set anchors.
    for (auto a : anchors) {
        Mesh::VertexHandle va = mesh.vertex_handle(a);
        arap.setConstraint(va.idx(), deform::convert::toEigen(mesh.point(va)));
    }
    
    // Setup trajectory
    Trajectory trajectory;
    
    Trajectory::Transform prev;
    prev = trajectory.addKeyPose(Trajectory::Transform::Identity());
    prev = trajectory.addKeyPose(prev * Eigen::Translation3f(1.f, 0.f, 0.f));
    prev = trajectory.addKeyPose(prev * Eigen::Translation3f(2.f, 0.f, 0.f));
    prev = trajectory.addKeyPose(prev * Eigen::AngleAxisf((float)M_PI / 2.f, Eigen::Vector3f::UnitX()));   
    
    DeformationUtil dutil(ma, handles.begin(), handles.end(), trajectory(0.0));
    
    // Create an OSG viewer to visualize incremental deformation.
    double duration = 10.f;

    deform::example::OSGViewer viewer(mesh, anchors, handles);
    viewer.onFrame([&](Mesh &mesh, double time) {
        

        double t = time / duration;
        if (t > 1.f)
            return false;
        
        dutil.updateConstraints(trajectory((float)t), arap);
        
        arap.deform(20);
        
        return true;
    });
    viewer.run();
    
    return 0;

}
