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

#include "osg_viewer.h"
#include "example_config.h"

#include <Eigen/Geometry>

#include <string>
#include <iostream>


const std::vector<int> handles = { 8,9,10,11,15,18,19,20,21,35,38,39,40,41,51,54,55,56,57,74,75,76,77,78,79,103,104,107,108,109,110,111,112,113,114,116,123,126,127,128,129,135,136,141,142,144,145,146,147,148,149,150,151,162,163,179,182,183,184,185,195,198,199,200,201,203,205,218,219,220,221,222,223,261,262,296,297,298,299,300,301,302,303,304,305,306,307,308,309,310,311,312,313,315,316,324,325,342,343,351,352,357,358,360,361,362,363,364,365,366,367,378,379 };


const std::vector<int> anchors = {0,1,2,3,4,5,6,7,12,26,27,28,29,30,33,34,46,47,48,92,93,94,95,96,97,98,99,100,101,102,105,106,115,117,118,119,120,137,152,155,158,164,165,166,167,168,169,170,171,172,173,174,177,178,190,191,192,193,194,196,197,202,204,211,226,229,232,236,237,238,239,240,241,242,243,244,245,246,247,248,249,250,251,252,253,254,255,256,257,258,259,269,271,275,277,278,280,284,286,287,289,293,295,332,334,338,340,380,381,382,383,384,385 };


typedef OpenMesh::TriMesh_ArrayKernelT<> Mesh;
typedef deform::OpenMeshAdapter<> Adapter;
typedef deform::AsRigidAsPossibleDeformation<Adapter, double> ARAP;

/** 
    Interpolate two rotations based on parameter in interval [0, 1] 
**/
Eigen::Matrix3f slerpRotation(Eigen::Quaternionf start, Eigen::Quaternionf end, float t) {
    return start.slerp(t, end).matrix();
}

void applyRotationToHandles(const Eigen::Matrix3f &incrementalRotation, Mesh &mesh, ARAP &arap) {
    for (auto h : handles) {
        Mesh::VertexHandle vh = mesh.vertex_handle(h);
        Eigen::Vector3f v = incrementalRotation * deform::convert::toEigen(mesh.point(vh));
        arap.setConstraint(vh.idx(), v);
    }
}


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
    
    // Motion related variables
    Eigen::Matrix3f prev = Eigen::Matrix3f::Identity();
    Eigen::Quaternionf start(Eigen::AngleAxisf(0, Eigen::Vector3f::UnitX()));
    Eigen::Quaternionf end(Eigen::AngleAxisf(M_PI, Eigen::Vector3f::UnitX()));
    
    // Create an OSG viewer to visualize incremental deformation.
    float inc = 0.001;
    float t = 0.f;
    deform::example::OSGViewer viewer(mesh, anchors, handles);
    viewer.onFrame([&](Mesh &mesh, double time) {
        
        t+= inc;
        if (t > 1.f)
            return false;

        Eigen::Matrix3f r = slerpRotation(start, end, t);
        applyRotationToHandles(r * prev.transpose(), mesh, arap);
        prev = r;
        
        arap.deform(20);
        
        return true;
    });
    viewer.run();
    
    return 0;

}
