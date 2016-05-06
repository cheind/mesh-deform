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


const std::vector<int> handles = { 12,13,14,15,40,41,42,43,66,94,95,96,97,120,148,149,150,151,174,203,204,205,257,258,259,310,311,312,313,336,364,365,366,367,390,419,420,421,473,474,475,526,527,528,529,552,580,581,582,583,606,635,636,637,689,690,691,812,813,814,815,816,817,818,819,820,974,975,976,977,978,979,980,981,982,1067,1068,1069,1121,1122,1123,1174,1175,1176,1177,1200,1228,1229,1230,1231,1254,1283,1284,1285,1337,1338,1339,1390,1391,1392,1393,1416,1444,1445,1446,1447,1470,1499,1500,1501,1553,1554,1555,1676,1677,1678,1679,1680,1681,1682,1683,1684,1838,1839,1840,1841,1842,1843,1844,1845,1846,1931,1932,1933,1985,1986,1987,2038,2039,2040,2041,2064,2092,2093,2094,2095,2118,2147,2148,2149,2201,2202,2203,2254,2255,2256,2257,2280,2308,2309,2310,2311,2334,2363,2364,2365,2417,2418,2419,2540,2541,2542,2543,2544,2545,2546,2547,2548,2702,2703,2704,2705,2706,2707,2708,2709,2710,2795,2796,2797,2849,2850,2851,2972,2973,2974,2975,2976,2977,2978,2979,2980,3134,3135,3136,3137,3138,3139,3140,3141,3142,3395,3396,3397,3398,3399,3400,3401,3402,3403,3404,3405,3406,3407,3408,3409,3410,3411,3412,3413,3414,3415,3416,3417,3418,3419,3420,3421,3881,3882,3883,3884,3885,3886,3887,3888,3889,3890,3891,3892,3893,3894,3895,3896,3897,3898,3899,3900,3901,3902,3903,3904,3905,3906,3907,4268,4269,4270,4271,4272,4273,4274,4275,4276,4430,4431,4432,4433,4434,4435,4436,4437,4438,4523,4524,4525,4577,4578,4579 };


const std::vector<int> anchors = {16,17,18,19,48,49,50,51,70,102,103,104,105,124,156,157,158,159,178,191,192,193,245,246,247,318,319,320,321,340,372,373,374,375,394,407,408,409,461,462,463,534,535,536,537,556,588,589,590,591,610,623,624,625,677,678,679,848,849,850,851,852,853,854,855,856,1010,1011,1012,1013,1014,1015,1016,1017,1018,1055,1056,1057,1109,1110,1111,1182,1183,1184,1185,1204,1236,1237,1238,1239,1258,1271,1272,1273,1325,1326,1327,1398,1399,1400,1401,1420,1452,1453,1454,1455,1474,1487,1488,1489,1541,1542,1543,1712,1713,1714,1715,1716,1717,1718,1719,1720,1874,1875,1876,1877,1878,1879,1880,1881,1882,1919,1920,1921,1973,1974,1975,2046,2047,2048,2049,2068,2100,2101,2102,2103,2122,2135,2136,2137,2189,2190,2191,2262,2263,2264,2265,2284,2316,2317,2318,2319,2338,2351,2352,2353,2405,2406,2407,2576,2577,2578,2579,2580,2581,2582,2583,2584,2738,2739,2740,2741,2742,2743,2744,2745,2746,2783,2784,2785,2837,2838,2839,3008,3009,3010,3011,3012,3013,3014,3015,3016,3170,3171,3172,3173,3174,3175,3176,3177,3178,3287,3288,3289,3290,3291,3292,3293,3294,3295,3296,3297,3298,3299,3300,3301,3302,3303,3304,3305,3306,3307,3308,3309,3310,3311,3312,3313,3773,3774,3775,3776,3777,3778,3779,3780,3781,3782,3783,3784,3785,3786,3787,3788,3789,3790,3791,3792,3793,3794,3795,3796,3797,3798,3799,4304,4305,4306,4307,4308,4309,4310,4311,4312,4466,4467,4468,4469,4470,4471,4472,4473,4474,4511,4512,4513,4565,4566,4567 };


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
    
    std::string pathToSource = std::string(DEFORM_ETC_DIR) + std::string("/bar2.obj");
    
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
    float inc = 0.01;
    float t = 0.f;
    deform::example::OSGViewer viewer(argc, argv, mesh, [&](Mesh &mesh, double time) {
        
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
