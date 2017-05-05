#ifndef ACQ_GESTION_H
#define ACQ_GESTION_H

#include "acq/cloudPrimitive.h"
#include "acq/cloudManager.h"
#include "acq/typedefs.h"

#include <Eigen/Dense>
#include <iostream>

namespace acq {

    //  ****** ============ Helper Functions  =============== ******* 
    Eigen::MatrixXi sample(int cloudSize) ; // sample the point cloud
    Eigen::Matrix3d computeVariance(Eigen::MatrixXd V); // compute the variance 

    //  ****** ============ Functions to handle a sphere =============== ******* 
    // determine if the 3 points gives a sphere or no 
    bool isSphere(Eigen::Matrix3d vertices, Eigen::Matrix3d normals, double threshold, double alpha) ;
    
    // create and store the sphere when it exists 
    void computeSphere(Eigen::Matrix<int, 3,1> sample_idx, Eigen::Matrix3d variance, DecoratedCloud& cloud, CloudPrimitive& primitives, double threshold, double alpha); 

    // compute the attributs using as many points as we have 
    Eigen::Matrix<double, 1,3> computerCenter(Eigen::MatrixXd vertices, Eigen::MatrixXd normals) ;
    double computerRadius(Eigen::MatrixXd thisVertices, Eigen::Matrix<double, 1,3> thisCenter) ;

    /********* ============= Functions to handle PLANE =============== *********/
    void computePlane(Eigen::Matrix<int, 3,1> sample_idx, Eigen::Matrix3d variance, DecoratedCloud &cloud, CloudPrimitive& primitives, double thresh, double alpha);
    bool isPlane(Eigen::MatrixXd V, Eigen::MatrixXd N, double T, double alpha);
    Eigen::Matrix<double, 1,3> computeNormal(Eigen::MatrixXd V, Eigen::Matrix<double, 1,3> _N);

    /********* ============= Functions to handle the final cloud =============== *********/
    void fuse(CloudPrimitive& best_primitives, CloudManager& clouds, double T_rad, double T_cent, double T_norm, double T_refPt);
    DecoratedCloud& gatherClouds(CloudManager& cloudManager) ;
    void cleanCloud(DecoratedCloud& cloudRansac, CloudManager& cloudManager, Eigen::Matrix3i inliers_idx) ;
}

#endif

