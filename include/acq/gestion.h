#ifndef ACQ_GESTION_H
#define ACQ_GESTION_H

#include "acq/cloudPrimitive.h"
#include "acq/cloudManager.h"
#include "acq/typedefs.h"

#include <Eigen/Dense>

namespace acq {

    //  ****** ============ Helper Functions  =============== ******* 
    Eigen::MatrixXi sample(int cloudSize) ; // sample the point cloud
    double computeVariance(Eigen::MatrixXd V); // compute the variance 

    //  ****** ============ Functions to handle a sphere =============== ******* 
    // determine if the 3 points gives a sphere or no 
    bool isSphere(Eigen::Matrix3d vertices, Eigen::Matrix3d normals, double threshold, double alpha) ;
    
    // create and store the sphere when it exists 
    bool computeSphere(Eigen::Matrix<int, 3,1> sample_idx, Eigen::Matrix3d variance, DecoratedCloud& cloud, CloudPrimitive& primitives, double threshold, double alpha); 

    // compute the attributs using as many points as we have 
    Eigen::Matrix<double, 1,3> computerCenter(Eigen::MatrixXd vertices, Eigen::MatrixXd normals) ;
    double computerRadius(Eigen::MatrixXd thisVertices, Eigen::Matrix<double, 1,3> thisCenter) ;

    /********* ============= Functions to handle PLANE =============== *********/
    bool computePlane(Eigen::Matrix3i sample_idx, Eigen::Matrix3d variance, DecoratedCloud &cloud, CloudPrimitive& primitives, double T);
    bool isPlane(Eigen::MatrixXd V, Eigen::MatrixXd N, double T, double alpha);
    Eigen::RowVector3d computeNormal(Eigen::MatrixXd V, Eigen::MatrixXd _N);

    /********* ============= Functions to handle the final cloud =============== *********/
     
    DecoratedCloud fuse(CloudPrimitive& best_primitives, CloudManager& clouds) ;
    bool cleanCloud(DecoratedCloud& oldCloud, DecoratedCloud& newCloud, Eigen::Matrix3i inliers_idx) ;
}

#endif

