#ifndef ACQ_GESTION_H
#define ACQ_GESTION_H

#include "acq/cloudPrimitive.h"
#include "acq/typedefs.h"

#include <Eigen/Dense>

namespace acq {

    // ****** Sample the point cloud and compute its variance *******
    Eigen::MatrixXi sample(int cloudSize) ; // Renvoie matrice d'index de vertices
    double computeVariance(Eigen::MatrixXd V);

    //  ****** Functions to handle a sphere ******* 
    // determine if the 3 points gives a sphere or no 
    bool isSphere(Eigen::Matrix3d vertices, Eigen::Matrix3d normals, double threshold, double alpha) ;
    
    // create and store the sphere when it exists 
    bool computeSphere(Eigen::Matrix<int, 3,1> sample_idx, Eigen::Matrix3d variance, DecoratedCloud& cloud, CloudPrimitive& primitives, double threshold, double alpha); 

    // compute the attributs using as many points as we have 
    Eigen::Matrix<double, 1,3> computerCenter(Eigen::MatrixXd vertices, Eigen::MatrixXd normals) ;
    double computerRadius(Eigen::MatrixXd thisVertices, Eigen::Matrix<double, 1,3> thisCenter) ;

    //  ****** Functiosn to handle a plane ******* 
    bool computePlane(Eigen::Matrix3i sample_idx, Eigen::Matrix3d variance, DecoratedCloud &cloud, CloudPrimitive& primitives, double T);
    bool isPlane(Eigen::MatrixXd V, Eigen::MatrixXd N, double T, double alpha);
    Eigen::RowVector3d computeNormal(Eigen::MatrixXd V, Eigen::MatrixXd _N);

// pour la matrix de color a la fin :

// pour la matrix de color a la fin :
   // DecoratedCloud removeCloud(DecoratedCloud cloud, )
}

#endif

