#ifndef ACQ_GESTION_H
#define ACQ_GESTION_H

#include "acq/typedefs.h"
#include "decoratedCloud.h"
#include <Eigen/Dense>

namespace acq {

    Eigen::MatrixXi sample(int cloudSize, int numberPoint); // Renvoie matrice d'index de vertices
    double computeVariance(Eigen::MatrixXd V);

    // try to compute a sphere with this sample points 
    bool computeSphere(Eigen::Matrix3i sample_idx, Eigen::Matrix3d variance, DecoratedCloud& cloud, cloudPrimitive primitives, double threshold, double alpha); 
    bool computePlane(Eigen::Matrix3i sample_idx, Eigen::Matrix3d variance, DecoratedCloud& cloud, cloudPrimitive primitives, double threshold, double alpha);

    // determine if the 3 points gives a sphere or no 
    bool isSphere(Eigen::Matrix3d vertices, Eigen::Matrix3d normals, double threshold, double alpha) ;
    Eigen::Matrix<double, 1,3> computerCenter(Eigen::Matrix3d vertices, Eigen::Matrix3d normals) ;
    double computerRadius(Eigen::MatrixXd thisVertices, Eigen::Matrix<double, 1,3> thisCenter) ;

// pour la matrix de color a la fin :

   // DecoratedCloud removeCloud(DecoratedCloud cloud, )
}

#endif

