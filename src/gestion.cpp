#include "gestion.h"

namespace acq {

// Renvoie matrice d'index de vertices
Eigen::MatrixXi sample(int cloudSize) {
    Eigen::Matrix<int, 3,1> sampleInd(numberPoint,1) ;
    // add a random indices between 0 and sizeMatrix in a numberPoint sized vector 
    for (int i=0; i<3; i++) {
        int newIndex = rand() % (cloudSize + 1) ;
        sampleInd(i) = newIndex ;
    }
    return sampleInd ;
}

// Return the variance of the point cloud ()
Eigen::Matrix3d computeVariance(Eigen::MatrixXd V) {
    // compute the matrix minus the min of each coordinate
    Eigen::MatrixXd centered = V.rowwise() - V.colwise().mean();

    // compute the covariance matrix 
    Eigen::Matrix3d cov = (centered.adjoint() * centered) / double(V.rows() - 1);
    return cov ;
}

// if the 3 points create a sphere, we add it to the primitives 
bool computeSphere(Eigen::Matrix<int, 3,1> sample_idx, Eigen::Matrix3d variance, DecoratedCloud& cloud, cloudPrimitive primitives, double threshold, double alpha) {
    Eigen::MatrixXd vertices = cloud.getVertices() ;
    Eigen::MatrixXd normals = cloud.getNormals() ;
    int cloudSize = vertices.rows() ;

    Eigen::Matrix3d thisVertices ;
    Eigen::Matrix3d thisNormals ;

    // extract the 3 points sampled by the indices
    for(int i =0 ; i< 3; i++) {
        thisVertices.row(i) = vertices.row(sample_idx(i,1)) ;
        thisNormals.row(i) = normals.row(sample_idx(i,1)) ;
    }

    // test if it's a sphere
    bool isSphere = isSphere(thisVertices, thisNormals, threshold, alpha) ;

    if (isSphere) {
        // compute the attribut for the object 
        Eigen::Matrix<double, 1,3> thisCenter = computerCenter(thisVertices, thisNormals) ;
        double thisRadius = computerRadius(thisVertices, thisCenter) ;

        // create the object and compute its score 
        Sphere thisSphere = Sphere(thisRadius, thisCenter) ;
        thisSphere.computeScore(variance, cloud, threshold, alpha) ;

        // store it in the cloud primitive 
        primitives.addPrimitive(thisSphere) ;

        // the sphere has been accepted : we return true 
        return true ;
    }

    else {
        return false ; 
    }
}

// return true if the 3 points create a valid sphere 
bool isSphere(Eigen::Matrix3d vertices, Eigen::Matrix3d normals, double threshold, double alpha) {
    // estimate the center and the radius using 2 points 
    Eigen::Matrix<double, 1,3> thisCenter = computerCenter(vertices.topRows(2), normals.topRows(2)) ;
    double estimatedRadius = computerRadius(vertices.topRows(2), thisCenter) ;

    // compute the estimated normal for the least point  
    Eigen::Matrix<double, 1,3> estimatedNormal = vertices.row(2) - thisCenter ;
    estimatedNormal = estimatedNormal.normalized() ;

    // test for the radius 
    double test1 = computerRadius(vertices.row(2), thisCenter) - estimatedRadius ;
    double test2 = estimatedNormal.dot(normals.row(2)) ;

    if (test1.abs() < threshold )
        if ( test2 < alpha ) {
            // if the 2 test are true, the 3 points form a sphere  
            return true ; 
        }
        else 
            return false ; 
    else 
        return false ; 
}

// compute the center of a shpere by finding the better intersection possible using least square computation
Eigen::Matrix<double, 1,3> computerCenter(Eigen::MatrixXd vertices, Eigen::MatrixXd normals) {
    Eigen::MatrixXd::Zero(3,3) R ;
    Eigen::Matrix<double, 3,1> q ;
    q = setZero(3,1) ;

    int numberPoint = vertices.rows() ;
    Eigen::Matrix<double, 1,3> thisNormal ;
    Eigen::Matrix<double, 1,3> thisNormal ;


    // fill the system 
    for (int i = 0; i< numberPoint; i++) {
        thisNormal = normals.row(i) ;
        thisPosition = vertices.row(i) ;

        R += eye(3) - thisNormal.transpose()*thisNormal ; 

        q += (eye(3) - thisNormal.transpose()*thisNormal) * thisPosition.transpose() ;
    }

    // solve the system using least jacobi decomposition 
    Eigen::Matrix<double, 3,1> thisCenter ;
    thisCenter = R.jacobiSvd(ComputeThinU | ComputeThinV).solve(q) ;

    return thisCenter.transpose() ;
}

double computerRadius(Eigen::MatrixXd thisVertices, Eigen::Matrix<double, 1,3> thisCenter) {
    // compute the distance between each point and the center
    int numberPoint = thisVertices.rows() ;
    Eigen::Matrix3d centerArray = thisCenter.replicate(numberPoint,1) ;
    Eigen::Matrix<double, numberPoint,1> distances = (thisVertices-centerArray).rowwise().norm() ;

    // compute the mean and return it 
    double meanRadius = distances.colwise().mean() ;
    return meanRadius ;
}

}

