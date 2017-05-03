#include "acq/impl/gestion.hpp"

namespace acq {

//  ****** ============ Helper Functions  =============== ******* 
// Sample the point cloud 
Eigen::MatrixXi sample(int cloudSize) {
    // we want to have 3 points sampled 
    Eigen::Matrix<int, 3,1> sampleInd(3,1) ;
    // add a random indices between 0 and sizeMatrix in a numberPoint sized vector 
    for (int i=0; i<3; i++) {
        int newIndex = rand() % (cloudSize + 1) ;
        sampleInd(i) = newIndex ;
    }
    return sampleInd ;
}

// Return the variance of the point cloud 
Eigen::Matrix3d computeVariance(Eigen::MatrixXd V) {
    // compute the matrix minus the min of each coordinate
    Eigen::MatrixXd centered = V.rowwise() - V.colwise().mean();

    // compute the covariance matrix 
    Eigen::Matrix3d cov = (centered.adjoint() * centered) / double(V.rows() - 1);
    return cov ;
}

//  ****** ============ Functions to handle a sphere =============== ******* 

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

    if (std::abs(test1) < threshold )
        if ( test2 < alpha ) {
            // if the 2 test are true, the 3 points form a sphere  
            return true ; 
        }
        else 
            return false ; 
    else 
        return false ; 
}

// if the 3 points create a sphere, we add it to the primitives 
bool computeSphere(Eigen::Matrix<int, 3,1> sample_idx, Eigen::Matrix3d variance, DecoratedCloud& cloud, CloudPrimitive primitives, double threshold, double alpha) {
    Eigen::MatrixXd vertices = cloud.getVertices() ;
    Eigen::MatrixXd normals = cloud.getNormals() ;
    int cloudSize = vertices.rows() ;

    Eigen::Matrix3d thisSampledVertices ;
    Eigen::Matrix3d thisSampledNormals ;

    // extract the 3 points sampled by the indices
    for(int i =0 ; i< 3; i++) {
        thisSampledVertices.row(i) = vertices.row(sample_idx(i,1)) ;
        thisSampledNormals.row(i) = normals.row(sample_idx(i,1)) ;
    }

    // test if it's a sphere
    bool is_sphere = isSphere(thisSampledVertices, thisSampledNormals, threshold, alpha) ;

    if (is_sphere) {
        // compute the attribut for the object 
        Eigen::Matrix<double, 1,3> thisCenter = computerCenter(thisSampledVertices, thisSampledNormals) ;
        double thisRadius = computerRadius(thisSampledVertices, thisCenter) ;

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

// compute the center of a shpere by finding the better intersection possible using least square computation
Eigen::Matrix<double, 1,3> computerCenter(Eigen::MatrixXd vertices, Eigen::MatrixXd normals) {
    Eigen::Matrix3d R = Eigen::Matrix3d::Zero(3,3) ;
    Eigen::Matrix<double, 3,1> q = Eigen::MatrixXd::Zero(3,1) ;
    Eigen::Matrix3d I = Eigen::Matrix3d::Identity(3,3) ;

    int numberPoint = vertices.rows() ;
    Eigen::Matrix<double, 1,3> thisNormal ;
    Eigen::Matrix<double, 1,3> thisPosition ;

    // fill the system 
    for (int i = 0; i< numberPoint; i++) {
        thisNormal = normals.row(i) ;
        thisPosition = vertices.row(i) ;

        R += I - thisNormal.transpose()*thisNormal ; 

        q += (I - thisNormal.transpose()*thisNormal) * thisPosition.transpose() ;
    }

    // solve the system using least jacobi decomposition 
    Eigen::Matrix<double, 3,1> thisCenter ;


    thisCenter = R.jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(q) ;

    return thisCenter.transpose() ;
}

double computerRadius(Eigen::MatrixXd thisVertices, Eigen::Matrix<double, 1,3> thisCenter) {
    // compute the distance between each point and the center
    int numberPoint = thisVertices.rows() ;
    Eigen::Matrix3d centerArray = thisCenter.replicate(numberPoint,1) ;
    Eigen::MatrixXd distances(numberPoint,1) ; 
    distances = (thisVertices-centerArray).rowwise().norm() ;

    // compute the mean and return it 
    double meanRadius = distances.mean() ;
    return meanRadius ;
}

/********* ============= Functions to handle PLANE =============== *********/
    bool computePlane(Eigen::Matrix3i sample_idx,
                      Eigen::Matrix3d variance,
                      DecoratedCloud &cloud,
                      CloudPrimitive &primitives,
                      double T) {
        Eigen::MatrixXd V = cloud.getVertices(), N = cloud.getNormals();
        const int cloudSize = V.rows();
        const int nSamples = sample_idx.rows();

        // ---- Retrieve the N vertices and their normals ----
        Eigen::Matrix3d thisVertex, thisNormal;
        for (int i = 0; i < 3; i++) {
            thisVertex.row(i) = V.row(sample_idx(i, 1));
            thisNormal.row(i) = N.row(sample_idx(i, 1));
        }

        if (isPlane(thisVertex, thisNormal, T, alpha)) {

            // ---- Create a new plane and compute its score ----
            Eigen::RowVector3d planeNormal = computeNormal(thisVertex, thisNormal);
            Eigen::RowVector3d planeRefPoint = V.colwise.mean();
            Plane newPlane = Plane(planeRefPoint, planeNormal);
            newPlane.computeScore(variance, cloud, threshold, alpha);

            // ---- Store it in the cloudPrimitive ----
            primitives.addPrimitive(newPlane);
            return 1;
        } else {
            return 0;
        }
    }

    /*** ----------- isPlane() ----------- ***/
    bool isPlane(Eigen::MatrixXd V, Eigen::MatrixXd N, double T, double alpha) {
        Eigen::RowVector3d planeNormal = computeNormal(thisVertex, thisNormal);
        bool isPlane = true;
        for (int i = 0; i < N.rows(); i++) {
            if (N.row(i).cross(planeNormal) < T) isPlane = false;
        }
        return isPlane;
    }

    /*** ----------- computeNormal() ----------- ***/
    Eigen::RowVector3d computeNormal(Eigen::MatrixXd V, Eigen::MatrixXd _N) {
        Eigen::RowVector3d N = Eigen::RowVector3d::Zero();
        for (int i = 0; i < V.rows() - 2; i++) {
            Eigen::RowVector3d P01 = V.row(1 + i) - V.row(i);
            Eigen::RowVector3d P02 = V.row(2 + i) - V.row(i);
            N += P02.dot(P01) / (V.rows() - 2);
        }
        if (_N.row(0).cross(N) < 0) N = -N;
        return N;
    }


}


