#include "acq/primitive.h"

namespace acq {

    /// ------- computeScore() ------
    void Plane::computeScore(Eigen::Matrix3d var, DecoratedCloud& cloud, double T, double alpha) {

        // --- Compute the Plane Inliers ---
        Eigen::MatrixXi inliers_idx =  this->computeInliers(cloud,T,alpha) ;
        this->setInliers_idx(inliers_idx) ;

        // --- Estimate the density of our plane ---
        int optimalInliers = this->findBestNumberPoints(var, cloud, inliers_idx) ;

        // --- Compute the plane score ---
        double score = (double(inliers_idx.rows())/double(optimalInliers))*100.0 ;

        // --- Set the score for this primitive ---
        this->setScore(score) ;
    }

    /// ------- findBestNumberPoints() ------
    int Plane::findBestNumberPoints(Eigen::Matrix3d var, DecoratedCloud& cloud, const Eigen::MatrixXi inliers_idx) {
        /* Compute the optimal number of points for this plane, given the cloud variance and the inliers */

        // --- Find an orthonormal basis of the plane ---
        Eigen::MatrixXd N = this->getNormal();
        Eigen::MatrixXd u; u<< -N(0,1),N(0,0),0.0; u.normalize();
        Eigen::MatrixXd v = u.cross(N); v.normalize();

        // --- Retrieve Inliers and project on (u,v) basis ---
        const int n = inliers_idx.rows();
        Eigen::MatrixXd inliers2D(n, 2);

        for(int i = 0; i<inliers_idx.rows(); i++){
            const int idx = inliers_idx(i,0);
            Eigen::MatrixXd this_vertex = cloud.getVertices().row(idx);
            inliers2D(i,0) = u.row(0).dot(this_vertex.row(0));
            inliers2D(i,1) = v.row(0).dot(this_vertex.row(0));
        }

        // --- Compute Plane Area ---
        double x_min = inliers2D.col(0).minCoeff();
        double x_max = inliers2D.col(0).maxCoeff();
        double y_min = inliers2D.col(1).minCoeff();
        double y_max = inliers2D.col(1).maxCoeff();

        double thisArea = std::abs(y_max-y_min)*std::abs(x_max-x_min);

        // --- Estimate Optimal Number of Points ---
        double meanVariance = (var(1,1)+var(2,2)+var(3,3))/3.0 ;
        int numberPoints = floor(thisArea/pow(meanVariance, 2.0)) ;

        return numberPoints ;
    }

    /// ------- computeInliers() ------
    Eigen::MatrixXi Plane::computeInliers(DecoratedCloud& cloud, double T, double alpha) {

        Eigen::MatrixXd N = this->getNormal();
        Eigen::MatrixXd P = this->getRefPoint();

        double d = N.row(0).dot(P.row(0));
        const long n = cloud.getVertices().rows() ;

        Eigen::MatrixXi inliers_idx;
        Eigen::MatrixXd _V, _N;
        int idx_counter = 0;

        for(int i=0; i<n; i++){
            _V = cloud.getVertices().row(i);
            _N = cloud.getVertices().row(i);

            // --- Check if in range and if normals match ---
            double dist = (_V.row(0).dot(N.row(0)) + d) / N.row(0).norm();
            if(dist < T && _N.row(0).dot(N.row(0)) < alpha){
                inliers_idx(idx_counter) = i;
                idx_counter++;
            }
        }

        return inliers_idx ;
    }
}