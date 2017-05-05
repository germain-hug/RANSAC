#include "acq/impl/primitive.hpp"

namespace acq {

    /// ------- computeScore() ------
    void Plane::computeScore(Eigen::Matrix3d var, DecoratedCloud& cloud, double T, double alpha) {
        std::cout << "Entered in compute score "<< std::endl ;

        // --- Compute the Plane Inliers ---
        Eigen::MatrixXi inliers_idx =  this->computeInliers(cloud,T,alpha) ;
        std::cout << "Inliers : "<< inliers_idx << std::endl ;
        this->setInliers_idx(inliers_idx) ;
        std::cout << "Inliers set " << std::endl ;


        // --- Estimate the density of our plane ---
        int optimalInliers = this->findBestNumberPoints(var, cloud, inliers_idx) ;

        std::cout << "Best number of point : " << optimalInliers << std::endl ;

        // --- Compute the plane score ---
        double score = (double(inliers_idx.rows())/double(optimalInliers))*100.0 ;
        std::cout << "Score finam : " << score << std::endl ;

        // --- Set the score for this primitive ---
        this->setScore(score) ;
        std::cout << "Score set : " << score << std::endl ;
    }

    /// ------- findBestNumberPoints() ------
    int Plane::findBestNumberPoints(Eigen::Matrix3d var, DecoratedCloud& cloud, Eigen::MatrixXi inliers_idx) {
        /* Compute the optimal number of points for this plane, given the cloud variance and the inliers */
        std::cout << "Entered in find best number points "<< std::endl ;

        // --- Find an orthonormal basis of the plane ---
        Eigen::Matrix<double, 1,3> N = this->getNormal();
        Eigen::Matrix<double, 1,3> u(-N(0,1),N(0,0),0.0) ;

        std::cout << "u : " << u << std::endl ;
        u.normalize();
        Eigen::Matrix<double, 1,3> v = u.cross(N); 
        v.normalize();

        std::cout << "Basis for the point u : " << u << "and v : " << v << std::endl ;
        
        // --- Retrieve Inliers and project on (u,v) basis ---
        const int n = inliers_idx.rows();
        Eigen::MatrixXd inliers2D(n, 2);
        Eigen::Matrix<double, 1,3> this_vertex ;
        int idx ;

        for(int i = 0; i<n; i++){
            idx = inliers_idx(i,0);
            this_vertex = cloud.getVertices().row(idx);
            inliers2D(i,0) = u.dot(this_vertex);
            inliers2D(i,1) = v.dot(this_vertex);
        }
        std::cout << "Inliers 2D OK " << std::endl ;

        // --- Compute Plane Area ---
        double x_min = inliers2D.col(0).minCoeff();
        double x_max = inliers2D.col(0).maxCoeff();
        double y_min = inliers2D.col(1).minCoeff();
        double y_max = inliers2D.col(1).maxCoeff();

        double thisArea = std::abs(y_max-y_min)*std::abs(x_max-x_min);

        std::cout << "this Area : " << thisArea << std::endl ;

        // --- Estimate Optimal Number of Points ---
        double meanVariance = (var(1,1)+var(2,2)+var(3,3))/3.0 ;
        int numberPoints = floor(thisArea/pow(meanVariance, 2.0)) ;

        return numberPoints ;

        std::cout << "final number of point : " << numberPoints << std::endl ;

    }

    /// ------- computeInliers() ------
    Eigen::MatrixXi Plane::computeInliers(DecoratedCloud& cloud, double T, double alpha) {
        int numberPoint = cloud.getVertices().rows() ;
        Eigen::RowVector3d N = this->getNormal();
        Eigen::RowVector3d P = this->getRefPoint();

        double d = N.dot(P);
        const long n = cloud.getVertices().rows() ;

        Eigen::MatrixXi inliers_idx(numberPoint, 1);

        Eigen::Matrix<double, 1,3> _V, _N;

        int idx_counter = 0;

        for(int i=0; i<n; i++){
            _V = cloud.getVertices().row(i);
            _N = cloud.getVertices().row(i);

            // --- Check if in range and if normals match ---
            double dist = (_V.dot(N) + d) / N.norm();
            if(dist < T && _N.dot(N) < alpha){
                inliers_idx(idx_counter,0) = i;
                idx_counter++;
            }
        }

        inliers_idx = inliers_idx.topRows(idx_counter-1) ;
        return inliers_idx ;
    }
}