#include "acq/impl/primitive.hpp"

namespace acq {

    /// ------- computeScore() ------
    void Plane::computeScore(Eigen::Matrix3d var, DecoratedCloud& cloud, double T, double alpha) {
        std::cout << "Entered in compute score "<< std::endl ;

        // --- Compute the Plane Inliers ---
        Eigen::MatrixXi inliers_idx =  this->computeInliers(cloud,T,alpha) ;
        this->setInliers_idx(inliers_idx) ;
        std::cout << "Inliers set " << std::endl ;


        // --- Estimate the density of our plane ---
        std::cout << "inliers rows = " << inliers_idx.rows() << std::endl;
        int optimalInliers = this->findBestNumberPoints(var, cloud, inliers_idx) ;

        std::cout << "optimalInliers : " << optimalInliers << " | inliers_idx.rows(): "<< inliers_idx.rows() << std::endl ;

        // --- Compute the plane score ---
        //if the number are closed, the score is high, otherwise it's low
        double score = 100.0 - double(std::abs(inliers_idx.rows()-optimalInliers)) /
                double(std::max(int(inliers_idx.rows()),optimalInliers))*100;

        std::cout << "Score final : " << score << "%" << std::endl ;

        // --- Set the score for this primitive ---
        this->setScore(score) ;
        std::cout << "Score set : " << score << std::endl ;
    }

    /// ------- findBestNumberPoints() ------
    int Plane::findBestNumberPoints(Eigen::Matrix3d var, DecoratedCloud& cloud, Eigen::MatrixXi inliers_idx) {
        /* Compute the optimal number of points for this plane, given the cloud variance and the inliers */
        std::cout << "Entered in find best number points "<< std::endl ;

        // --- Find an orthonormal basis of the plane ---
        Eigen::Matrix<double, Eigen::Dynamic, 3> N(1, 3); N = this->getNormal().normalized();
        Eigen::Matrix<double, Eigen::Dynamic, 3> u(1, 3);

        // Compute u so that it's orthogonal to N
        if(N(0,1)!=0 || N(0,0)!=0)
            u << -N(0,1), N(0,0), 0.0;
        else if(N(0,2)!=0 || N(0,0)!=0)
            u << -N(0,2), 0.0, N(0,0);
        else u << 0.0, -N(0,2), N(0,1);
        if(u.norm()!=0.0) u.normalize();

        // Compute v to make a basis (u,v,N)
        Eigen::MatrixXd v = u.row(0).cross(N.row(0));
        if(v.norm()!=0.0) v.normalize();

        std::cout << " u " << u << " N " << N << std::endl;


        std::cout << "Basis for the point u : " << u << "and v : " << v << std::endl ;
        
        // --- Retrieve Inliers and project on (u,v) basis ---
        const int n = inliers_idx.rows();
        Eigen::MatrixXd inliers2D(n, 2);
        Eigen::MatrixXd this_vertex ;
        int idx;


        for(int i = 0; i<n; i++){
            idx = inliers_idx(i,0);
            this_vertex = cloud.getVertices().row(idx); /// ---- NaN !!!!!
            inliers2D(i,0) = u.row(0).dot(this_vertex.row(0));
            inliers2D(i,1) = v.row(0).dot(this_vertex.row(0));
        }
        std::cout << "Inliers 2D OK " << std::endl ;

        // --- Compute Plane Area ---
        double x_min = inliers2D.col(0).minCoeff();
        double x_max = inliers2D.col(0).maxCoeff();
        double y_min = inliers2D.col(1).minCoeff();
        double y_max = inliers2D.col(1).maxCoeff();

        double thisArea = std::abs(y_max-y_min)*std::abs(x_max-x_min);

        std::cout << "this Area : " << thisArea << std::endl;

        // --- Estimate Optimal Number of Points ---
        double var_x = u.row(0).dot(var.diagonal());
        double var_y = v.row(0).dot(var.diagonal());
        //double meanVariance = sqrt(pow(var_x,2.0)+pow(var_y,2.0));
        double meanVariance = std::abs((var_x + var_y)/2.0);

        std::cout << "mean var" << meanVariance << std::endl;
        int numberPoints = floor(thisArea/(meanVariance));
        std::cout << "final number of point : " << numberPoints << std::endl;

        return numberPoints ;

    }

    /// ------- computeInliers() ------
    Eigen::MatrixXi Plane::computeInliers(DecoratedCloud& cloud, double T, double alpha) {
        std::cout << "ah coucou debut func plaaaan" << std::endl ;

        int numberPoint = cloud.getVertices().rows();
        Eigen::MatrixXi inliers_idx(numberPoint, 1);

        Eigen::Matrix<double, 1, 3> N = this->getNormal().normalized();
        Eigen::Matrix<double, 1, 3> P = this->getRefPoint();

        std::cout << "N : " << N << std::endl;

        if( N.norm() > 0 && numberPoint > 0) {
            double d = N.dot(P);
            const long n = cloud.getVertices().rows();

            Eigen::Matrix<double, 1, 3> _V, _N;

            int idx_counter = 0;
            for (int i = 0; i < n; i++) {
                _V = cloud.getVertices().row(i);
                _N = cloud.getVertices().row(i).normalized();

                // --- Check if in range and if normals match ---
                double dist = std::abs((_V.dot(N) + d) / N.norm());
               // std::cout << "dist : " << dist << std::endl;

                if (dist < T && std::abs(_N.dot(N)) < alpha) {
                    inliers_idx(idx_counter, 0) = i;
                    idx_counter++;
                }
            }
            std::cout << "Number index : " << idx_counter << std::endl;
            
            if (idx_counter == 0) {
                inliers_idx = inliers_idx.topRows(1);
            }
            else {
                inliers_idx = inliers_idx.topRows(idx_counter - 1) ;
            }
        }

        return inliers_idx ;
    }
}