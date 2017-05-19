#include "acq/impl/primitive.hpp"

namespace acq {

    /// ------- computeScore() ------
    void Plane::computeScore(Eigen::Matrix3d var, DecoratedCloud& cloud, double T, double alpha) {

        // --- Compute the Plane Inliers ---
        Eigen::MatrixXi inliers_idx =  this->computeInliers(cloud,T,alpha) ;
        if(inliers_idx.rows() > 0) {
            this->setInliers_idx(inliers_idx);


            // --- Estimate the density of our plane ---
            //int optimalInliers = this->findBestNumberPoints(var, cloud, inliers_idx);

            //std::cout << "optimalInliers : " << optimalInliers << " | inliers_idx.rows(): " << inliers_idx.rows()
            //          << std::endl;

            // --- Compute the plane score ---
            //if the number are closed, the score is high, otherwise it's low
            //double score = 100.0 - double(std::abs(inliers_idx.rows() - optimalInliers)) /
             //                      double(std::max(int(inliers_idx.rows()), optimalInliers)) * 100;

            double score = 100.0 - double(std::abs(inliers_idx.rows() - 121)) /
                                                         double(std::max(int(inliers_idx.rows()), 121)) * 100;

            // --- Set the score for this primitive ---
            this->setScore(score);
            std::cout << " inliers_idx.rows(): " << inliers_idx.rows() << " Score set : " << score << std::endl;
        } else {
            std::cout << "No Inliers found for this primitive" << std::endl;
        }
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

        std::cout << "Basis u " << u << " v " << v << " N " << N << std::endl ;
        
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

        // --- Compute Plane Area ---
        double x_min = inliers2D.col(0).minCoeff();
        double x_max = inliers2D.col(0).maxCoeff();
        double y_min = inliers2D.col(1).minCoeff();
        double y_max = inliers2D.col(1).maxCoeff();

        double thisArea = std::abs(y_max-y_min)*std::abs(x_max-x_min);


        // --- Estimate Optimal Number of Points ---
        double var_x = u.row(0).dot(var.diagonal());
        double var_y = v.row(0).dot(var.diagonal());
        double meanVariance = (pow(var_x,2.0)+pow(var_y,2.0));
        //double meanVariance = std::abs((var_x + var_y)/2.0);
        std::cout << " Area : " << thisArea << " Mean var : " << meanVariance << std::endl;
        int numberPoints = floor(thisArea/(pow(meanVariance,2.0)));

        return numberPoints ;

    }

    /// ------- computeInliers() ------
    Eigen::MatrixXi Plane::computeInliers(DecoratedCloud& cloud, double T, double alpha) {

        int numberPoint = cloud.getVertices().rows();
        Eigen::MatrixXi inliers_idx(numberPoint, 1);

        Eigen::Matrix<double, 1, 3> N = this->getNormal().normalized();
        Eigen::Matrix<double, 1, 3> P = this->getRefPoint();

        if( N.norm() > 0 && numberPoint > 0) {
            double d = -std::abs(N.dot(P));
            const long n = cloud.getVertices().rows();

            Eigen::Matrix<double, 1, 3> _V, _N;

            int idx_counter = 0; double dist = 0;
            for (int i = 0; i < n; i++) {
                _V = cloud.getVertices().row(i);
                _N = cloud.getVertices().row(i).normalized();

                // --- Check if in range and if normals match ---
                //dist = std::abs((_V.dot(N) + d) / N.norm());
                dist = std::abs((N.dot(_V - P)) / N.norm());
                //std::cout << " dist " << dist << " T " << T << " std::abs(_N.dot(N)) " << std::abs(_N.dot(N)) << std::endl;
                if (dist < T && std::abs(_N.dot(N)) > alpha) {
                    inliers_idx(idx_counter, 0) = i;
                    idx_counter++;
                }
            }

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