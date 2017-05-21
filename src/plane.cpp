#include "acq/impl/primitive.hpp"

namespace acq {

    /// ------- computeScore() ------
    void Plane::computeScore(Eigen::Matrix3d var, DecoratedCloud& cloud, double T, double alpha) {

        // --- Compute the Plane Inliers ---
        Eigen::MatrixXi inliers_idx =  this->computeInliers(cloud,T,alpha) ;
        if(inliers_idx.rows() > 0) {
            this->setInliers_idx(inliers_idx);


            // --- Estimate the density of our plane ---
            //double inliersDensity = inliers_idx.rows()/this->findInliersBoundingBox(var, cloud, inliers_idx);

            // --- Compute the plane score ---
            double density_max = 130, score = 0;
            int inliers_min = 40;
            const int n = inliers_idx.rows();

            if(n > inliers_min){
                score = 80 + 0.2*(100.0 - (std::abs(density_max - n)) /
                                double(std::max(density_max, double(n))) * 100.0);
                //std::cout << " inliers_idx.rows(): " << inliers_idx.rows() << " Score set : " << score << std::endl;
            }

            // --- Set the score for this primitive ---
            this->setScore(score);

        } else {
           // std::cout << "No Inliers found for this primitive" << std::endl;
        }
    }

    /// ------- findInliersBoundingBox() ------
    double Plane::findInliersBoundingBox(Eigen::Matrix3d var, DecoratedCloud& cloud, Eigen::MatrixXi inliers_idx) {
        /* Compute the optimal number of points for this plane, given the cloud variance and the inliers */

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

        // --- Retrieve Inliers and project on (u,v) basis ---
        const int n = inliers_idx.rows();
        Eigen::MatrixXd inliers2D(n, 2);
        Eigen::MatrixXd this_vertex ;
        int idx;

        for(int i = 0; i<n; i++){
            idx = inliers_idx(i,0);
            this_vertex = cloud.getVertices().row(idx);
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
        return thisArea;

    }

    /// ------- computeInliers() ------
    Eigen::MatrixXi Plane::computeInliers(DecoratedCloud& cloud, double T, double alpha) {

        int numberPoint = cloud.getVertices().rows();
        Eigen::MatrixXi inliers_idx(numberPoint, 1);

        Eigen::Matrix<double, 1, 3> N = this->getNormal().normalized();
        Eigen::Matrix<double, 1, 3> P = this->getRefPoint();


        if( N.norm() > 0 && numberPoint > 0) {
            const long n = cloud.getVertices().rows();

            Eigen::Matrix<double, 1, 3> _V, _N;

            int idx_counter = 0; double dist = 0;
            for (int i = 0; i < n; i++) {
                _V = cloud.getVertices().row(i);
                _N = cloud.getNormals().row(i).normalized();

                if(_N.dot(N) < 0) _N = -_N;

                // --- Check if in range and if normals match ---
                dist = std::abs((N.dot(_V - P)) / N.norm());
              //  std::cout << " dist " << dist << " T " << T << " std::abs(_N.dot(N)) " << std::abs(_N.dot(N)) << std::endl;
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

    Primitive* Plane::clone(){
        Primitive* thisPlane = new Plane(this->getRefPoint(), this->getNormal()) ;
        thisPlane->setType(2) ;
        return thisPlane;
    }
}