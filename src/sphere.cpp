#include "acq/primitive.h"

namespace acq {
    void Sphere::computeScore(Eigen::Matrix3d variance, DecoratedCloud& cloud, double threshold, double alpha) {
        // compute the inliers 
        Eigen::MatrixXi inliers_idx =  this->computeInliers(cloud,threshold, alpha) ;

        // set the inliers for this primitive 
        this->setInliers_idx(inliers_idx) ;

        // Look how many point are supposed to be there  
        int bestNumber = this->findBestNumberPoints(variance) ;

        // compute the score 
        int numberInliers = inliers_idx.rows() ;

            // --- Compute the plane score ---
            double density_max = 35, score = 0;
            int inliers_min = 20;
            const int n = inliers_idx.rows();

            if(n > inliers_min){
                score = 80 + 0.2*(100.0 - (std::abs(density_max - n)) /
                                double(std::max(density_max, double(n))) * 100.0);
                //std::cout << " inliers_idx.rows(): " << inliers_idx.rows() << " Score set : " << score << std::endl;
            }

        //double score = (double(numberInliers)/double(bestNumber))*100.0 ;

        // set the score for this primitive 
        this->setScore(score) ;
    }

    // find the best number of point for this sphere accordingly to the radius and the variance of points 
    int Sphere::findBestNumberPoints(Eigen::Matrix3d variance) {
        double thisArea = M_PI*4.0*pow(_radius, 2.0) ;
        Eigen::Matrix<double, 1,3> varianceVector = variance.diagonal() ; 
    
        double meanVariance = varianceVector.norm() ;
        double areaAroundPoint = M_PI*pow(meanVariance/4.7, 2.0) ;

        //double areaAroundPoint = M_PI*pow(meanVariance, 2.0) ;

        int numberPoints = floor(thisArea/areaAroundPoint) ;

        return numberPoints ;
    }

    Eigen::MatrixXi Sphere::computeInliers(DecoratedCloud& cloud, double threshold, double alpha) {        
        int numberPoint = cloud.getVertices().rows() ;
        Eigen::Matrix<double, 1,3> thisVertice, thisNormal, estimatedNormal ;
        int index_inliers = 0 ;
        double thisRadius, test1, test2 ;

        Eigen::MatrixXi inliers_idx(numberPoint, 1) ;

        // test for each point if it is in the sphere or not 
        for (int i=0; i < numberPoint; i++) {
            thisVertice = cloud.getVertices().row(i) ;
            thisNormal = cloud.getNormals().row(i) ;

            // compute the estimated normal and radius for this point  
            thisRadius = (thisVertice - _center).norm() ;
            estimatedNormal = (thisVertice - _center).normalized() ;

            // test between the distance and the radius  
            test1 = thisRadius - _radius ;
            test2 = estimatedNormal.dot(thisNormal) ;

            if (std::abs(test1) < threshold ) {
                if ( std::abs(test2) > alpha ) {
                    // if the 2 test are true, the point is an inlier 
                    inliers_idx(index_inliers,0) = i ;
                    index_inliers += 1 ; 
                }
            }   
        }

        // only get back the important part 
        if (index_inliers == 0) {
            inliers_idx = inliers_idx.topRows(1);
        }
        else {
            inliers_idx = inliers_idx.topRows(index_inliers - 1);
        }
        return inliers_idx ;
    }

    Primitive* Sphere::clone(){
        Primitive* thisSphere = new Sphere(this->getRadius(), this->getCenter()) ;
        thisSphere->setType(1) ;
        return thisSphere;
    }

}