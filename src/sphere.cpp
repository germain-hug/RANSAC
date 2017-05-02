#include "acq/primitive.h"

namespace acq {
    void Sphere::computeScore(Eigen::Matrix3d variance, DecoratedCloud& cloud, double threshold, double alpha) {
        // Look how many point are supposed to be there  
        int bestNumber = this->findBestNumberPoints(variance) ;

        // compute the inliers 
        Eigen::MatrixXi inliers_idx =  this->computeInliers(cloud,threshold, alpha) ;

        // set the inliers for this primitive 
        this->setInliers_idx(inliers_idx) ;

        // compute the score 
        int numberInliers = inliers_idx.rows() ;
        double score = (double(numberInliers)/double(bestNumber))*100.0 ;

        // set the score for this primitive 
        this->setScore(score) ;
    }

    // find the best number of point for this sphere accordingly to the radius and the variance of points 
    int Sphere::findBestNumberPoints(Eigen::Matrix3d variance) {
        double thisArea = M_PI*4.0*pow(_radius, 2.0) ;
        double meanVariance = (variance(1,1)+variance(2,2)+variance(3,3))/3.0 ;

        double areaAroundPoint = M_PI*pow(meanVariance, 2.0) ;
        int numberPoints = floor(thisArea/areaAroundPoint) ;

        return numberPoints ;
    }

    Eigen::MatrixXi Sphere::computeInliers(DecoratedCloud& cloud, double threshold, double alpha) {
        int numberPoint = cloud.getVertices().rows() ;
        Eigen::Matrix<double, 1,3> thisVertice, thisNormal, estimatedNormal ;
        Eigen::Matrix<double, 1,3> thisNormal ;
        int index_inliers = 0 ;
        double thisRadius, test1, test2 ;

        Eigen::Matrix<int,numberPoint, 1> inliers ;

        // test for each point if it is in the sphere or not 
        for (int i=0; i < numberPoint; i++) {
            thisVertice = cloud.getVertices().row(i) ;
            thisNormal = cloud.getNormals().row(i) ;

            // compute the estimated normal and radius for this point  
            thisRadius = computerRadius(thisVertice, _center) ;
            estimatedNormal = thisVertice - _center ;
            estimatedNormal = estimatedNormal.normalized() ;

            // test between the distance and the radius  
            test1 = thisRadius - _radius ;
            test2 = estimatedNormal.dot(thisNormal) ;

            if (test1.abs() < threshold ) {
                if ( test2 < alpha ) {
                    // if the 2 test are true, the point is an inlier 
                    inliers.row(index_inliers) = i ;
                    index_inliers += 1 ; 
                }
            }   
        }
        // only get back the important part 
        inliers = inliers.topRows(index_inliers) ;

        return inliers ;
    }
}