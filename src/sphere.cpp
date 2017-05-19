#include "acq/primitive.h"

namespace acq {
    void Sphere::computeScore(Eigen::Matrix3d variance, DecoratedCloud& cloud, double threshold, double alpha) {
        // compute the inliers 
        Eigen::MatrixXi inliers_idx =  this->computeInliers(cloud,threshold, alpha) ;

        // set the inliers for this primitive 
        this->setInliers_idx(inliers_idx) ;

        // Look how many point are supposed to be there  
        int bestNumber = this->findBestNumberPoints(variance) ;

        std::cout << "Best number " << bestNumber << std::endl ;
        // compute the score 
        int numberInliers = inliers_idx.rows() ;

        double score = (double(numberInliers)/double(bestNumber))*100.0 ;

        std::cout << "Inliers " << numberInliers << std::endl ;
        std::cout << "score " << score << std::endl ;

        // set the score for this primitive 
        this->setScore(score) ;
    }

    // find the best number of point for this sphere accordingly to the radius and the variance of points 
    int Sphere::findBestNumberPoints(Eigen::Matrix3d variance) {
        double thisArea = M_PI*4.0*pow(_radius, 2.0) ;
        Eigen::Matrix<double, 1,3> varianceVector = variance.diagonal() ; 
    
std::cout << "varianceVector : "<< varianceVector << std::endl ;

        double meanVariance = varianceVector.norm() ;

std::cout << "mean Variance : "<< meanVariance << std::endl ;

        double areaAroundPoint = M_PI*pow(meanVariance/4.7, 2.0) ;

        std::cout << "small area : "<< areaAroundPoint << std::endl ;

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
            estimatedNormal = thisVertice - _center ;
            estimatedNormal = estimatedNormal.normalized() ;

            // test between the distance and the radius  
            test1 = thisRadius - _radius ;
            test2 = estimatedNormal.dot(thisNormal) ;

            if (std::abs(test1) < threshold ) {
                if ( test2 > alpha ) {
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
}