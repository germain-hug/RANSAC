#include "acq/impl/ransac.hpp"

namespace acq {

    bool ransac(DecoratedCloud &cloud, CloudPrimitive &best_primitives, CloudManager &cloudManager, 
                double thresh, double alpha, double thresh_best, int iterationsTotal, int numberSample) {

            // ****************** INITIALISATION   ***********
            int numberOfPoint = cloud.getVertices().rows() ; 
            Eigen::Matrix<int, 3,1> thisSample ;
            Eigen::MatrixXi thisInliers ; 
            bool prim_detected = false, test_thisSphere, test_thisPlane, test ;
            int bestPrim_idx, nbAllPrim ;
            double best_score ;

            // compute the variance 
            Eigen::Matrix3d variance = computeVariance(cloud.getVertices()) ;
            /*std::cout << "thresh : "<< thresh  << std::endl ; 
            std::cout << "alpha : "<< alpha  << std::endl ; 
            std::cout << "thresh_best : "<< thresh_best  << std::endl ; 
            std::cout << "iterationsTotal : "<< iterationsTotal  << std::endl ; 
            std::cout << "numberSample : "<< numberSample  << std::endl ; */

            // will contain all the primitives created 
            CloudPrimitive allPrimitive ;

            // create the primitive for this iteration 
            for (int i=0 ; i<iterationsTotal; i++) {
                // sample the right amount of point 
                for (int j=0; j<numberSample; j++) {
                    // sample the point 
                    thisSample = sample(numberOfPoint) ;

                    // test for the primitive, if they exist : add them in the cloud primitive 
                     computeSphere(thisSample, variance, cloud, allPrimitive, thresh, alpha) ;
                    //computePlane(thisSample, variance, cloud, allPrimitive, thresh, alpha) ;
                }

                nbAllPrim = allPrimitive.getCloudSize() ;
                        std::cout << "nb prim : "<< nbAllPrim << std::endl ; 

                // if a primitive has been created in the turn 
                if (nbAllPrim>0) {
                    // get back the best primitive 
                    bestPrim_idx = allPrimitive.findBestScore() ;
                    Primitive* best_prim = allPrimitive.getPrimitive(bestPrim_idx) ;

                    // test for the score 
                    best_score = best_prim->getScore() ;

                    // store the results both in primitives and clou
                    if (best_score > thresh_best) {
                         std::cout << " Enter if  " << std::endl ; 

                        best_primitives.addPrimitive(best_prim) ;

                        std::cout << "cloud size : " <<best_primitives.getCloudSize() << std::endl ; 
                        int thisType = best_prim->getType() ;

                        std::cout << " type " << thisType << std::endl ;

                        std::cout << "ah coucou" << std::endl ;
                        thisInliers = best_prim->computeInliers(cloud, thresh, alpha) ;
                        std::cout << " Inliers computed  " << thisInliers << std::endl ;

                        cleanCloud(cloud, cloudManager, thisInliers) ;

                        std::cout << "Cloud clean  " << std::endl ; 
                        throw std::exception() ;
                    }
                }
            }

            std::cout << "sortie de RANSAC" << std::endl ;

            // cloudManager and cloudPrimitive contains the result of the function 
    };

}


