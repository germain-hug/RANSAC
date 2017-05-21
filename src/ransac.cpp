#include "acq/impl/ransac.hpp"

namespace acq {

    bool ransac(DecoratedCloud &cloud, CloudPrimitive &best_primitives, CloudManager &cloudManager, 
                double thresh, double alpha, double thresh_best, int iterationsTotal, int numberSample) {

            // ****************** INITIALISATION   ***********
            int numberOfPoint = cloud.getVertices().rows() ; 
            Eigen::Matrix<int, 3,1> thisSample ;
            Eigen::MatrixXi thisInliers ; 
            bool prim_detected = false, test_thisSphere, test_thisPlane, test ;
            int bestPrim_idx, nbAllPrim, i ;
            double best_score ;
            int newSize = 4, n_inliers ;
            bool primitiveFound = false ;
            // compute the variance 
            Eigen::Matrix3d variance = computeVariance(cloud.getVertices()) ;

            // will contain all the primitives created 
            CloudPrimitive allPrimitive ;

            // create the primitive for this iteration 
            for (i=0 ; i<iterationsTotal; i++) {
                // sample the right amount of point 
                for (int j=0; j<numberSample; j++) {
                    // sample the point 
                    thisSample = sample(numberOfPoint) ;

                    // test for the primitive, if they exist : add them in the cloud primitive 
                    computeSphere(thisSample, variance, cloud, allPrimitive, thresh, alpha) ;
                    //computePlane(thisSample, variance, cloud, allPrimitive, thresh, alpha);
                }

                int test = i ;
               // std::cout << "iter : " << test  << std::endl;

                nbAllPrim = allPrimitive.getCloudSize() ;

                // if a primitive has been created in the turn 
                if (nbAllPrim>0) {
                    // get back the best primitive 
                    bestPrim_idx = allPrimitive.findBestScore() ;
                    Primitive* best_prim = allPrimitive.getPrimitive(bestPrim_idx) ;

                    // test for the score
                    best_score = best_prim->getScore() ;

                    // store the results both in primitives and cloud
                    if (best_score > thresh_best) {
                        thisInliers = best_prim->computeInliers(cloud, thresh, alpha) ;                     
                        n_inliers = thisInliers.rows();

                        if(n_inliers > 1) {
                            // copy the primitive to store and add it to the newCloud                           
                            Primitive* prim_Storage = best_prim->clone() ;
                            best_primitives.addPrimitive(prim_Storage) ;

                            cleanCloud(cloud, cloudManager, thisInliers) ;
                            newSize = cloud.getVertices().rows() ;
                            primitiveFound = true ;
                        }
                        allPrimitive.deletePrimitive(bestPrim_idx) ;
                    }
                    else {
                        // if the primitive isn't good enough, not take into account
                        allPrimitive.deletePrimitive(bestPrim_idx) ;
                    }

                    if (newSize < 3) {
                        break ;
                    }
                }
            }            
            // free the memory allocated with all the primitives not used 
            allPrimitive.clearAllPrimitives() ;
std::cout << "size best prim  in RANSAC : " << best_primitives.getCloudSize() << std::endl ;
            // cloudManager and cloudPrimitive contains the result of the function
            return primitiveFound ; // Just return a bool
    };

}


