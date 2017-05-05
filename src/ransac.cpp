#include "acq/impl/ransac.hpp"

namespace acq {

    void ransac(DecoratedCloud &cloud, CloudPrimitive &best_primitives, CloudManager &cloudManager, 
                double thresh, double alpha, double thresh_best, int iterationsTotal, int numberSample) {

            std::cout << "Inter into ransac" << std::endl ;

            // ****************** INITIALISATION   ***********
            int numberOfPoint = cloud.getVertices().rows() ; 
            Eigen::Matrix<int, 3,1> thisSample ;
            Eigen::MatrixXi thisInliers ; 
            bool prim_detected = false, test_thisSphere, test_thisPlane, test ;
            int bestPrim_idx, nbAllPrim ;
            double best_score ;

            std::cout << "Variable initialized" << std::endl ;

            // compute the variance 
            Eigen::Matrix3d variance = computeVariance(cloud.getVertices()) ;

            std::cout << "Variance : " << variance << std::endl ;

            // will contain all the primitives created 
            CloudPrimitive allPrimitive ;

            std::cout << "CloudPrimitive created " << std::endl ;

            // create the primitive for this iteration 
            for (int i=0 ; i<iterationsTotal; i++) {
                // sample the right amount of point 
                for (int j=0; i<numberSample; j++) {
                    // sample the point 
                    thisSample = sample(numberOfPoint) ;

                     std::cout << "Sample idx :  " << thisSample << std::endl ;

                    // test for the primitive, if they exist : add them in the cloud primitive 
                    computeSphere(thisSample, variance, cloud, allPrimitive, thresh, alpha) ;

                    std::cout << "Sphere test OK pour i= "<< i << "et j "<< j << std::endl ;

                    computePlane(thisSample, variance, cloud, allPrimitive, thresh, alpha) ;
                    
                    std::cout << "Plan test OK  pour i= "<< i << "et j "<< j <<  std::endl ;

                }
                nbAllPrim = allPrimitive.getCloudSize() ;
                // if a primitive has been created in the turn 
                if (nbAllPrim>0) {
                    // get back the best primitive 
                    bestPrim_idx = allPrimitive.findBestScore() ;
                    Primitive& best_prim = allPrimitive.getPrimitive(bestPrim_idx) ;

                    // test for the score 
                    best_score = best_prim.getScore() ;

                    // store the results both in primitives and clou
                    if (best_score > thresh_best) {
                        best_primitives.addPrimitive(best_prim) ;
                        thisInliers = best_prim.computeInliers(cloud, thresh, alpha) ;

                        cleanCloud(cloud, cloudManager, thisInliers) ;
                    }
                }
            }

            // cloudManager and cloudPrimitive contains the result of the function 
    };

}


