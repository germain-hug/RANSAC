#include "acq/impl/ransac.hpp"

namespace acq {

    void ransac(DecoratedCloud &cloud, CloudPrimitive &best_primitives, CloudManager &cloudManager, 
                double thresh, double alpha, double thresh_best, int iterationsTotal, int numberSample) {

            // ****************** INITIALISATION   ***********
            int numberOfPoint = cloud.getVertices.rows() ; 
            Eigen::Matrix<int, 3,1> thisSample
            bool prim_detected = false, test_thisSphere, test_thisPlane, test ;
            int bestPrim_idx, nbAllPrim ;
            Primitive& best_prim ;
            double best_score ;
            // compute the variance 
            Eigen::Matrix3d variance = computeVariance(cloud.getVertices()) ;
            // will contain all the primitives created 
            CloudPrimitive allPrimitive ;

            // create the primitive for this iteration 
            for (int i=0 ; i<M; i++) {
                // sample the right amount of point 
                for (int j=0; i<N; j++) {
                    // 
                    thisSample = sample(numberOfPoint) ;

                    // test for the primitive, if they exist : add them in the cloud primitive 
                    test_thisSphere = computeSphere(thisSample, variance,cloud, allPrimitive, thresh, alpha) ;
                    test_thisPlane = computeSphere(thisSample, variance,cloud, allPrimitive, thresh, alpha) ;
                }
                nbAllPrim = allPrimitive.getCloudSize() ;
                // if a primitive has been created in the turn 
                if (nbAllPrim>0) {
                    // get back the best primitive 
                    bestPrim_idx = allPrimitive.findBestScore() ;
                    best_prim = allPrimitive.getPrimitive(bestPrim_idx) ;

                    // test for the score 
                    best_score = best_prim->getScore() ;

                    if (best_score > thresh_best) {
                        best_primitives.addPrimitive(best_prim) ;
                    }

                }

            }

    };

}


