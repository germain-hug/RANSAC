#include "acq/evaluation.h"

namespace acq {

void computeBoundingBox(float &Xmax,float & Xmin,float & Ymax,float & Ymin,float &Zmax, float & Zmin, DecoratedCloud& cloud) {
    Eigen::MatrixXd maximum(1, 3)  ;
    maximum = cloud.getVertices().colwise().maxCoeff() ;
    Xmax = maximum(0) ;
    Ymax= maximum(1) ;
    Zmax  = maximum(2) ;

    Eigen::MatrixXd minimum(1, 3)  ;
    minimum = cloud.getVertices().colwise().minCoeff() ;
    Xmin = minimum(0) ;
    Ymin = minimum(1) ;
    Zmin = minimum(2) ;
}

Eigen::MatrixXd addNoise(float noise, DecoratedCloud& cloud, int typeMatrix) {
        // compute the value of the boundingBox for the second cloud 
        float Xmax, Xmin, Ymax, Ymin, Zmax, Zmin ;
        computeBoundingBox(Xmax, Xmin, Ymax, Ymin,Zmax, Zmin, cloud) ;

        // set the variance to sigma = noise% of the bouding box size in each direction
        float sigmaX =(Xmax-Xmin)*noise  ;
        float sigmaY = (Ymax-Ymin)*noise ;
        float sigmaZ = (Zmax-Zmin)*noise ;

        int M ;

        if (typeMatrix==1) 
            M = cloud.getVertices().rows() ;
        else if (typeMatrix==2)
            M = cloud.getNormals().rows() ;


        // initialization 
        Eigen::MatrixXd random(M, 3)  ;

         // construct noise 
          for(int i = 0; i< M; i++) {
             random(i, 0) = std::rand()*sigmaX /RAND_MAX ; 
             random(i, 1) = std::rand()*sigmaY /RAND_MAX ;
             random(i, 2) = std::rand()*sigmaZ /RAND_MAX ;
         }
        
        Eigen::MatrixXd newMatrix(M,3) ;

        // add noise to the matrix 
        if (typeMatrix==1) 
            newMatrix = cloud.getVertices() + random ;
        else if (typeMatrix==2)
            newMatrix = cloud.getNormals() + random ;

        return newMatrix ;
}


}


