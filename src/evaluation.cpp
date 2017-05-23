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

void connectedComponentManager(CloudManager& thisCloudManager, double threshold) {
    int sizeCloud = thisCloudManager.getCloudSize() ;

    for (int i = 0; i< sizeCloud; i++) {
        connectedComponent(thisCloudManager.getCloud(i), threshold) ;
    } 
}


void connectedComponent(DecoratedCloud& cloud, double threshold) {
    int nbPoints = cloud.getVertices().rows() ;
    int dim = cloud.getVertices().cols() ;
    Eigen::MatrixXd colors(nbPoints,dim) ;
    Eigen::MatrixXd visited(nbPoints,1) ;
    visited = Eigen::MatrixXd::Zero(nbPoints, 1) ;
    Eigen::MatrixXd vertices(nbPoints,dim) ;
    vertices = cloud.getVertices() ;

    Eigen::RowVector3d thisLabel = Eigen::RowVector3d(std::rand()/double(RAND_MAX),
                                       std::rand()/double(RAND_MAX),
                                       std::rand()/double(RAND_MAX)) ;

    // construct the kdTree for this cloud 
    ANNpointArray verticesArray ;  
    verticesArray = annAllocPts(nbPoints, dim);
    verticesArray = matrixToANNArray(cloud.getVertices()) ; 
        
    ANNkd_tree*	kdTree = new ANNkd_tree(					// build search structure
					verticesArray,					// the data points
					nbPoints,						// number of points
					dim);						// dimension of space

    for (int i = 0 ; i < nbPoints; i++) {
        if (visited(i,0) == 0) {
            labelVertices(thisLabel, verticesArray, colors, i, visited, kdTree, threshold) ;

            // change the label 
            thisLabel = Eigen::RowVector3d(std::rand()/double(RAND_MAX),
                                       std::rand()/double(RAND_MAX),
                                       std::rand()/double(RAND_MAX)) ;
        }
    }
    cloud.setColors(colors) ;

    // free the memory used     
    delete kdTree ;
    annDeallocPts(verticesArray) ;
}

void labelVertices(Eigen::RowVector3d thisColor, ANNpointArray verticesArray, Eigen::MatrixXd& colors, 
                    int this_idx, Eigen::MatrixXd& visited, ANNkd_tree*	kdTree, double threshold) {

        // store the result for this vertices, visited and the color  
        visited(this_idx, 0) = 1 ; 
        colors.row(this_idx) = thisColor ;

        // chose to selectione the closest 4 points 
        // allocate variable for the kdTree
        ANNidxArray	nnIdx;	
        nnIdx = new ANNidx[4]; // the index 

	    ANNdistArray dists;		
	    dists = new ANNdist[4]; //the distance
        ANNpoint queryPt;
        queryPt = annAllocPt(3);
        queryPt = verticesArray[this_idx] ; // the point

        kdTree->annkSearch(				// search
				queryPt,						// query point
				4,								// number of near neighbors
				nnIdx,							// nearest neighbors (returned)
				dists,							// distance (returned)
				0.0);

        // take the results
        int indexFind ;
        double distance ;

        for (int thisNeigh = 0; thisNeigh< 4; thisNeigh++) {
            // take the result 
            indexFind =  nnIdx[thisNeigh] ;
            distance =  dists[thisNeigh] ;

            if (visited(indexFind, 0) == 0) {
                if (distance < threshold) {
                    labelVertices(thisColor, verticesArray, colors, indexFind,visited, kdTree, threshold) ; 
                }
                else 
                    std::cout << "distance : " << distance <<std::endl ;
            }
        }

        // free memory 
        delete [] nnIdx ;
        delete [] dists;
        annDeallocPt(queryPt) ;
 }

// function to transform a matrix into a ANNpointArray
ANNpointArray matrixToANNArray(Eigen::MatrixXd const& points) {
    unsigned int M = points.rows() ;
    int dim = 3 ;

	ANNpointArray dataPts;	

    dataPts = annAllocPts(M, dim);			// allocate data points

    for (int i=0; i<M; i++) {
        ANNpoint point ;
        point = annAllocPt(dim) ;

        for (int j=0; j<dim; j++) {
            point[j] = points(i,j) ;
        }
        dataPts[i] = point ;
    }
    return dataPts ;
}

}


