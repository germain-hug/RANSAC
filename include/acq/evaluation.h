#ifndef ACQ_EVALUATION_H
#define ACQ_EVALUATION_H

#include "acq/cloudManager.h"
#include "acq/cloudPrimitive.h"

#include <ANN/ANN.h>					// ANN declarations

namespace acq {
// ***********************$ Function to add noise and test smoothing ******************
void computeBoundingBox(float &Xmax,float & Xmin,float & Ymax,float & Ymin,float &Zmax, float & Zmin, DecoratedCloud& cloud) ;

// with 1 add noise on the vertices, with 2 add noise on the normals 
Eigen::MatrixXd addNoise(float noise, DecoratedCloud& cloud, int typeMatrix) ;


// ***********************$ Function to perform connected component algorithm ******************
void connectedComponentManager(CloudManager& thisCloudManager, CloudPrimitive& best_primitives, double threshold) ;

void connectedComponent(DecoratedCloud& cloud, double threshold) ;

void labelVertices(Eigen::RowVector3d thisColor, ANNpointArray verticesArray, Eigen::MatrixXd& colors, 
                    int this_idx, Eigen::MatrixXd& visited, ANNkd_tree*	kdTree, double threshold, int connectivity) ;

ANNpointArray matrixToANNArray(Eigen::MatrixXd const& points) ;

}

#endif

// do connected component 