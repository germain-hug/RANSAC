#ifndef ACQ_EVALUATION_H
#define ACQ_EVALUATION_H

#include "acq/decoratedCloud.h"
#include <ANN/ANN.h>					// ANN declarations

namespace acq {
// ***********************$ Function to add noise and test smoothing ******************
void computeBoundingBox(float &Xmax,float & Xmin,float & Ymax,float & Ymin,float &Zmax, float & Zmin, DecoratedCloud& cloud) ;

// with 1 add noise on the vertices, with 2 add noise on the normals 
Eigen::MatrixXd addNoise(float noise, DecoratedCloud& cloud, int typeMatrix) ;

void connectedComponent(DecoratedCloud& cloud, double threshold) ;

void labelVertices(Eigen::RowVector3d thisColor, ANNpointArray verticesArray, Eigen::MatrixXd& colors, 
                    int this_idx, Eigen::MatrixXd& visited, ANNkd_tree*	kdTree, double threshold) ;

ANNpointArray matrixToANNArray(Eigen::MatrixXd const& points) ;

}

#endif

// do connected component 