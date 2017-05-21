#ifndef ACQ_EVALUATION_H
#define ACQ_EVALUATION_H

#include "acq/decoratedCloud.h"

namespace acq {
// ***********************$ Function to add noise and test smoothing ******************
void computeBoundingBox(float &Xmax,float & Xmin,float & Ymax,float & Ymin,float &Zmax, float & Zmin, DecoratedCloud& cloud) ;

// with 1 add noise on the vertices, with 2 add noise on the normals 
Eigen::MatrixXd addNoise(float noise, DecoratedCloud& cloud, int typeMatrix) ;



}

#endif

// do connected component 