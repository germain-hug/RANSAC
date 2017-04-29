#ifndef ACQ_GESTION_H
#define ACQ_GESTION_H

#include "acq/typedefs.h"
#include "decoratedCloud.h"

namespace acq {
    Eigen::MatrixXi sample(int cloudSize, int N); // Renvoie matrice d'index de vertices
    double computeVariance(Eigen::MatrixXd V);
    void computeSphere(Eigen::Matrix3i sample_idx, DecoratedCloud cloud); // test if sphere -> create sphere primitive
    void computePlane(Eigen::Matrix3i sample_idx, DecoratedCloud cloud);


// pour la matrix de color a la fin :

   // DecoratedCloud removeCloud(DecoratedCloud cloud, )
}

#endif

