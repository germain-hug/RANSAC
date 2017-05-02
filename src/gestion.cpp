#include "gestion.h"

// Renvoie matrice d'index de vertices
Eigen::MatrixXi sample(int cloudSize, int N) {
    Eigen::MatrixXi sampleInd(numberPoint,1) ;
    // add a random indices between 0 and sizeMatrix in a numberPoint sized vector 
    for (int i=0; i<numberPoint; i++) {
        int newIndex = rand() % (cloudSize + 1) ;
        sampleInd(i) = newIndex ;
    }
    return sampleInd ;
}


