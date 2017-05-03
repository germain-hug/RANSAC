#include "acq/impl/ransac.hpp"

namespace acq {

    void ransac(DecoratedCloud &cloud, CloudPrimitive &best_primitives, CloudManager &cloudManager, double thresh,
                double alpha, double thresh_best) {
        
                    /* 
                    1) definir un number d'iteration totale M
                    2) definir un number de sample par iteration N

-> creer un CloudPrimitive& primitives pour l'ensemble des primitive'

    2bis) boucle for pour M
                    3) boucle for -> faire les sample N
                      4) tester pour chaque sample les 2 primitives 
                      5) si ok les rentrer dans cloud primitive 


    6) prendre  le best prim 
    7) si score suffisant > best_primitives + cloud dans


                    */

    };

}


