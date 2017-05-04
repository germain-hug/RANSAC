#ifndef IGLFRAMEWORK_RANSAC_H
#define IGLFRAMEWORK_RANSAC_H

#include "acq/cloudPrimitive.h"
#include "acq/decoratedCloud.h"
#include "acq/gestion.h"
#include "acq/cloudManager.h"

namespace acq {
    void ransac(DecoratedCloud& cloud, CloudPrimitive& best_primitives, CloudManager& cloudManager, 
                double thresh, double alpha, double thresh_best, int iterationsTotal, int numberSample);
}

#endif

