#ifndef ACQ_CLOUDPRIMITIVE_H
#define ACQ_CLOUDPRIMITIVE_H

#include "acq/primitive.h"
#include <vector>

namespace acq {

/** Small class to keep track of multiple primitives */
    class CloudPrimitive {
    public:

        // test constructor 
        CloudPrimitive() {} ;
        ~CloudPrimitive() {} ;

        /** add a primitive to the vector */
        void addPrimitive(Primitive* primitive);

        /** set the primitive to a fixed place */
        void setPrimitive(Primitive*, int index);

        /** get back the primitive from a vector  */
        Primitive* getPrimitive(int index);
        int getCloudSize() {return _primitives.size();};

        // find the primitive with the best score : return the index
        int findBestScore() ;

        // delete primitive
        void deletePrimitive(int index) ;

        void clearAllPrimitives() ;

        void deleteCloudFromIndex(int indexStart) ;
    protected:
        std::vector<Primitive*> _primitives ;

    public:
        // See https://eigen.tuxfamily.org/dox-devel/group__TopicStructHavingEigenMembers.html
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    }; 

} 

#endif 
