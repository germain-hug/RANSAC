#ifndef ACQ_CLOUDPRIMITIVE_H
#define ACQ_CLOUDPRIMITIVE_H

#include "acq/primitive.h"
#include <vector>

namespace acq {

/** \brief Small class to keep track of multiple point clouds. */
    class CloudPrimitive {
    public:
        /** \brief Append a cloud to the list of clouds. */
        void addPrimitive(Primitive const& primitive);

        /** \brief Overwrite a cloud at a specific index. */
        void setPrimitive(Primitive const& primitive, int index);

        /** \brief Get cloud with specific index. */
        Primitive& getPrimitive(int index);

        /** \brief Get cloud with specific index (const version). */
        Primitive const& getPrimitive(int index) const;

        // find the primitive with the best score : return the index
        int findBestScore() ;

        // delete primitive
        void deletePrimitive(int index) ;

    protected:
        std::vector<Primitive> _primitives ; //!< List of clouds possibly with normals and faces.

    public:
        // See https://eigen.tuxfamily.org/dox-devel/group__TopicStructHavingEigenMembers.html
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    }; //...class CloudManager

} //...ns acq

#endif // ACQ_CLOUDMANAGER_H
