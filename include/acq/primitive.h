#ifndef ACQ_PRIMITIVE_H
#define ACQ_PRIMITIVE_H

#include "acq/typedefs.h"

namespace acq {

    /** ---- PRIMITIVE ---- */
    class Primitive {

    public:
        double getScore()const{return _score;}
        Eigen::MatrixXi getInliers_idx()const{return _inliers_idx;}

        void setScore(double score){_score = score;}
        void setInliers_idx(Eigen::MatrixXi inliers_idx){_inliers_idx = inliers_idx;}

        void computeScore(double variance, Eigen::MatrixXd pointCloud); // Appelle computeInliers
        Eigen::MatrixXi computeInliers(Eigen::MatrixXd pointCloud);


    protected:
        double _score;
        Eigen::MatrixXi _inliers_idx;
    };

    /** ---- SPHERE ---- */
    class Sphere : public Primitive {

    public:
        // constructor and destructor 
        Sphere(double radius, Eigen::Matrix3d center) : _radius(radius), _center(center) {} ;
        ~Sphere(){};

        // getters/setters 
        double getRadius()const{return _radius;}
        Eigen::Matrix3d getCenter()const{return _center;}
        void setRadius(double radius){_radius = radius;}
        void setCenter(Eigen::Matrix3d center){_center = center;}

        void computeScore(double variance, DecoratedCloud& cloud);
        Eigen::MatrixXi computeInliers(DecoratedCloud& cloud);

    private:
        // radius and center of the sphere 
        double _radius;
        Eigen::Matrix3d _center;
    };

    /** ---- PLANE ---- */
    class Plane : public Primitive {

    public:
        Plane(Eigen::Matrix3d refPoint, Eigen::Matrix3d normal) : _refPoint(refPoint), _normal(normal) {} ;
        ~Plane(){};

        Eigen::Matrix3d getNormal()const{return _normal;}
        Eigen::Matrix3d getRefPoint()const{return _refPoint;}
        void setNormal(Eigen::Matrix3d normal){_normal = normal;}
        void setRefPoint(Eigen::Matrix3d refPoint){_refPoint = refPoint;}

        void computeScore(double variance, DecoratedCloud& cloud);
        Eigen::MatrixXi computeInliers(DecoratedCloud& cloud);

    private:
        Eigen::Matrix3d _refPoint;
        Eigen::Matrix3d _normal;
    };
}

#endif
