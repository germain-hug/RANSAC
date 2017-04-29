#ifndef ACQ_PRIMITIVE_H
#define ACQ_PRIMITIVE_H

#include "acq/typedefs.h"

namespace acq {

    /** ---- PRIMITIVE ---- */
    class Primitive {

    public:
        double getScore()const{return _score;}
        Eigen::MatrixXd getInliers_idx()const{return _inliers_idx;}

        void setScore(double score){_score = score;}
        void setInliers_idx(Eigen::MatrixXi inliers_idx){_inliers_idx = inliers_idx;}

        void computeScore(double variance, Eigen::MatrixXd pointCloud){} // Appelle computeInliers
        Eigen::MatrixXd computeInliers(Eigen::MatrixXd pointCloud){}


    protected:
        double _score;
        Eigen::MatrixXi _inliers_idx;
    };

    /** ---- SPHERE ---- */
    class Sphere : public Primitive {

    public:
        Sphere(double radius, Eigen::Matrix3d center) : radius(_radius), center(_center) {} ;
        ~Sphere(){};

        double getRadius()const{return _radius;}
        double getCenter()const{return _center;}
        void setRadius(double radius){_radius = radius;}
        void setCenter(Eigen::Matrix3d center){_center = center;}

        void computeScore(double variance, Eigen::MatrixXd pointCloud){}
        Eigen::MatrixXd computeInliers(Eigen::MatrixXd pointCloud){}

    private:
        double _radius;
        Eigen::Matrix3d _center;
    };

    /** ---- PLANE ---- */
    class Plane : public Primitive {

    public:
        Plane(Eigen::Matrix3d refPoint, Eigen::Matrix3d normal) : refPoint(_refPoint), normal(_normal) {} ;
        ~Plane(){};

        double getNormal()const{return _normal;}
        double getRefPoint()const{return _refPoint;}
        void setNormal(Eigen::Matrix3d normal){_normal = normal;}
        void setNormal(Eigen::Matrix3d refPoint){_refPoint = refPoint;}

        void computeScore(double variance, Eigen::MatrixXd pointCloud){}
        Eigen::MatrixXd computeInliers(Eigen::MatrixXd pointCloud){}

    private:
        Eigen::Matrix3d _refPoint;
        Eigen::Matrix3d _normal;
    };
}

#endif
