#ifndef ACQ_PRIMITIVE_H
#define ACQ_PRIMITIVE_H

#include "acq/typedefs.h"
#include "acq/decoratedCloud.h"

#include <math.h>
#include <cmath>

namespace acq {

    /** ---- PRIMITIVE ---- */
    class Primitive {

    public:
        // constructor/destructor
        Primitive() {} ;
        virtual ~Primitive() {} ;

        // getters.setters 
        double getScore()const{return _score;}
        int getType()const{return _type;}
        Eigen::MatrixXi getInliers_idx()const{return _inliers_idx;}
        void setScore(double score){_score = score;}
        void setInliers_idx(Eigen::MatrixXi inliers_idx){_inliers_idx = inliers_idx;}

        //
        virtual double getRadius(){};
        virtual Eigen::Matrix<double, 1,3> getCenter(){};
        virtual Eigen::RowVector3d getNormal(){};
        virtual Eigen::RowVector3d getRefPoint(){};
        virtual Eigen::MatrixXi computeInliers(DecoratedCloud& cloud, double threshold, double alpha){} ;

        void computeScore(Eigen::Matrix3d variance, DecoratedCloud& pointCloud, double threshold, double alpha); 
        int findBestNumberPoints(Eigen::Matrix3d variance) ;
    
    protected:
        double _score; int _type; // 1: Sphere, 2: Plane
        Eigen::MatrixXi _inliers_idx;
    };

    /** ---- SPHERE ---- */
    class Sphere : public Primitive {

    public:
        // constructor and destructor 
        Sphere(double radius, Eigen::Matrix<double, 1,3> center, int _type = 1) : _radius(radius), _center(center) {} ;
        ~Sphere(){};

        // getters/setters 
        double getRadius()const{return _radius;}
        Eigen::Matrix<double, 1,3> getCenter()const{return _center;}
        void setRadius(double radius){_radius = radius;}
        void setCenter(Eigen::Matrix<double, 1,3> center){_center = center;}

        void computeScore(Eigen::Matrix3d variance, DecoratedCloud& cloud, double threshold, double alpha);
        Eigen::MatrixXi computeInliers(DecoratedCloud& cloud, double threshold, double alpha);
        int findBestNumberPoints(Eigen::Matrix3d variance) ;

    private:
        // radius and center of the sphere 
        double _radius;
        Eigen::Matrix<double, 1,3> _center;
    };

    /** ---- PLANE ---- */
    class Plane : public Primitive {

    public:
        Plane(Eigen::RowVector3d refPoint, Eigen::RowVector3d normal, int _type = 2) : _refPoint(refPoint), _normal(normal) {} ;
        ~Plane(){};

        Eigen::RowVector3d getNormal()const{return _normal;}
        Eigen::RowVector3d getRefPoint()const{return _refPoint;}
        void setNormal(Eigen::RowVector3d normal){_normal = normal;}
        void setRefPoint(Eigen::RowVector3d refPoint){_refPoint = refPoint;}

        void computeScore(Eigen::Matrix3d variance, DecoratedCloud& cloud, double T, double alpha);
        Eigen::MatrixXi computeInliers(DecoratedCloud& cloud, double T, double alpha);
        int findBestNumberPoints(Eigen::Matrix3d var, DecoratedCloud& cloud, Eigen::MatrixXi inliers_idx);

    private:
        Eigen::RowVector3d _refPoint;
        Eigen::RowVector3d _normal;
    };
}

#endif
