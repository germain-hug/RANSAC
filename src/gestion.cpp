#include "acq/impl/gestion.hpp"

namespace acq {

//  ****** ============ Helper Functions  =============== ******* 
// Sample the point cloud 
Eigen::MatrixXi sample(int cloudSize) {
    Eigen::Matrix<int, 3,1> sampleInd(numberPoint,1) ;
    // add a random indices between 0 and sizeMatrix in a numberPoint sized vector 
    for (int i=0; i<3; i++) {
        int newIndex = rand() % (cloudSize + 1) ;
        sampleInd(i) = newIndex ;
    }
    return sampleInd ;
}

// Return the variance of the point cloud 
Eigen::Matrix3d computeVariance(Eigen::MatrixXd V) {
    // compute the matrix minus the min of each coordinate
    Eigen::MatrixXd centered = V.rowwise() - V.colwise().mean();

    // compute the covariance matrix 
    Eigen::Matrix3d cov = (centered.adjoint() * centered) / double(V.rows() - 1);
    return cov ;
}

//  ****** ============ Functions to handle a sphere =============== ******* 

// return true if the 3 points create a valid sphere 
bool isSphere(Eigen::Matrix3d vertices, Eigen::Matrix3d normals, double threshold, double alpha) {
    // estimate the center and the radius using 2 points 
    Eigen::Matrix<double, 1,3> thisCenter = computerCenter(vertices.topRows(2), normals.topRows(2)) ;
    double estimatedRadius = computerRadius(vertices.topRows(2), thisCenter) ;

    // compute the estimated normal for the least point  
    Eigen::Matrix<double, 1,3> estimatedNormal = vertices.row(2) - thisCenter ;
    estimatedNormal = estimatedNormal.normalized() ;

    // test for the radius 
    double test1 = computerRadius(vertices.row(2), thisCenter) - estimatedRadius ;
    double test2 = estimatedNormal.dot(normals.row(2)) ;

    if (test1.abs() < threshold )
        if ( test2 < alpha ) {
            // if the 2 test are true, the 3 points form a sphere  
            return true ; 
        }
        else 
            return false ; 
    else 
        return false ; 
}

// if the 3 points create a sphere, we add it to the primitives 
bool computeSphere(Eigen::Matrix<int, 3,1> sample_idx, Eigen::Matrix3d variance, DecoratedCloud& cloud, cloudPrimitive primitives, double threshold, double alpha) {
    Eigen::MatrixXd vertices = cloud.getVertices() ;
    Eigen::MatrixXd normals = cloud.getNormals() ;
    int cloudSize = vertices.rows() ;

    Eigen::Matrix3d thisSampledVertices ;
    Eigen::Matrix3d thisSampledNormals ;

    // extract the 3 points sampled by the indices
    for(int i =0 ; i< 3; i++) {
        thisSampledVertices.row(i) = vertices.row(sample_idx(i,1)) ;
        thisSampledNormals.row(i) = normals.row(sample_idx(i,1)) ;
    }

    // test if it's a sphere
    bool isSphere = isSphere(thisSampledVertices, thisSampledNormals, threshold, alpha) ;

    if (isSphere) {
        // compute the attribut for the object 
        Eigen::Matrix<double, 1,3> thisCenter = computerCenter(thisSampledVertices, thisSampledNormals) ;
        double thisRadius = computerRadius(thisSampledVertices, thisCenter) ;

        // create the object and compute its score 
        Sphere thisSphere = Sphere(thisRadius, thisCenter) ;
        thisSphere.computeScore(variance, cloud, threshold, alpha) ;

        // store it in the cloud primitive 
        primitives.addPrimitive(thisSphere) ;

        // the sphere has been accepted : we return true 
        return true ;
    }

    else {
        return false ; 
    }
}

// compute the center of a shpere by finding the better intersection possible using least square computation
Eigen::Matrix<double, 1,3> computerCenter(Eigen::MatrixXd vertices, Eigen::MatrixXd normals) {
    Eigen::MatrixXd::Zero(3,3) R ;
    Eigen::Matrix<double, 3,1> q ;
    q = setZero(3,1) ;

    int numberPoint = vertices.rows() ;
    Eigen::Matrix<double, 1,3> thisNormal ;
    Eigen::Matrix<double, 1,3> thisNormal ;


    // fill the system 
    for (int i = 0; i< numberPoint; i++) {
        thisNormal = normals.row(i) ;
        thisPosition = vertices.row(i) ;

        R += eye(3) - thisNormal.transpose()*thisNormal ; 

        q += (eye(3) - thisNormal.transpose()*thisNormal) * thisPosition.transpose() ;
    }

    // solve the system using least jacobi decomposition 
    Eigen::Matrix<double, 3,1> thisCenter ;
    thisCenter = R.jacobiSvd(ComputeThinU | ComputeThinV).solve(q) ;

    return thisCenter.transpose() ;
}

double computerRadius(Eigen::MatrixXd thisVertices, Eigen::Matrix<double, 1,3> thisCenter) {
    // compute the distance between each point and the center
    int numberPoint = thisVertices.rows() ;
    Eigen::Matrix3d centerArray = thisCenter.replicate(numberPoint,1) ;
    Eigen::Matrix<double, numberPoint,1> distances = (thisVertices-centerArray).rowwise().norm() ;

    // compute the mean and return it 
    double meanRadius = distances.colwise().mean() ;
    return meanRadius ;
}

/********* ============= Functions to handle PLANE =============== *********/
    bool computePlane(Eigen::Matrix3i sample_idx,
                      Eigen::Matrix3d variance,
                      DecoratedCloud &cloud,
                      CloudPrimitive &primitives,
                      double thresh, double alpha) {
        Eigen::MatrixXd V = cloud.getVertices(), N = cloud.getNormals();
        const int cloudSize = V.rows();
        const int nSamples = sample_idx.rows();

        // ---- Retrieve the N vertices and their normals ----
        Eigen::Matrix3d thisVertex, thisNormal;
        for (int i = 0; i < 3; i++) {
            thisVertex.row(i) = V.row(sample_idx(i, 1));
            thisNormal.row(i) = N.row(sample_idx(i, 1));
        }

        if (isPlane(thisVertex, thisNormal, thresh, alpha)) {

            // ---- Create a new plane and compute its score ----
            Eigen::RowVector3d planeNormal = computeNormal(thisVertex, thisNormal);
            Eigen::RowVector3d planeRefPoint = V.colwise().mean();
            Plane newPlane = Plane(planeRefPoint, planeNormal);
            newPlane.computeScore(variance, cloud, thresh, alpha);

            // ---- Store it in the cloudPrimitive ----
            primitives.addPrimitive(newPlane);
            return 1;
        } else {
            return 0;
        }
    }

    /*** ----------- isPlane() ----------- ***/
    bool isPlane(Eigen::MatrixXd V, Eigen::MatrixXd N, double T, double alpha) {
        Eigen::RowVector3d planeNormal = computeNormal(V, N);
        bool isPlane = true;
        for (int i = 0; i < N.rows(); i++) {
            if (N.row(i).dot(planeNormal) < T) isPlane = false;
        }
        return isPlane;
    }

    /*** ----------- computeNormal() ----------- ***/
    Eigen::RowVector3d computeNormal(Eigen::MatrixXd V, Eigen::MatrixXd _N) {
        Eigen::RowVector3d N = Eigen::RowVector3d::Zero();
        for (int i = 0; i < V.rows() - 2; i++) {
            Eigen::RowVector3d P01 = V.row(1 + i) - V.row(i);
            Eigen::RowVector3d P02 = V.row(2 + i) - V.row(i);
            N += P02.cross(P01) / (V.rows() - 2);
        }
        if (_N.row(0).dot(N) < 0) N = -N;
        return N;
    }



    /********* ============= Functions to handle the final cloud =============== *********/
    DecoratedCloud fuse(CloudPrimitive& best_primitives,
                        CloudManager& clouds,
                        double T_rad,  // Radius   Distance Threshold (Sphere)
                        double T_cent, // Center   Distance Threshold (Sphere)
                        double T_norm, // Normals  Distance Threshold (Plane)
                        double T_refPt // RefPoint Distance Threshold (Plane)
    ){

        std::vector<std::pair<int,int>> fuses;
        const int n = best_primitives.getCloudSize();

        // === Consider every pair of primitive for merging ===
        Primitive first_prim, second_prim;

        for(int i=0; i<n-1; i++){ // First Primitive
            first_prim = best_primitives.getPrimitive(i);
            for(int j=i+1; j<n; j++){ // Second Primitive
                second_prim = best_primitives.getPrimitive(j);

                // --- Both primitives have the same type ---
                if(first_prim.getType()==second_prim.getType()){

                    // ---- They are both Spheres ---
                    if(first_prim.getType()==1){
                        double d1 = (first_prim.getCenter() - second_prim.getCenter()).norm();
                        double d2 = first_prim.getRadius() - second_prim.getRadius();
                        if(d1 < T_cent && d2 < T_rad){
                            fuses.push_back(std::make_pair(i,j));
                        }
                    }// ---- They are both Planes ---
                    else{
                        double d1 = (first_prim.getNormal().dot(second_prim.getNormal()));
                        double d2 = (first_prim.getRefPoint() - second_prim.getRefPoint()).dot(first_prim.getNormal());
                        if(d1 < T_norm && d2 < T_refPt){
                            fuses.push_back(std::make_pair(i,j));
                        }
                    }
                }
            }
        }

        // === Merge all the primitives in the cloudManager ===
        for(int i=0; i<fuses.size(); i++){

            int pr_1 = fuses.at(i).first;
            int pr_2 = fuses.at(i).second;

            // Merge both primitives
            Eigen::MatrixXd V1 = clouds.getCloud(pr_1).getVertices();
            Eigen::MatrixXd V2 = clouds.getCloud(pr_2).getVertices();
            Eigen::MatrixXd N1 = clouds.getCloud(pr_1).getNormals();
            Eigen::MatrixXd N2 = clouds.getCloud(pr_2).getNormals();
            Eigen::MatrixXi F1 = clouds.getCloud(pr_1).getFaces();
            Eigen::MatrixXi F2 = clouds.getCloud(pr_2).getFaces();
            V1 << V1, V2;
            F1 << F1, F2;
            N1 << N1, N2;
            // Update corresponding Meshes
            clouds.getCloud(pr_1).setVertices(V1);
            clouds.getCloud(pr_1).setFaces(F1);
            clouds.getCloud(pr_1).setNormals(N1);
        }
        for(int i=0; i<fuses.size(); i++) { // Delete redundant clouds
            clouds.deleteCloud(fuses.at(i).second);
        }

        // === Build new cloud with color attributes ===
        DecoratedCloud newCloud;
        Eigen::MatrixXd V, C, N;
        Eigen::MatrixXi F;

        for(int i=0; i< clouds.getCloudSize(); i++){
            V << V, clouds.getCloud(i).getVertices();
            N << N, clouds.getCloud(i).getNormals();
            F << F, clouds.getCloud(i).getFaces();
            C << C, Eigen::RowVector3d(std::rand()/double(RAND_MAX),
                                       std::rand()/double(RAND_MAX),
                                       std::rand()/double(RAND_MAX)).replicate(
                    clouds.getCloud(i).getFaces().rows(), 1);
        }


        { // ---- Update Cloud ---
            newCloud.setVertices(V);
            newCloud.setColors(C);
            newCloud.setFaces(F);
            newCloud.setNormals(N);
        }
        return newCloud;

    };


    bool cleanCloud(DecoratedCloud& oldCloud, DecoratedCloud& newCloud, Eigen::Matrix3i inliers_idx){
            // New Cloud = Old cloud without the detected Inliers


    };


}


