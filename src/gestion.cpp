#include "acq/impl/gestion.hpp"

namespace acq {

//  ****** ============ Helper Functions  =============== ******* 
// Sample the point cloud 
Eigen::MatrixXi sample(int cloudSize) {
    // we want to have 3 points sampled 
    Eigen::Matrix<int, 3,1> sampleInd(3,1) ;
    int newIndex ;
    // add a random indices between 0 and sizeMatrix in a numberPoint sized vector 
    for (int i=0; i<3; i++) {
        bool isUnique = false;
        while(!isUnique){
            newIndex = rand() % (cloudSize + 1) ;
            isUnique = true;
            for(int j=0; j<i; j++) {
                if(sampleInd(j)==newIndex){
                    isUnique = false;
                }
            }
        }
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
int isSphere(Eigen::Matrix3d vertices, Eigen::Matrix3d normals, double threshold, double alpha) {
    // estimate the center and the radius using 2 points 
    Eigen::Matrix<double, 1,3> thisCenter = computerCenter(vertices.topRows(2), normals.topRows(2)) ;
    double estimatedRadius = computerRadius(vertices.topRows(2), thisCenter) ;

    // compute the estimated normal for the least point  
    Eigen::Matrix<double, 1,3> estimatedNormal = vertices.row(2) - thisCenter ;
    estimatedNormal.normalize() ;

    // test for the radius 
    double test1 = computerRadius(vertices.row(2), thisCenter) - estimatedRadius ;
    double test2 = estimatedNormal.dot(normals.row(2).normalized()) ;

    int thisReturn = 0 ;

    if (std::abs(test1) < threshold ) {
        if ( test2 > alpha ) {
            // if the 2 test are true, the 3 points form a sphere  
            thisReturn = 1 ; 
        }
    }
    return thisReturn ;
}

// if the 3 points create a sphere, we add it to the primitives 
void computeSphere(Eigen::Matrix<int, 3,1> sample_idx, Eigen::Matrix3d variance, DecoratedCloud& cloud, CloudPrimitive& primitives, double threshold, double alpha) {
    std::cout << "Enter compute Sphere " << std::endl ;

    Eigen::MatrixXd vertices = cloud.getVertices() ;
    Eigen::MatrixXd normals = cloud.getNormals() ;
    int cloudSize = vertices.rows() ;

    Eigen::Matrix3d thisSampledVertices ;
    Eigen::Matrix3d thisSampledNormals ;

    // extract the 3 points sampled by the indices
    for(int i =0 ; i< 3; i++) {
        thisSampledVertices.row(i) = vertices.row(sample_idx(i,0)) ;
        thisSampledNormals.row(i) = normals.row(sample_idx(i,0)) ;
    }
    // test if it's a sphere
    int is_sphere = isSphere(thisSampledVertices, thisSampledNormals, threshold, alpha) ;

    if (is_sphere==1) {
        // compute the attribut for the object 
        Eigen::Matrix<double, 1,3> thisCenter = computerCenter(thisSampledVertices, thisSampledNormals) ;
        double thisRadius = computerRadius(thisSampledVertices, thisCenter) ;

        // create the object and compute its score 
        Primitive* thisSphere = new Sphere(thisRadius, thisCenter) ;
        thisSphere->computeScore(variance, cloud, threshold, alpha) ;

        // set the type 
        thisSphere->setType(1) ;
        
        // store it in the cloud primitive 
        primitives.addPrimitive(thisSphere) ;
    }
}

// compute the center of a shpere by finding the better intersection possible using least square computation
Eigen::Matrix<double, 1,3> computerCenter(Eigen::MatrixXd vertices, Eigen::MatrixXd normals) {
    Eigen::Matrix3d R = Eigen::Matrix3d::Zero(3,3) ;
    Eigen::Matrix<double, 3,1> q = Eigen::MatrixXd::Zero(3,1) ;
    Eigen::Matrix3d I = Eigen::Matrix3d::Identity(3,3) ;

    int numberPoint = vertices.rows() ;
    Eigen::Matrix<double, 1,3> thisNormal ;
    Eigen::Matrix<double, 1,3> thisPosition ;

    // fill the system 
    for (int i = 0; i< numberPoint; i++) {
        thisNormal = normals.row(i) ;
        thisPosition = vertices.row(i) ;

        R += I - thisNormal.transpose()*thisNormal ; 

        q += (I - thisNormal.transpose()*thisNormal) * thisPosition.transpose() ;
    }

    // solve the system using least jacobi decomposition 
    Eigen::Matrix<double, 3,1> thisCenter ;


    thisCenter = R.jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(q) ;

    return thisCenter.transpose() ;
}

double computerRadius(Eigen::MatrixXd thisVertices, Eigen::Matrix<double, 1,3> thisCenter) {
    // compute the distance between each point and the center
    int numberPoint = thisVertices.rows() ;
    Eigen::MatrixXd centerArray = thisCenter.replicate(numberPoint,1) ;
    Eigen::MatrixXd distances(numberPoint,1) ; 
    distances = (thisVertices-centerArray).rowwise().norm() ;

    // compute the mean and return it 
    double meanRadius = distances.mean() ;
    return meanRadius ;
}

/********* ============= Functions to handle PLANE =============== *********/
    void computePlane(Eigen::Matrix<int, 3,1> sample_idx,
                      Eigen::Matrix3d variance,
                      DecoratedCloud &cloud,
                      CloudPrimitive &primitives,
                      double thresh, double alpha) {        
        Eigen::MatrixXd V = cloud.getVertices(),
        N = cloud.getNormals();
        const int cloudSize = V.rows();
        const int nSamples = sample_idx.rows();


        // ---- Retrieve the N vertices and their normals ----
        Eigen::Matrix3d thisVertex, thisNormal;
        for (int i = 0; i < 3; i++) {
            thisVertex.row(i) = V.row(sample_idx(i, 0));
            thisNormal.row(i) = N.row(sample_idx(i, 0));
        }

        std::cout << "thisVertex(): " << thisVertex << std::endl;

        if (isPlane(thisVertex, thisNormal, thresh, alpha)) {
            std::cout << "Plane detected" << std::endl ;

            // ---- Create a new plane and compute its score ----
            Eigen::Matrix<double, 1,3> planeNormal = computeNormal(thisVertex, thisNormal.row(0));

            Eigen::Matrix<double, 1,3> planeRefPoint = V.colwise().mean();
            std::cout << "Ref point computed : "<< planeRefPoint << std::endl ;

            Primitive* newPlane = new Plane(planeRefPoint, planeNormal);
            std::cout << "Plane created "<< std::endl ;

            newPlane->computeScore(variance, cloud, thresh, alpha);
            std::cout << "Score computed "<< std::endl ;

            newPlane->setType(2) ;

            // ---- Store it in the cloudPrimitive ----
            primitives.addPrimitive(newPlane);
        } 
    }

    /*** ----------- isPlane() ----------- ***/
    bool isPlane(Eigen::MatrixXd V, Eigen::MatrixXd N, double T, double alpha) {
        Eigen::Matrix<double, 1,3> planeNormal = computeNormal(V, N.row(0));
        bool isPlane = true;
        for (int i = 0; i < N.rows(); i++) {
            if (N.row(i).dot(planeNormal) < T) isPlane = false;
        }
        return isPlane;
    }

    /*** ----------- computeNormal() ----------- ***/
    Eigen::Matrix<double, 1,3> computeNormal(Eigen::MatrixXd V, Eigen::Matrix<double, 1,3> _N) {
        Eigen::Matrix<double, 1,3> N; N << 0.0,0.0,0.0;
        Eigen::Matrix<double, 1,3> P01, P02 ;

        for (int i = 0; i < V.rows() - 2; i++) {
            P01 = V.row(1 + i) - V.row(i);
            P02 = V.row(2 + i) - V.row(i);
            std::cout << "P01: " << P01 << " | P02: " << P02 <<  std::endl;
            std::cout << " V: " << V << std::endl;
            N += P02.cross(P01) / (V.rows() - 2);
        }

        // Check for normal orientation
        if (_N.dot(N) < 0) N = -N;
        // Normalize
        if(N(0,0)!=0 && N(0,1)!=0 && N(0,2)!=0) N = N.normalized();

        std::cout << "Normal Computed : " << N << std::endl ;
        return N;
    }



    /********* ============= Functions to handle the final cloud =============== *********/
    // fuse the cloud with the same primitives 
    void fuse(CloudPrimitive& best_primitives,
              CloudManager& clouds,
              double T_rad,  // Radius   Distance Threshold (Sphere)
              double T_cent, // Center   Distance Threshold (Sphere)
              double T_norm, // Normals  Distance Threshold (Plane)
              double T_refPt // RefPoint Distance Threshold (Plane)
    ) {

        std::vector<std::pair<int,int>> fuses;
        const int n = best_primitives.getCloudSize();

        // === Consider every pair of primitive for merging ===
        Primitive* first_prim, *second_prim;

        for(int i=0; i<n-1; i++){ // First Primitive
            first_prim = best_primitives.getPrimitive(i);
            for(int j=i+1; j<n; j++){ // Second Primitive
                second_prim = best_primitives.getPrimitive(j);

                // --- Both primitives have the same type ---
                if(first_prim->getType()==second_prim->getType()){

                    // ---- They are both Spheres ---
                    if(first_prim->getType()==1){
                        double d1 = (first_prim->getCenter() - second_prim->getCenter()).norm();
                        double d2 = std::abs(first_prim->getRadius() - second_prim->getRadius());
                        if(d1 < T_cent && d2 < T_rad){
                            fuses.push_back(std::make_pair(i,j));
                        }
                    }// ---- They are both Planes ---
                    else{
                        double d1 = (first_prim->getNormal().dot(second_prim->getNormal()));
                        double d2 = (first_prim->getRefPoint() - second_prim->getRefPoint()).dot(first_prim->getNormal());
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
            V1 << V1, V2;
            N1 << N1, N2;

            // Update corresponding Meshes
            clouds.getCloud(pr_1).setVertices(V1);
            clouds.getCloud(pr_1).setNormals(N1);
        }
        for(int i=0; i<fuses.size(); i++) { // Delete redundant clouds
            clouds.deleteCloud(fuses.at(i).second);
        }
    }

    // take a cloudManager and gather all the cloud in one 
    DecoratedCloud gatherClouds(CloudManager& cloudManager) {
        // === Build new cloud with color attributes ===
        int numberOfCloud = cloudManager.getCloudSize() ;
        int numberOfVertices =0 ;

        std::cout << "enter in gathercloud " << numberOfCloud << " to merge " << std::endl ;

        // determine the size of the new cloud 
        for(int i=0; i< numberOfCloud; i++){
            numberOfVertices += cloudManager.getCloud(i).getVertices().rows() ;
        }

        // create the matrix to store the result
        Eigen::MatrixXd V(numberOfVertices,3) ;
        Eigen::MatrixXd C(numberOfVertices,3) ;
        Eigen::MatrixXd N(numberOfVertices,3);

        std::cout << "test normals : " << cloudManager.getCloud(0).getNormals() << std::endl ;

        int indiceStart = 0 ;
        int nbVertCloud ;

        for(int i=0; i< numberOfCloud; i++){
            // number of vertex for this cloud 
            nbVertCloud = cloudManager.getCloud(i).getVertices().rows() ;

            // fill each block 
            V.block(indiceStart,0,nbVertCloud,3) = cloudManager.getCloud(i).getVertices();

        std::cout << "Vertex OK" << std::endl ;

            N.block(indiceStart,0,nbVertCloud,3) = cloudManager.getCloud(i).getNormals();

        std::cout << "Normals OK" << std::endl ;

            C.block(indiceStart,0,nbVertCloud,3) = Eigen::RowVector3d(std::rand()/double(RAND_MAX),
                                       std::rand()/double(RAND_MAX),
                                       std::rand()/double(RAND_MAX)).replicate(
                    cloudManager.getCloud(i).getVertices().rows(), 1);

        std::cout << "colors OK" << std::endl ;

            // update the indice to start filling 
            indiceStart += nbVertCloud  ;
        }

        std::cout << "Vertex test : " << V << std::endl ;

        // ---- Create new Cloud ---
        DecoratedCloud newCloud = DecoratedCloud(V,N,C) ;
        return newCloud;
    }


    void cleanCloud(DecoratedCloud& cloudRansac, CloudManager& cloudManager, Eigen::Matrix3i inliers_idx){
        std::cout << "Enter Clean Cloud" << std::endl;
        // ---- We remove inliers from cloudRansac -----
        const int n_inliers = inliers_idx.rows();
        if(n_inliers > 0) {
            const int n_cloud = cloudRansac.getVertices().rows();
            Eigen::MatrixXd V_in(n_cloud, 3), V_out(n_cloud, 3);
            int inliers_valid = 0, outliers_valid = 0;

            // ---- For every vertex, search if is an inlier ---
            for(int i=0; i<n_cloud; i++){
                bool isValid = true;
                for(int j=0; j<n_inliers; j++){
                    if(inliers_idx(j,0)==i){isValid = false; break;}
                }

                // Vertex is valid, save it to the current cloud
                if(isValid){
                    V_in << V_in, cloudRansac.getVertices().row(i);
                    inliers_valid++;
                } else{
                    // Vertex is non valid, add it to cloud of inliers
                    V_out << V_out, cloudRansac.getVertices().row(i);
                    outliers_valid++;
                }
             }
            cloudManager.addCloud(DecoratedCloud(V_out.topRows(inliers_valid-1))); // Store cloud of inliers
            cloudRansac.setVertices(V_in.topRows(outliers_valid-1)); // Keep cloud deprived from inliers
        }
        std::cout << "Exit Clean Cloud" << std::endl;
    }



}


