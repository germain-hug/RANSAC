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
       // Compute the matrix minus the min of each coordinate
        Eigen::MatrixXd centered = V.rowwise() - V.colwise().mean();

        // Compute the covariance matrix
        Eigen::Matrix3d cov = (centered.adjoint() * centered) / double(V.rows());
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
        Eigen::MatrixXd V = cloud.getVertices(), N = cloud.getNormals();

        const int cloudSize = V.rows();
        const int nSamples = sample_idx.rows();


        // ---- Retrieve the N vertices and their normals ----
        Eigen::Matrix3d thisVertex, thisNormal;
        for (int i = 0; i < 3; i++) {
            thisVertex.row(i) = V.row(sample_idx(i, 0));
            thisNormal.row(i) = N.row(sample_idx(i, 0));
        }

        if (isPlane(thisVertex, thisNormal, thresh, alpha)) {
            // ---- Create a new plane and compute its score ----
            Eigen::Matrix<double, 1,3> planeNormal = computeNormal(thisVertex, thisNormal.row(0));
            Eigen::Matrix<double, 1,3> planeRefPoint = thisVertex.colwise().mean();

            Primitive* newPlane = new Plane(planeRefPoint, planeNormal);
            newPlane->computeScore(variance, cloud, thresh, alpha);
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
            if (std::abs(N.row(i).dot(planeNormal)) < alpha) isPlane = false;
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
            N += P02.cross(P01) / (V.rows() - 2);
        }

        // Check for normal orientation
        //if (_N.dot(N) < 0) N = -N;
        // Normalize
        if(N(0,0)!=0 && N(0,1)!=0 && N(0,2)!=0) N = N.normalized();
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

        // initialization 
        int nbCloudInitial = clouds.getCloudSize();
        Eigen::MatrixXi visited = Eigen::MatrixXi::Zero(1, nbCloudInitial) ;
        double d1,d2 ;
        int current_label = 0 ;
        // === Consider every pair of primitive for merging ===
        Primitive* first_prim, *second_prim;

        for(int i=0; i < nbCloudInitial; i++){ // First Primitive
            if (visited(0,i) == 0) {
                current_label = current_label+1 ;
                visited(0,i) = current_label ;

                first_prim = best_primitives.getPrimitive(i);

                for(int j=i+1; j<nbCloudInitial; j++){ // Second Primitive
                second_prim = best_primitives.getPrimitive(j);

                // --- Both primitives have the same type ---
                if(first_prim->getType()==second_prim->getType()){

                    // ---- They are both Spheres ---
                    if(first_prim->getType()==1){
                        d1 = (static_cast<Sphere*>(first_prim)->getCenter() - static_cast<Sphere*>(second_prim)->getCenter()).norm();
                        d2 = std::abs(static_cast<Sphere*>(first_prim)->getRadius() - static_cast<Sphere*>(second_prim)->getRadius());

                        if(d1 < T_cent && d2 < T_rad){
                            visited(0,j) = current_label ;
                        }
                    }// ---- They are both Planes ---
                    else{
                        double d1 = std::abs((static_cast<Plane*>(first_prim)->getNormal().dot(static_cast<Plane*>(second_prim)->getNormal())));
                        double d2 = std::abs((static_cast<Plane*>(first_prim)->getRefPoint() - static_cast<Plane*>(second_prim)->getRefPoint()).dot(static_cast<Plane*>(second_prim)->getNormal()));
                        double d3 = std::abs((static_cast<Plane*>(first_prim)->getRefPoint() - static_cast<Plane*>(second_prim)->getRefPoint()).dot(static_cast<Plane*>(first_prim)->getNormal()));

                        if(d1 > T_norm && (d2 < T_refPt || d3 < T_refPt)){
                            visited(0,j) = current_label ;
                        }
                    }
                }
            }
            }
        }
        // need to fuse the meshes

        int numberOfVertices = 0 ;
        // === Merge all the primitives in the cloudManager ===
        for(int thisLabel=1; thisLabel<=current_label; thisLabel++) {

            // determine the size of the new vertices 
            for(int j=0; j< nbCloudInitial; j++) {
                if (visited(0,j) == thisLabel) {
                    numberOfVertices += clouds.getCloud(j).getVertices().rows() ;
                }
            }

            // create the new matrices 
            Eigen::MatrixXd V(numberOfVertices, 3) ;
            Eigen::MatrixXd N(numberOfVertices, 3) ;
            Eigen::MatrixXd C(numberOfVertices,3) ;

            int indiceStart = 0 ;
            int nbVertCloud = 0;

            for(int j=0; j< nbCloudInitial; j++){
                if (visited(0,j) == thisLabel) {
                // number of vertex for this cloud
                nbVertCloud = clouds.getCloud(j).getVertices().rows() ;

                // fill each block
                V.block(indiceStart,0,nbVertCloud,3) = clouds.getCloud(j).getVertices();

                N.block(indiceStart,0,nbVertCloud,3) = clouds.getCloud(j).getNormals();

                C.block(indiceStart,0,nbVertCloud,3) = Eigen::RowVector3d(std::rand()/double(RAND_MAX),
                                       std::rand()/double(RAND_MAX),
                                       std::rand()/double(RAND_MAX)).replicate(
                    clouds.getCloud(j).getVertices().rows(), 1);

                // update the indice to start filling
                indiceStart += nbVertCloud  ;
                }
            }
            clouds.setCloud(DecoratedCloud(V,N,C), thisLabel-1) ;
            numberOfVertices = 0;
        }
        if(clouds.getCloudSize() - current_label > 0) {
            clouds.deleteCloudFromIndex(current_label+1);
        }

    }


    // take a cloudManager and gather all the cloud in one 
    DecoratedCloud* gatherClouds(CloudManager& cloudManager, int colorExit) {
        // === Build new cloud with color attributes ===
        int numberOfCloud = cloudManager.getCloudSize() ;
        int numberOfVertices =0 ;

        // determine the size of the new cloud 
        for(int i=0; i< numberOfCloud; i++){
            numberOfVertices += cloudManager.getCloud(i).getVertices().rows() ;
        }

        // create the matrix to store the result
        Eigen::MatrixXd V(numberOfVertices,3) ;
        Eigen::MatrixXd C(numberOfVertices,3) ;
        Eigen::MatrixXd N(numberOfVertices,3);

        int indiceStart = 0 ;
        int nbVertCloud ;

        for(int i=0; i< numberOfCloud; i++){
            // number of vertex for this cloud 
            nbVertCloud = cloudManager.getCloud(i).getVertices().rows() ;

            // fill each block 
            V.block(indiceStart,0,nbVertCloud,3) = cloudManager.getCloud(i).getVertices();

            N.block(indiceStart,0,nbVertCloud,3) = cloudManager.getCloud(i).getNormals();

            if (colorExit==1) {
                C.block(indiceStart,0,nbVertCloud,3) = cloudManager.getCloud(i).getColors();
            }
            else { 
                C.block(indiceStart,0,nbVertCloud,3) = Eigen::RowVector3d(std::rand()/double(RAND_MAX),
                                       std::rand()/double(RAND_MAX),
                                       std::rand()/double(RAND_MAX)).replicate(
                    cloudManager.getCloud(i).getVertices().rows(), 1);
            }

            // update the indice to start filling 
            indiceStart += nbVertCloud  ;
        }

        // ---- Create new Cloud ---
        DecoratedCloud* newCloud = new DecoratedCloud(V,N,C) ;
        return newCloud ;
    }


    void cleanCloud(DecoratedCloud& cloudRansac, CloudManager& cloudManager, Eigen::MatrixXi inliers_idx){
        // ---- We remove inliers from cloudRansac -----
        int n_inliers = inliers_idx.rows();

        if(n_inliers > 0) {
            int n_cloud = cloudRansac.getVertices().rows();
            Eigen::MatrixXd V_in(n_cloud, 3);
            Eigen::MatrixXd V_out(n_cloud, 3) ;
            Eigen::MatrixXd N_in(n_cloud, 3) ;
            Eigen::MatrixXd N_out(n_cloud, 3);

            int inliers_valid = 0, outliers_valid = 0;

            // ---- For every vertex, search if is an inlier ---
            for(int i=0; i<n_cloud; i++){
                bool isValid = true;
                for(int j=0; j<n_inliers; j++){
                    if(inliers_idx(j,0)==i){isValid = false;}
                }

                // Vertex is valid, save it to the current cloud
                if(isValid){
                    V_in.row(inliers_valid) =  cloudRansac.getVertices().row(i);
                    N_in.row(inliers_valid) =  cloudRansac.getNormals().row(i);
                    inliers_valid++;
                } else{
                    // Vertex is non valid, add it to cloud of inliers
                    V_out.row(outliers_valid) = cloudRansac.getVertices().row(i);
                    N_out.row(outliers_valid) = cloudRansac.getNormals().row(i);
                    outliers_valid++;
                }
             }
            Eigen::MatrixXd N_out_top = N_out.topRows(outliers_valid-1);

            cloudManager.addCloud(DecoratedCloud(V_out.topRows(outliers_valid-1),N_out_top)); // Store cloud of inliers
            cloudRansac.setVertices(V_in.topRows(inliers_valid-1)); // Keep cloud deprived from inliers
            cloudRansac.setNormals(N_in.topRows(inliers_valid-1));
        }

    }




    /// ---- RECONSTRUCTION ----
    void reconstruct(CloudPrimitive& best_primitives, DecoratedCloud& cloud, int nbSamples, double T, double alpha) {
        const int n = best_primitives.getCloudSize();
        for (int i = 0; i < n; i++) { // For every primitive

            if (best_primitives.getPrimitive(i)->getType() == 2) { // If we found a plane
                Primitive *this_prim = best_primitives.getPrimitive(i);
                Eigen::MatrixXi inliers_idx = this_prim->computeInliers(cloud, T, alpha); // Retrieve the inliers
                std::cout << "cloud Size " << cloud.getVertices().rows() << std::endl;

                sampleFromPrimitive(cloud, inliers_idx, *this_prim, nbSamples);
                std::cout << "cloud Size " << cloud.getVertices().rows() << std::endl;

            }

        }
    }



    void sampleFromPrimitive(DecoratedCloud& cloud, Eigen::MatrixXi inliers_idx, Primitive plane, int nbSample) {
        /* Compute the optimal number of points for this plane, given the cloud variance and the inliers */

        Eigen::MatrixXd V = cloud.getVertices();
        Eigen::MatrixXd C = cloud.getColors();
        std::cout << "V.rows()" << V.rows() << std::endl;


        // --- Find an orthonormal basis of the plane ---
        Eigen::Matrix<double, Eigen::Dynamic, 3> N(1, 3); N = plane.getNormal().normalized();
        Eigen::Matrix<double, Eigen::Dynamic, 3> u(1, 3);

        // Compute u so that it's orthogonal to N
        if(N(0,1)!=0 || N(0,0)!=0)
            u << -N(0,1), N(0,0), 0.0;
        else if(N(0,2)!=0 || N(0,0)!=0)
            u << -N(0,2), 0.0, N(0,0);
        else u << 0.0, -N(0,2), N(0,1);
        if(u.norm()!=0.0) u.normalize();

        // Compute v to make a basis (u,v,N)
        Eigen::MatrixXd v = u.row(0).cross(N.row(0));
        if(v.norm()!=0.0) v.normalize();


        // --- Retrieve 3D Inliers and project on (u,v) basis ---
        const int n = inliers_idx.rows();
        Eigen::MatrixXd inliers2D(n, 2);
        Eigen::MatrixXd this_vertex ;
        int idx;

        for(int i = 0; i<n; i++){
            idx = inliers_idx(i,0);
            this_vertex = cloud.getVertices().row(idx);
            inliers2D(i,0) = u.row(0).dot(this_vertex.row(0));
            inliers2D(i,1) = v.row(0).dot(this_vertex.row(0));
        }

        // --- Compute Bounding Box Area ---
        double x_min = inliers2D.col(0).minCoeff();
        double x_max = inliers2D.col(0).maxCoeff();
        double y_min = inliers2D.col(1).minCoeff();
        double y_max = inliers2D.col(1).maxCoeff();



        // --- Generate new 2D samples in (u,v) coordinates ---
        srand(time(NULL));
        double this_x, this_y, f_x, f_y;

        for(int i = 0; i < nbSample; i++) {
            f_x = (double) rand() / RAND_MAX;
            f_y = (double) rand() / RAND_MAX;

            this_x = x_min + f_x * (x_max - x_min);
            this_y = y_min + f_y * (y_max - y_min);

            bool isValid = true;
            for (int k = 0; k < inliers2D.rows(); k++) {
                if (sqrt(pow(inliers2D(k, 0) - this_x, 2.0) + pow(inliers2D(k, 1) - this_y, 2.0)) < 0.05) {
                    isValid = false;
                    break;
                }
            }


            // We found a valid sample, compute its world coordinates
            if (isValid) {

                double x = this_x*u(0,0) +  this_y*v(0,0);
                double y = this_x*u(0,1) +  this_y*v(0,1);
                double z = this_x*u(0,2) +  this_y*v(0,2);
                std::cout << "x " << x << " y " << y << " z " << z << std::endl;
                Eigen::MatrixXd new_vertex(1,3); new_vertex << x, y, z;
                Eigen::MatrixXd new_color(1,3); new_color << 1.0, 0.0, 0.0;
                std::cout << "V.rows()" << V.rows() << std::endl;

                V << V, new_vertex;
                C << C, new_color;
                std::cout << "V.rows()" << V.rows() << std::endl;

            }
        }
        std::cout << "V.rows()" << V.rows() << std::endl;
        cloud.setVertices(V);
        cloud.setColors(C);
    }

}


