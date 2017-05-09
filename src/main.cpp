#include "acq/normalEstimation.h"
#include "acq/ransac.h"

#include "nanogui/formhelper.h"
#include "nanogui/screen.h"

#include "igl/readOBJ.h"
#include "igl/viewer/Viewer.h"

#include <iostream>

namespace acq {

/** \brief                      Re-estimate normals of cloud \p V fitting planes
 *                              to the \p kNeighbours nearest neighbours of each point.
 * \param[in ] kNeighbours      How many neighbours to use (Typiclaly: 5..15)
 * \param[in ] vertices         Input pointcloud. Nx3, where N is the number of points.
 * \param[in ] maxNeighbourDist Maximum distance between vertex and neighbour.
 * \param[out] viewer           The viewer to show the normals at.
 * \return                      The estimated normals, Nx3.
 */
NormalsT
recalcNormals(
    int                 const  kNeighbours,
    CloudT              const& vertices,
    float               const  maxNeighbourDist
) {
    NeighboursT const neighbours =
        calculateCloudNeighbours(
            /* [in]        cloud: */ vertices,
            /* [in] k-neighbours: */ kNeighbours,
            /* [in]      maxDist: */ maxNeighbourDist
        );

    // Estimate normals for points in cloud vertices
    NormalsT normals =
        calculateCloudNormals(
            /* [in]               Cloud: */ vertices,
            /* [in] Lists of neighbours: */ neighbours
        );

    return normals;
} //...recalcNormals()

void setViewerNormals(
    igl::viewer::Viewer      & viewer,
    CloudT              const& vertices,
    NormalsT            const& normals
) {
    // [Optional] Set viewer face normals for shading
    //viewer.data.set_normals(normals);

    // Clear visualized lines (see Viewer.clear())
    viewer.data.lines = Eigen::MatrixXd(0, 9);

    // Add normals to viewer
    viewer.data.add_edges(
        /* [in] Edge starting points: */ vertices,
        /* [in]       Edge endpoints: */ vertices + normals * 0.01, // scale normals to 1% length
        /* [in]               Colors: */ Eigen::Vector3d::Zero()
    );
}

} //...ns acq

int main(int argc, char *argv[]) {

    // How many neighbours to use for normal estimation, shown on GUI.
    int kNeighbours = 10;
    // Maximum distance between vertices to be considered neighbours (FLANN mode)
    float maxNeighbourDist = 0.15; //TODO: set to average vertex distance upon read

    // Dummy enum to demo GUI
    enum Orientation { Up=0, Down, Left, Right } dir = Up;
    // Dummy variable to demo GUI
    bool boolVariable = true;
    // Dummy variable to demo GUI
    float floatVariable = 0.1f;

    // ********* VARIABLES FOR THE ALGORITHM  ********* 
    int nbIteration = 5 ; 
    int samplePerIt = 5 ;
    double thresh = 0.002 ;
    double alpha = 0.60 ;
    double thresh_best = 50.0 ;
    // will store the current primitives and the point cloud per primitives
    acq::CloudPrimitive best_primitives ;
    acq::CloudManager cloudManagerParts ;

    // deals with several meshes 
    enum MeshType { test1=0, test2, test3} typeMesh = test1 ;
    //************************************
    
    // Load a mesh in OFF format
    std::string meshPath = "../models/sphere.off";
    if (argc > 1) {
        meshPath = std::string(argv[1]);
        if (meshPath.find(".obj") == std::string::npos) {
            std::cerr << "Only ready for  OBJ files for now...\n";
            return EXIT_FAILURE;
        }
    } else {
        std::cout << "Usage: iglFrameWork <path-to-off-mesh.obj>." << "\n";
    }

    // Visualize the mesh in a viewer
    igl::viewer::Viewer viewer;
    viewer.core.show_lines = false;
    viewer.core.show_overlay = false;

    // Store cloud so we can store normals later
    acq::CloudManager cloudManagerOldMesh;
    // Read mesh from meshPath
    {
        Eigen::MatrixXd V, TC, N;
        Eigen::MatrixXi F, FTC, FN;
        igl::readOFF(meshPath, V, F);

        //if(FTC.size() == 0) FTC = F;

        if (V.rows() <= 0) {
            std::cerr << "Could not read mesh at " << meshPath
                      << "...exiting...\n";
            return EXIT_FAILURE;
        }


        // Store read vertices and faces
        //N.rowwise().normalize();

        cloudManagerOldMesh.addCloud(acq::DecoratedCloud(V, F));

        // Set Normals from OBJ file
        //cloudManagerOldMesh.getCloud(typeMesh).setNormals(N);
        //std::cout << N.size() << std::endl;

        // Update viewer
        //acq::setViewerNormals(viewer, cloudManagerOldMesh.getCloud(typeMesh).getVertices(), N);
    
        // set the mesh 
        viewer.data.clear() ;

        // Show mesh
        viewer.data.set_mesh(
            cloudManagerOldMesh.getCloud(typeMesh).getVertices(),
            cloudManagerOldMesh.getCloud(typeMesh).getFaces()
        );
    
    }

    // Extend viewer menu using a lambda function
    viewer.callback_init =
        [
            &cloudManagerOldMesh, &kNeighbours, &maxNeighbourDist,
            &floatVariable, &boolVariable, &dir, &nbIteration, &samplePerIt, 
            &best_primitives, &cloudManagerParts, &thresh, &alpha, &thresh_best, 
            &typeMesh 
        ] (igl::viewer::Viewer& viewer)
    {
        // Add an additional menu window
        viewer.ngui->addWindow(Eigen::Vector2i(900,10), "Acquisition3D");

        // ***** TODO : add different meshes *****
        viewer.ngui->addGroup("Choose your mesh");

        viewer.ngui->addVariable<MeshType>("Which mesh do you want ?",typeMesh)->setItems(
            {"Mesh1","Mesh2","Mesh3"}
        );

        viewer.ngui->addButton("Show the original mesh",
                               [&]() {
            viewer.data.clear() ;

            // Show mesh
            viewer.data.set_mesh(
                             cloudManagerOldMesh.getCloud(typeMesh).getVertices(),
                             cloudManagerOldMesh.getCloud(typeMesh).getFaces()
                             ); 

            // clean all the primitives   
            best_primitives.clearAllPrimitives() ;
            // clear the cloudManager 
            cloudManagerParts.clearCloud() ;
        });        

       viewer.ngui->addButton("Compute Normals",
                               [&]() {

        cloudManagerOldMesh.getCloud(typeMesh).setNormals(
                acq::recalcNormals(
                        /* [in]      K-neighbours for FLANN: */ kNeighbours,
                        /* [in]             Vertices matrix: */ cloudManagerOldMesh.getCloud(typeMesh).getVertices(),
                        /* [in]      max neighbour distance: */ maxNeighbourDist
                )
        );

        // Estimate neighbours using FLANN
        acq::NeighboursT const neighbours =
                acq::calculateCloudNeighboursFromFaces(
                        /* [in] Faces: */ cloudManagerOldMesh.getCloud(typeMesh).getFaces()
                );

        // Estimate normals for points in cloud vertices
        cloudManagerOldMesh.getCloud(typeMesh).setNormals(
                acq::calculateCloudNormals(
                        /* [in]               Cloud: */ cloudManagerOldMesh.getCloud(typeMesh).getVertices(),
                        /* [in] Lists of neighbours: */ neighbours
                )
        );

        int nFlips =
                acq::orientCloudNormalsFromFaces(
                        /* [in    ] Lists of neighbours: */ cloudManagerOldMesh.getCloud(typeMesh).getFaces(),
                        /* [in,out]   Normals to change: */ cloudManagerOldMesh.getCloud(typeMesh).getNormals()
                );

        viewer.data.clear() ;

        // Show mesh
        viewer.data.set_mesh(
            cloudManagerOldMesh.getCloud(typeMesh).getVertices(),
            cloudManagerOldMesh.getCloud(typeMesh).getFaces()
        );

        // Update viewer
        acq::setViewerNormals(
                /* [in, out] Viewer to update: */ viewer,
                /* [in]            Pointcloud: */ cloudManagerOldMesh.getCloud(typeMesh).getVertices(),
                /* [in] Normals of Pointcloud: */ cloudManagerOldMesh.getCloud(typeMesh).getNormals()
        );

        });


        viewer.ngui->addGroup("Choose the parameters");

        // ask for the number of global iteration  
        viewer.ngui->addVariable<int>("Number of iteration : ",[&] (int val) {
                nbIteration = val; }, 
                [&]() { 
                    return nbIteration; 
        } ); 

        viewer.ngui->addVariable<int>("Sample per iteration : ",[&] (int val) {
                samplePerIt = val; }, 
                [&]() { 
                    return samplePerIt; 
        } );  
        
        viewer.ngui->addVariable<double>("Threshold for distance : ",[&] (double val) {
                thresh = val; }, 
                [&]() { 
                    return thresh; 
        } );

        viewer.ngui->addVariable<double>("Threshold for angles : ",[&] (double val) {
                alpha = val; }, 
                [&]() { 
                    return alpha; 
        } );                    
        
        viewer.ngui->addVariable<double>("Threshold for best primitive (%): ",[&] (double val) {
                thresh_best = val; }, 
                [&]() { 
                    return thresh_best; 
        } );   

          viewer.ngui->addButton("RANSAC",
                               [&]() {
            // set the vertices
            acq::DecoratedCloud& thisCloud = cloudManagerOldMesh.getCloud(typeMesh) ;

             // apply RANSAC 
             bool ransacSuccess = ransac(thisCloud, best_primitives, cloudManagerParts, 
                thresh, alpha, thresh_best, nbIteration, samplePerIt) ;

                if (ransacSuccess) {
            // fuse the result in the new cloud 
            acq::DecoratedCloud& newCloud = gatherClouds(cloudManagerParts) ;

             viewer.data.clear() ;

            // Show mesh
            viewer.data.set_mesh(
                             newCloud.getVertices(),
                             newCloud.getFaces()
            );

             viewer.data.set_colors(newCloud.getColors());
                }
                else {
                    std::cout << "RANSAC didn't find any primitive" << std::endl ;
                }
        });

          viewer.ngui->addButton("Primitive fusion",
                               [&]() {
            // ******** find values for the threshold *******$
            double T_rad = 0.02 ;
            double T_cent = 0.02 ;
            double T_norm = 0.02 ;
            double T_refPt = 0.02 ;

            // fuse the similar primitive in cloud manager 
            fuse(best_primitives, cloudManagerParts, T_rad, T_cent, T_norm, T_refPt) ;

           // fuse the result in the new cloud 
            acq::DecoratedCloud& newCloud = gatherClouds(cloudManagerParts) ;

            // visualisation 
            viewer.data.clear() ;

            // Show mesh
            viewer.data.set_mesh(
                             newCloud.getVertices(),
                             newCloud.getFaces()
            );

             viewer.data.set_colors(newCloud.getColors());

        });         

        // Generate menu
        viewer.screen->performLayout();

        return false;
    }; //...viewer menu


    // Start viewer
    viewer.launch();

    return 0;
} //...main()
