#include "acq/normalEstimation.h"
#include "acq/decoratedCloud.h"
#include "acq/cloudManager.h"
#include "acq/gestion.h"

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

    /* ====================================== */
    /* ===== TODO 1) : LIRE FICHIER OBJ ===== */
    /* ====================================== */


    // Load a mesh in OFF format
    std::string meshPath = "../models/sphere_cube.obj";
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
    acq::CloudManager cloudManager;
    // Read mesh from meshPath
    {
        Eigen::MatrixXd V, TC, N;
        Eigen::MatrixXi F, FTC, FN;
        igl::readOBJ(meshPath, V, TC, N, F, FTC, FN);

        if(FTC.size() == 0) FTC = F;

        if (V.rows() <= 0) {
            std::cerr << "Could not read mesh at " << meshPath
                      << "...exiting...\n";
            return EXIT_FAILURE;
        }


        // Store read vertices and faces
        N.rowwise().normalize();
        cloudManager.addCloud(acq::DecoratedCloud(V, F, N));

        // Show mesh
        viewer.data.set_mesh(
            cloudManager.getCloud(0).getVertices(),
            cloudManager.getCloud(0).getFaces()
        );

        // Set Normals from OBJ file
        cloudManager.getCloud(0).setNormals(N);
        std::cout << N.size() << std::endl;

        // Update viewer
        acq::setViewerNormals(viewer, cloudManager.getCloud(0).getVertices(), N);
    }

    // Extend viewer menu using a lambda function
    viewer.callback_init =
        [
            &cloudManager, &kNeighbours, &maxNeighbourDist,
            &floatVariable, &boolVariable, &dir
        ] (igl::viewer::Viewer& viewer)
    {
        // Add an additional menu window
        viewer.ngui->addWindow(Eigen::Vector2i(900,10), "Acquisition3D");



        // Generate menu
        viewer.screen->performLayout();

        return false;
    }; //...viewer menu


    // Start viewer
    viewer.launch();

    return 0;
} //...main()
