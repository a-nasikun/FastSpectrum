#include "FastSpectrum.h"
#include "FastSpectrumGUI.h"
#include "App_MeshFilter.h"
#include "App_DiffusionDistance.h"
#include "App_VariableOperator.h"

/* [MAIN FUNCTION THAT CALLS EVERYTHING ELSE] */
int main(int argc, char *argv[])
{
	Eigen::SparseMatrix<double> Basis;
	Eigen::MatrixXd redEigVects, V;
	Eigen::VectorXd redEigVals;
	Eigen::MatrixXi F;
	string meshFile; 
	
	/* GUI FUNCTIONALITY */	
	igl::opengl::glfw::Viewer			viewer;			// Set up main viewer
	FastSpectrum						fastSpectrum;
	igl::opengl::glfw::imgui::ImGuiMenu menu;			// Attach a menu plugin
	viewer.plugins.push_back(&menu);									
	//IMGUI_CHECKVERSION();								// Setup Dear ImGui binding
	ImGui::CreateContext();
	ImGui::StyleColorsDark();							// Set style

	meshFile = "../Models/AIM894_Chinese Dragon/894_Chinese Dragon.obj";
	//fastSpectrum.computeEigenPairs(meshFile, 1000, Basis, redEigVects, redEigVals);
	fastSpectrum.computeEigenPairs(meshFile, 102, Sample_Farthest_Point, Basis, redEigVects, redEigVals);

	fastSpectrum.getV(V);
	fastSpectrum.getF(F);

	/* MESH FILTER */
	Eigen::MatrixXd Vnew;
	Eigen::SparseMatrix<double> M;
	fastSpectrum.getMassMatrix(M);
	constructMeshFilter(V, M, Basis, redEigVects, Filter_LowPass, 10, 10, Vnew);
	V = Vnew;
	printf("Mesh filter with Low Pass filter produced %d-by-%d vertices\n", V.rows(), V.cols());

	/* DIFFUSION DISTANCE */
	//vector<vector<Eigen::VectorXd>> DiffTensor;
	//vector<double> t;
	//t.push_back(0.10);
	//t.push_back(0.25);
	//t.push_back(0.50);
	//constructDiffusionTensor(Basis, redEigVects, redEigVals, V, 500, t, DiffTensor);
	//printf("Diffusion tensor with %d different values of t and 500 eigenpairs. Size=[%d, %d, %d].\n", t.size(), DiffTensor.size(), DiffTensor[0].size(), DiffTensor[0][0].size());
			
	/* VARIABLE OPERATOR */
	VariableOperator	varOp;
	Eigen::MatrixXd		approxEigVects;

	varOp.constructVariableOperator(V, F, 100, Sample_Farthest_Point, 0.5, approxEigVects);	


	menu.callback_draw_viewer_window = [&]()
	{
		showMenu(viewer, menu, fastSpectrum);		
	};
	/* END OF GUI FUNCTIONALITY */

	/* User interaction via keyboard */
	const auto &key_down = [](igl::opengl::glfw::Viewer &viewer, unsigned char key, int mod)->bool
	{
		switch (key)
		{
		default:
			return false;
		}
		//viewer.data.set_vertices(V);
		viewer.data().compute_normals();
		//viewer.core.align_camera_center(V, F);
		return true;
	};

	viewer.callback_init			= [&](igl::opengl::glfw::Viewer &viewer){return false;};
	viewer.data().set_mesh(V, F);

	Eigen::VectorXd Z = approxEigVects.col(2);
	Eigen::MatrixXd vColor;
	igl::jet(Z, true, vColor);
	printf("approx size=%dx%d, Z=%d, color=%dx%d, V=%dx%d\n", approxEigVects.rows(), approxEigVects.cols(), Z.size(), vColor.rows(), vColor.cols(), V.rows(), V.cols());
	viewer.data().set_colors(vColor);

	viewer.callback_key_down		= key_down;

	Eigen::Vector4f					bgCol(0.5f, 0.4f, 0.9f, 0.2f);
	//Eigen::Vector4f				bgCol(0.95f, 0.95f, 0.95f, 1.0f);
	viewer.core.background_color	= bgCol;
	viewer.data().point_size		= 5.0f;
	viewer.data().line_width		= 0.5f; 	
	return viewer.launch();
}
