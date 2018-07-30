#include "FastSpectrum.h"
#include "FastSpectrumGUI.h"

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
	fastSpectrum.computeEigenPairs(meshFile, 1000, Basis, redEigVects, redEigVals);
	fastSpectrum.getV(V);
	fastSpectrum.getF(F);
		
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
	viewer.callback_key_down		= key_down;

	Eigen::Vector4f					bgCol(0.5f, 0.4f, 0.9f, 0.2f);
	//Eigen::Vector4f				bgCol(0.95f, 0.95f, 0.95f, 1.0f);
	viewer.core.background_color	= bgCol;
	viewer.data().point_size		= 5.0f;
	viewer.data().line_width		= 0.5f; 	
	return viewer.launch();
}
