#include "FastSpectrum.h"


#include <imgui/imgui.cpp>
#include <imgui/imgui.h>
#include <imgui/imgui_demo.cpp>
#include <imgui/imgui_draw.cpp>
//#include <imgui_impl_glfw_gl3.h>
#include <../imgui/imgui_impl_glfw_gl3.h>
#include <../imgui/imgui_impl_glfw_gl3.cpp>

#include <igl/opengl/glfw/Viewer.h>
#include <igl/opengl/glfw/ViewerPlugin.h>
#include <igl/igl_inline.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>


/* For Visulization/Trivial stuffs */
int				eigenToShow = 0;



/* [MAIN FUNCTION THAT CALLS EVERYTHING ELSE] */
int main(int argc, char *argv[])
{
	igl::opengl::glfw::Viewer viewer;
	FastSpectrum fastSpectrum;
	Eigen::MatrixXd redEigVects, V;
	Eigen::VectorXd redEigVals;
	Eigen::MatrixXi F;
	string meshFile = "../Models/AIM894_Chinese Dragon/894_Chinese Dragon.obj";


	/* VISUALIZATION */
	// Variables For Viewer visualization
	bool boolDiffDist = false, boolMeshFilter = false, boolVarOperator = false;
	float menuWindowLeft = 0.0f, menuWindowWidth = 200.0f;

	int basisToShow = 0;
	int eigToShow = 0;
	Eigen::MatrixXd		vColor;
	Eigen::VectorXd		Z;
	bool boolShowSamples = false;

	// Attach a menu plugin
	igl::opengl::glfw::imgui::ImGuiMenu menu;
	viewer.plugins.push_back(&menu);
	// Setup Dear ImGui binding
	IMGUI_CHECKVERSION();
	ImGui::CreateContext();
	// Setup style
	ImGui::StyleColorsDark();

	/* END OF VISUALIZATION */

	fastSpectrum.computeEigenPairs(meshFile, 1000, redEigVects, redEigVals);
	fastSpectrum.getV(V);
	fastSpectrum.getF(F);


	cout << "===============================================" << endl;
	cout << "======== Welcome to Fast Spectrum =============" << endl;
	cout << "========== Visualization tools ================" << endl;
	cout << "Press \"B\" or \"b\" \t: Show 1 (random index) basis." << endl;
	cout << "Press \"S\" or \"s\" \t: Show sample points." << endl;
	cout << "Press \"1\" \t\t: Show a lower eigenvector (i=i-1)" << endl;
	cout << "Press \"2\" \t\t: Show a larger eigenvector (i=i+1)" << endl;
	cout << "===============================================" << endl;
		
	menu.callback_draw_viewer_window = [&]()
	{
		// MAIN WINDOW 	
		{
			ImGui::SetNextWindowPos(ImVec2(menuWindowLeft * menu.menu_scaling(), 10), ImGuiSetCond_FirstUseEver);
			ImGui::SetNextWindowSize(ImVec2(menuWindowWidth, 220), ImGuiSetCond_FirstUseEver);
			ImGui::Begin(
				"Fast Spectrum", nullptr,
				ImGuiWindowFlags_NoSavedSettings
			);
			// Mesh
			if (ImGui::CollapsingHeader("Mesh", ImGuiTreeNodeFlags_DefaultOpen))
			{
				float w = ImGui::GetContentRegionAvailWidth();
				float p = ImGui::GetStyle().FramePadding.x;
				if (ImGui::Button("Load##Mesh", ImVec2((w) / 2.f, 0)))
				{
					// Load mesh + display the correct mesh
					viewer.open_dialog_load_mesh();
					viewer.data() = viewer.data_list.at(0);
					V = viewer.data().V;
					F = viewer.data().F;

					printf("Number of vertices: %d\n", V.rows());
				}
				ImGui::SameLine(0, p);
				if (ImGui::Button("Save##Mesh", ImVec2((w - p) / 2.f, 0)))
				{
					viewer.open_dialog_save_mesh();
				}
				if (ImGui::Button("Run Algorithm", ImVec2((w), 0)))
				{
					//fastSpectrum.computeEigenPairs(meshFile, 1000, redEigVects, redEigVals);
					fastSpectrum.computeEigenPairs(V, F, 1000, redEigVects, redEigVals);
				}

				if (ImGui::CollapsingHeader("Program Step by Step")) {
					/* Setup Laplace-Beltrami Matrix */
					if (ImGui::Button("[1] Construct Laplacian Matrix", ImVec2((w) / 2.f, 0))) {
						
					}

					/* For Sampling */
					if (ImGui::Button("[2] Sampling", ImVec2((w) / 2.f, 0))) {
						
					}

					/* For Basis */
					if (ImGui::Button("[3] Construct Basis", ImVec2((w - 4.0f*p) / 2.f, 0))) {
						
					}

					/* For Low-Dim Problem */
					if (ImGui::Button("[4] Construct Low-Dim Eigenproblem", ImVec2(w, 0))) {
						
					}
					/* Solving Low-Dim Eigenproblem */
					if (ImGui::Button("[5] Solve Low-Dim Eigenproblem", ImVec2(w, 0))) {
						
					}
				}
			}
			ImGui::End();
		}

		// VISUALIZATION WINDOW 
		{
			// Define next window position + size
			ImGui::SetNextWindowPos(ImVec2(menuWindowLeft * menu.menu_scaling(), 230), ImGuiSetCond_FirstUseEver);
			ImGui::SetNextWindowSize(ImVec2(menuWindowWidth, 100), ImGuiSetCond_FirstUseEver);
			ImGui::Begin(
				"Visualization", nullptr,
				ImGuiWindowFlags_NoSavedSettings
			);
			// Mesh

			float w = ImGui::GetContentRegionAvailWidth();
			float p = ImGui::GetStyle().FramePadding.x;

			/* For Sampling */
			if (ImGui::Button("Show Samples", ImVec2((w), 0))) {
				//boolShowSamples = !boolShowSamples;
				//if (boolShowSamples) {
				//	for (int i = 0; i < Sample.size(); i++) {
				//		viewer.data().add_points(V.row(Sample(i)), Eigen::RowVector3d(0.7, 0.1, 0.1));
				//	}
				//}
			}

			/* For Basis */

			ImGui::Text("Show Basis", ImVec2((w - 4.0f*p) / 2.f, 0));
			ImGui::SameLine(0, p);
			if (ImGui::Button("-##Basis", ImVec2((w - p) / 6.f, 0))) {
				//basisToShow = max(basisToShow - 1, 0);
				//Z = U.col(basisToShow);
				//igl::jet(Z, true, vColor);
				//viewer.data().set_colors(vColor);
			}

			ImGui::SameLine(0, p);

			if (ImGui::Button("+##Basis", ImVec2((w - p) / 6.f, 0))) {
				//basisToShow = min(basisToShow + 1, (int)Sample.size());
				//Z = U.col(basisToShow);
				//igl::jet(Z, true, vColor);
				//viewer.data().set_colors(vColor);
			}
			ImGui::SameLine(0, p);
			if (ImGui::Button("Rand", ImVec2((w - p) / 6.f, 0))) {
				//Z = U.col(rand() % Sample.size());
				//igl::jet(Z, true, vColor);
				//viewer.data().set_colors(vColor);
			}
			/* For Eigenvectors */
			ImGui::Text("Show Eigenvectors", ImVec2((w - 4.0f*p) / 2.f, 0));

			ImGui::SameLine(0, p);
			if (ImGui::Button("-##EigVec", ImVec2((w - p) / 4.f, 0))) {
				eigToShow = max(eigToShow - 1, 0);
				Eigen::SparseMatrix<double> U;
				fastSpectrum.getFunctionBasis(U);
				Z = U * redEigVects.col(eigToShow);
				igl::jet(Z, true, vColor);
				viewer.data().set_colors(vColor);
			}

			ImGui::SameLine(0, p);
			if (ImGui::Button("+##EigVec", ImVec2((w - p) / 4.f, 0))) {
				eigToShow = min(eigToShow + 1, (int)redEigVals.size());
				Eigen::SparseMatrix<double> U;
				fastSpectrum.getFunctionBasis(U);
				Z = U * redEigVects.col(eigToShow);
				igl::jet(Z, true, vColor);
				viewer.data().set_colors(vColor);
			}

			/* Overlays */
			ImGui::Checkbox("Wireframe", &(viewer.data().show_lines));
			ImGui::Checkbox("Fill Faces", &(viewer.data().show_faces));

			ImGui::End();
		}
	};


	/* User interaction via keyboard */
	const auto &key_down = [](igl::opengl::glfw::Viewer &viewer, unsigned char key, int mod)->bool
	{
		

		switch (key)
		{
		case 'S':
		case 's':
			break;
		case 'B':
		case 'b':
			break;
		case '1':
			break;
		case '2':
			break;

		default:
			return false;
		}
		//viewer.data.set_vertices(V);
		viewer.data().compute_normals();
		//viewer.core.align_camera_center(V, F);
		return true;
	};
	viewer.callback_init = [&](igl::opengl::glfw::Viewer &viewer){return false;};
	viewer.data().set_mesh(V, F);
	viewer.callback_key_down		= key_down;

	//Eigen::Vector4f bgCol(1.0, 1.0, 1.0, 1.0);
	//viewer.core.background_color	= bgCol;
	viewer.data().point_size = 5.0f;
	viewer.data().line_width = 0.5f; 	
	return viewer.launch();
}
