#pragma once
#ifndef FAST_SPECTRUM_GUI
#define FAST_SPECTRUM_GUI

#include <imgui/imgui.h>
#include <imgui/imgui.cpp>
#include <imgui/imgui_demo.cpp>
#include <imgui/imgui_draw.cpp>
#include <imgui_impl_glfw_gl3.h>
#include <imgui_impl_glfw_gl3.cpp>

#include <igl/igl_inline.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/opengl/glfw/ViewerPlugin.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>

#include "FastSpectrum.h"

int		basisToShow=0;
int		eigToShow=0;
bool	boolShowSamples = false;

void showMenu(igl::opengl::glfw::Viewer &viewer, igl::opengl::glfw::imgui::ImGuiMenu &menu, FastSpectrum &fastSpectrum) {
	// MAIN WINDOW 	
	{
		// Variable for Processing
		Eigen::MatrixXd redEigVects, V = viewer.data().V;
		Eigen::VectorXd redEigVals;
		Eigen::MatrixXi F = viewer.data().F;
		//string meshFile = "../Models/AIM894_Chinese Dragon/894_Chinese Dragon.obj";

		// Variables For Viewer visualization
		float	menuWindowLeft = 0.f,
				menuWindowWidth = 200.f,
				menuWindowHeight = 450.f;

		
		Eigen::MatrixXd		vColor;
		Eigen::VectorXd		Z;

		ImGui::SetNextWindowPos(ImVec2(1+menuWindowLeft * menu.menu_scaling(), 1), ImGuiSetCond_FirstUseEver);
		ImGui::SetNextWindowSize(ImVec2(menuWindowWidth, menuWindowHeight), ImGuiSetCond_FirstUseEver);
		//ImGui::SetNextWindowSize(ImGui::GetIO().DisplaySize, ImGuiSetCond_Always);
		ImGui::Begin(
			"Fast Spectrum##2", nullptr,
			ImGuiWindowFlags_NoSavedSettings
		);

		float w = ImGui::GetContentRegionAvailWidth();
		float p = ImGui::GetStyle().FramePadding.x;

		// Main Algorithm
		if (ImGui::CollapsingHeader("Executing Algorithm", ImGuiTreeNodeFlags_DefaultOpen)){
			if (ImGui::Button("Load##Mesh", ImVec2((w) / 2.f, 0))){
				// Load mesh + display the correct mesh
				viewer.open_dialog_load_mesh();
				viewer.data() = viewer.data_list.at(0);
				V.resize(viewer.data().V.rows(), viewer.data().V.cols());
				F.resize(viewer.data().F.rows(), viewer.data().F.cols());
				V = viewer.data().V;
				F = viewer.data().F;	
				fastSpectrum.setV(V);
				fastSpectrum.setF(F);

				viewer.data().show_lines = false;
				printf("A new mesh with %d vertices (%d faces).\n", V.rows(), F.rows());
			}
			ImGui::SameLine(0, p);
			if (ImGui::Button("Save##Mesh", ImVec2((w - p) / 2.f, 0))){
				viewer.open_dialog_save_mesh();
			}
			if (ImGui::Button("Run Algorithm", ImVec2((w), 30))){
				//fastSpectrum.computeEigenPairs(meshFile, 1000, redEigVects, redEigVals);
				fastSpectrum.computeEigenPairs( V, F, 1000, redEigVects, redEigVals);
				
				// Show the first non-zero eigenvectors
				eigToShow					= 1;
				Eigen::SparseMatrix<double>	U;
				fastSpectrum.getFunctionBasis(U);
				Z							= U * redEigVects.col(eigToShow);
				igl::jet(Z, true, vColor);
				viewer.data().set_colors(vColor);
			}

			ImGui::Spacing(); ImGui::Spacing(); ImGui::Spacing();

			/* STEP BY STEP CONSTRUCTION */
			if (ImGui::CollapsingHeader("Program Step by Step")) {
				/* Setup Laplace-Beltrami Matrix */
				if (ImGui::Button("[1] Construct Laplacian Matrix", ImVec2(w, 0))) {
					fastSpectrum.constructLaplacianMatrix();
				}
				/* For Sampling */
				if (ImGui::Button("[2] Sampling", ImVec2(w, 0))) {
					fastSpectrum.runSampling();
				}
				/* For Basis */
				if (ImGui::Button("[3] Construct Basis", ImVec2(w, 0))) {
					fastSpectrum.constructBasis();
				}
				/* For Low-Dim Problem */
				if (ImGui::Button("[4] Construct Low-Dim Eigenproblem", ImVec2(w, 0))) {
					fastSpectrum.constructRestrictedProblem();
				}
				/* Solving Low-Dim Eigenproblem */
				if (ImGui::Button("[5] Solve Low-Dim Eigenproblem", ImVec2(w, 0))) {
					fastSpectrum.solveRestrictedProblem();

					// Show the first non-zero eigenvectors
					eigToShow					= 1;
					Eigen::SparseMatrix<double>	U;
					fastSpectrum.getFunctionBasis(U);
					fastSpectrum.getReducedEigVects(redEigVects);
					printf("Basis %dx%d\n", U.rows(), U.cols());
					printf("EigVect %dx%d\n", redEigVects.rows(), redEigVects.cols());
					Z							= U * redEigVects.col(eigToShow);
					igl::jet(Z, true, vColor);
					viewer.data().set_colors(vColor);
				}
			}
		}

		ImGui::Spacing(); ImGui::Spacing(); ImGui::Spacing();

		if (ImGui::CollapsingHeader("Visualization##FastSpectrum", ImGuiTreeNodeFlags_CollapsingHeader)){
			/* For Sampling */
			if (ImGui::Button("Show Samples", ImVec2((w), 0))) {				
				viewer.data().clear();
				fastSpectrum.getV(V);
				fastSpectrum.getF(F);
				viewer.data().set_mesh(V, F);

				Eigen::VectorXi Sample;
				fastSpectrum.getSamples(Sample);
				boolShowSamples = !boolShowSamples;

				if (boolShowSamples) {	
					for (int i = 0; i < Sample.size(); i++) {
						viewer.data().add_points(V.row(Sample(i)), Eigen::RowVector3d(0.7, 0.1, 0.1));
					}
				}
				else {
				}
			}

			/* For Basis */			
			ImGui::Text("Show Basis ");
			ImGui::SameLine(0, p);
			if (ImGui::Button("- 1##Basis", ImVec2((w - 85.0f - 5.0*p) /4.f, 0))) {
				Eigen::SparseMatrix<double> U;
				fastSpectrum.getFunctionBasis(U);
				if (basisToShow > 0)
					basisToShow--;
				Z = U.col(basisToShow);
				igl::jet(Z, true, vColor);
				viewer.data().set_colors(vColor);
			}
			ImGui::SameLine(0, p);	
			ImGui::Text("<%d>", basisToShow);
			ImGui::SameLine(0, p);

			if (ImGui::Button("+ 1##Basis", ImVec2((w - 85.0f - 5.0*p) / 4.f, 0))) {
				Eigen::SparseMatrix<double> U;
				fastSpectrum.getFunctionBasis(U);
				if (basisToShow < U.cols())
					basisToShow++;

				Z = U.col(basisToShow);
				igl::jet(Z, true, vColor);
				viewer.data().set_colors(vColor);
			}
			ImGui::SameLine(0, p);
			if (ImGui::Button("Rand", ImVec2((w - 85.0f - 5.0*p)*2/4.f, 0))) {
				Eigen::SparseMatrix<double> U;
				fastSpectrum.getFunctionBasis(U);
				basisToShow = rand() % U.cols();
				Z = U.col(basisToShow);
				igl::jet(Z, true, vColor);
				viewer.data().set_colors(vColor);
			}
			/* For Eigenvectors */
			ImGui::Text("Show EigVec");
			ImGui::SameLine(0, p);
			if (ImGui::Button("- 1##EigVec1", ImVec2((w-90.0f-3*p)/2.f, 0))) {
				fastSpectrum.getReducedEigVects(redEigVects);
				Eigen::SparseMatrix<double> U;
				fastSpectrum.getFunctionBasis(U);
				if (eigToShow > 0)
					eigToShow--;

				Z = U * redEigVects.col(eigToShow);
				if(eigToShow>0)
					igl::jet(Z, true, vColor);
				else 
					igl::jet(Z, false, vColor);		// for zero-th eigenvector, it doesn't make sense to scale, since they're supposed to be constant
				viewer.data().set_colors(vColor);
			}
			ImGui::SameLine(0, p);
			ImGui::Text("[%d]", eigToShow);
			ImGui::SameLine(0, p);

			if (ImGui::Button("+ 1##EigVec1", ImVec2((w-90.0f-3*p)/2.f, 0))) {
				fastSpectrum.getReducedEigVects(redEigVects);
				Eigen::SparseMatrix<double> U;
				fastSpectrum.getFunctionBasis(U);
				if (eigToShow < redEigVects.cols())
					eigToShow++;

				Z = U * redEigVects.col(eigToShow);
				igl::jet(Z, true, vColor);
				viewer.data().set_colors(vColor);
			}

			/* Overlays */
			ImGui::Checkbox("Wireframe", &(viewer.data().show_lines));
			ImGui::Checkbox("Fill Faces", &(viewer.data().show_faces));
		}

		ImGui::End();
	}
}

#endif // !FAST_SPECTRUM_GUI

