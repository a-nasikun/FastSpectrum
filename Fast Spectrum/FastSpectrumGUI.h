#pragma once
#ifndef FAST_SPECTRUM_GUI
#define FAST_SPECTRUM_GUI

#include <imgui/imgui.h>
#include <imgui/imgui.cpp>
#include <imgui/imgui_demo.cpp>
#include <imgui/imgui_draw.cpp>
#include <imgui_impl_glfw_gl3.h>

#include <igl/igl_inline.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/opengl/glfw/ViewerPlugin.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>
#include <igl/isolines.h>

#include "FastSpectrum.h"
#include "App_MeshFilter.h"
#include "App_VariableOperator.h"
#include "App_DiffusionDistance.h"

int					basisToShow		= 0;
int					eigToShow		= 1;
bool				boolShowSamples = false;
static bool			boolDiffDist	= false, boolMeshFilter = false, boolVarOperator = false;
static int			numOfSample		= 1000;
static int			sampleType		= Sample_Poisson_Disk; 
static int			dataToShow		= 0;
static int			filterType		= Filter_LowPass;
static int			diffDistNumEigs = 100, numAllEigs;
static float		varOpT			= 0.25f;
static vector<double> diffDistT(1, 0.1);

static vector<vector<Eigen::VectorXd>> DiffusionTensor;
static Eigen::VectorXd diffVector;

static VariableOperator varOperator;

void showMenu(igl::opengl::glfw::Viewer &viewer, igl::opengl::glfw::imgui::ImGuiMenu &menu, FastSpectrum &fastSpectrum) {
	// MAIN WINDOW 	
	{
		// Variable for Processing
		Eigen::MatrixXd redEigVects, V = viewer.data().V;
		Eigen::VectorXd redEigVals;
		Eigen::MatrixXi F = viewer.data().F;
		
		// Variables For Viewer visualization
		float	menuWindowLeft		= 0.f,
				menuWindowWidth		= 200.f,
				menuWindowHeight	= 600.f;
		
		Eigen::MatrixXd		vColor;
		Eigen::VectorXd		Z;

		ImGui::SetNextWindowPos(ImVec2(5+menuWindowLeft * menu.menu_scaling(), 5), ImGuiSetCond_FirstUseEver);
		ImGui::SetNextWindowSize(ImVec2(menuWindowWidth, menuWindowHeight), ImGuiSetCond_FirstUseEver);
		ImGui::Begin("Fast Spectrum##2", nullptr, ImGuiWindowFlags_NoSavedSettings);

		float w = ImGui::GetContentRegionAvailWidth();
		float p = ImGui::GetStyle().FramePadding.x;

		// Diffusion Distance
		//diffDistT.reserve(1);
		//diffDistT.push_back(0.1);

		// Main Algorithm
		if (ImGui::CollapsingHeader("Executing Algorithm", ImGuiTreeNodeFlags_DefaultOpen)){
			if (ImGui::Button("Load##Mesh", ImVec2((w) / 2.f, 0))){
				// Load mesh + display the correct mesh
				int numOfModel = viewer.data_list.size();
				viewer.open_dialog_load_mesh();
				if (viewer.data_list.size() > numOfModel) {
					viewer.data() = viewer.data_list.at(++dataToShow);
					viewer.data_list.at(dataToShow - 1).clear();
				}

				V.resize(viewer.data().V.rows(), viewer.data().V.cols());
				F.resize(viewer.data().F.rows(), viewer.data().F.cols());
				V = viewer.data().V;
				F = viewer.data().F;	
				fastSpectrum.setV(V);
				fastSpectrum.setF(F);

				viewer.core.align_camera_center(V, F);

				viewer.data().show_lines = false;
				printf("A new mesh with %d vertices (%d faces).\n", V.rows(), F.rows());
				viewer.data().set_mesh(V, F);
				//for (int i = 0; i < viewer.data_list.size(); i++) {
				//	printf("Data %d: %d (%d) \n", i, viewer.data_list.at(i).V.rows(), dataToShow);
				//}
			}
			ImGui::SameLine(0, p);
			if (ImGui::Button("Save##Mesh", ImVec2((w - p) / 2.f, 0))){
				viewer.open_dialog_save_mesh();
			}
			if (ImGui::Button("Run Algorithm", ImVec2((w), 30))){
				Eigen::SparseMatrix<double>	U;
				fastSpectrum.computeEigenPairs(V, F, numOfSample, (SamplingType) sampleType, U, redEigVects, redEigVals);
				
				// Show the first non-zero eigenvectors
				eigToShow					= 1;
				Z							= U * redEigVects.col(eigToShow);
				igl::jet(Z, true, vColor);
				viewer.data().set_colors(vColor);
			}
			
			ImGui::Spacing(); ImGui::Separator(); ImGui::Spacing(); ImGui::Spacing();

			/* STEP BY STEP CONSTRUCTION */
			if (ImGui::CollapsingHeader("Program Step by Step")) {
				/* Setup Laplace-Beltrami Matrix */
				if (ImGui::Button("[1] Construct Laplacian Matrix", ImVec2(w, 0))) {
					fastSpectrum.constructLaplacianMatrix();
				}
				/* For Sampling */
				if (ImGui::Button("[2] Sampling", ImVec2(w, 0))) {
					Eigen::SparseMatrix<double> M;
					fastSpectrum.getMassMatrix(M);
					fastSpectrum.getV(V);
					if (M.rows() != V.rows() || V.rows() < 1) {
						cout << "Error! Make sure you construct Laplacian matrices before creating samples." << endl;
						return;
					}

					if (numOfSample < 1) {
						cout << "Error! Sample cannot be smaller than 1. The program will assign 1000 as the default sampling." << endl;
						numOfSample = 1000;
					}
					else if (numOfSample > V.rows()) {
						cout << "Error! Sample cannot be larger than the number of vertices in the mesh. The program will assign 1000 as the default sampling.\n" << endl;
						numOfSample = 1000;
					}

					fastSpectrum.setSample(numOfSample, (SamplingType)sampleType);
					fastSpectrum.constructSample();
				}
				ImGui::RadioButton("Poisson Disk", &sampleType, Sample_Poisson_Disk); ImGui::SameLine();
				ImGui::RadioButton("Farthest Point", &sampleType, Sample_Farthest_Point);
				//ImGui::RadioButton("Random", &sampleType, Sample_Random);

				ImGui::Text("Number of sample"); ImGui::SameLine();
				static char numSample[6] = "1000"; ImGui::InputText("", numSample, 6);
				numOfSample = atoi(numSample);

				ImGui::Spacing(); ImGui::Spacing();

				/* For Basis */
				if (ImGui::Button("[3] Construct Basis", ImVec2(w, 0))) {
					Eigen::SparseMatrix<double> M;
					Eigen::VectorXi samples;
					int numSample;
					double avgEdgeLen;
					fastSpectrum.getMassMatrix(M);
					fastSpectrum.getSamples(samples);
					fastSpectrum.getSampleSize(numSample);
					fastSpectrum.getAvgEdgeLength(avgEdgeLen);
					if (numSample < 1 || V.rows() != M.rows() || samples.maxCoeff() >= V.rows() || !(avgEdgeLen > 0.0)) {
						cout << "Error! Make sure you execute previous steps before constructing basis functions." << endl;
						return;
					}

					fastSpectrum.constructBasis();
				}
				/* For Low-Dim Problem */
				if (ImGui::Button("[4] Construct Low-Dim Eigenproblem", ImVec2(w, 0))) {
					Eigen::SparseMatrix<double> Basis;
					int numSample;
					fastSpectrum.getFunctionBasis(Basis);
					fastSpectrum.getSampleSize(numSample);
					if (Basis.cols() != numSample || Basis.rows()!=V.rows() || !(Basis.cols()>0)) {
						cout << "Error! Make sure you execute previous steps before constructing Low-Dim Eigenproblem." << endl;
						return;
					}

					fastSpectrum.constructRestrictedProblem();
				}
				/* Solving Low-Dim Eigenproblem */
				if (ImGui::Button("[5] Solve Low-Dim Eigenproblem", ImVec2(w, 0))) {
					Eigen::SparseMatrix<double> MBar, Basis;
					fastSpectrum.getReducedMassMatrix(MBar);
					fastSpectrum.getFunctionBasis(Basis);
					if (Basis.cols() != MBar.cols() || MBar.cols()<1) {
						cout << "Error! Make sure you execute previous steps before solving the Low-Dim Eigenproblem." << endl;
						return;
					}

					fastSpectrum.solveRestrictedProblem();

					// Show the first non-zero eigenvectors
					Eigen::SparseMatrix<double>	U;
					fastSpectrum.getFunctionBasis(U);
					fastSpectrum.getReducedEigVects(redEigVects);
					printf("Approximated (lifted) eigenvectors %dx%d\n", U.rows(), redEigVects.cols());
					eigToShow					= 1;
					Z							= U * redEigVects.col(eigToShow);
					igl::jet(Z, true, vColor);
					viewer.data().set_colors(vColor);
				}
			}
		}
				
		ImGui::Spacing(); ImGui::Separator(); ImGui::Spacing(); ImGui::Spacing();

		/* VISUALIZATION WINDOW */
		if (ImGui::CollapsingHeader("Visualization##FastSpectrum", ImGuiTreeNodeFlags_DefaultOpen)){
			/* For Sampling */
			if (ImGui::Button("Show Samples", ImVec2((w), 0))) {
				Eigen::VectorXi Sample;
				fastSpectrum.getSamples(Sample);
				boolShowSamples = !boolShowSamples;

				if (Sample.maxCoeff() > V.rows()) {
					cout << "Error! Make sure you compute new samples before viewing them." << endl;
					return;
				}

				if (boolShowSamples) {	
					for (int i = 0; i < Sample.size(); i++) {
						viewer.data().add_points(V.row(Sample(i)), Eigen::RowVector3d(0.7, 0.1, 0.1));
					}
				}
				else {
					viewer.data().clear();
					fastSpectrum.getV(V);
					fastSpectrum.getF(F);
					viewer.data().set_mesh(V, F);
				}

				viewer.data().point_size = 5.0f;
			}

			/* For Basis */			
			ImGui::Text("Show Basis ");
			ImGui::SameLine(0, p);
			if (ImGui::Button("- 1##Basis", ImVec2((w - 85.0f - 5.0*p) /4.f, 0))) {
				Eigen::SparseMatrix<double> U;
				fastSpectrum.getFunctionBasis(U);

				if (U.rows() != V.rows()) {
					cout << "Error! Make sure you have computed basis function on the current model." << endl;
					return;
				}

				if (basisToShow > 0)
					basisToShow--;

				if (boolShowSamples) {
					boolShowSamples = false;
					viewer.data().clear();
					viewer.data().set_mesh(V, F);
				}

				Z = U.col(basisToShow);
				igl::jet(Z, true, vColor);
				viewer.data().set_colors(vColor);
			}
			ImGui::SameLine(0, p);	
			ImGui::Text("[%d]", basisToShow);
			ImGui::SameLine(0, p);

			if (ImGui::Button("+ 1##Basis", ImVec2((w - 85.0f - 5.0*p) / 4.f, 0))) {
				Eigen::SparseMatrix<double> U;
				fastSpectrum.getFunctionBasis(U);

				if (U.rows() != V.rows()) {
					cout << "Error! Make sure you have computed basis function on the current model." << endl;
					return;
				}

				if (basisToShow < U.cols())
					basisToShow++;

				if (boolShowSamples) {
					boolShowSamples = false;
					viewer.data().clear();
					viewer.data().set_mesh(V, F);
				}

				Z = U.col(basisToShow);
				igl::jet(Z, true, vColor);
				viewer.data().set_colors(vColor);
			}
			ImGui::SameLine(0, p);
			if (ImGui::Button("Rand", ImVec2((w - 85.0f - 5.0*p)*2/4.f, 0))) {
				Eigen::SparseMatrix<double> U;
				fastSpectrum.getFunctionBasis(U);

				if (U.rows() != V.rows()) {
					cout << "Error! Make sure you have computed basis function on the current model." << endl;
					return;
				}
				
				basisToShow = rand() % U.cols();

				if (boolShowSamples) {
					boolShowSamples = false;
					viewer.data().clear();
					viewer.data().set_mesh(V, F);
				}

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

				if (U.rows() != V.rows()) {
					cout << "Error! Make sure you have computed the eigenfunctions on the current model." << endl;
					return;
				}

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
				
				if (U.rows() != V.rows()) {
					cout << "Error! Make sure you have computed the eigenfunctions on the current model." << endl;
					return;
				}

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
		
		ImGui::Spacing(); ImGui::Separator(); ImGui::Spacing(); ImGui::Spacing();

		/* APPLICATION LAUNCHER */
		if (ImGui::CollapsingHeader("Application Launcher")) {
			float w = ImGui::GetContentRegionAvailWidth();
			float p = ImGui::GetStyle().FramePadding.x;

			/* For Sampling */
			if (ImGui::Button("Launch \"Diffusion Distance\"", ImVec2((w), 0))) {
				boolDiffDist = !boolDiffDist;
				if (boolDiffDist) {

					fastSpectrum.getReducedEigVals(redEigVals);
					numAllEigs = redEigVals.size();

					boolMeshFilter = false;
					boolVarOperator = false;
				}
				std::cout << "Launching Diffusion Distance!" << std::endl;
			}

			if (ImGui::Button("Launch \"Mesh Filtering\"", ImVec2((w), 0))) {
				boolMeshFilter = !boolMeshFilter;
				if (boolMeshFilter) {
					boolDiffDist = false;
					boolVarOperator = false;
				}
				std::cout << "Launching Mesh Filtering" << std::endl;
			}

			if (ImGui::Button("Launch \"Variable-Operator\"", ImVec2((w), 0))) {
				boolVarOperator = !boolVarOperator;
				if (boolVarOperator) {
					boolDiffDist = false;
					boolMeshFilter = false;
				}
				std::cout << "Launching Variable-Operator" << std::endl;
			}
		}


		// WINDOW FOR DIFFUSION DISTANCE		
		if (boolDiffDist) {
			int windWidth, windHeight;
			glfwGetWindowSize(viewer.window, &windWidth, &windHeight);
			ImGui::SetNextWindowPos(ImVec2(((float)windWidth - 5.0f - menuWindowWidth)*menu.menu_scaling(), 5.0f), ImGuiSetCond_FirstUseEver);
			//ImGui::SetNextWindowPos(ImVec2((5.0f + menuWindowLeft + menuWindowWidth + 5.0f) * menu.menu_scaling(), 5.0f), ImGuiSetCond_FirstUseEver);
			ImGui::SetNextWindowSize(ImVec2(menuWindowWidth, 200), ImGuiSetCond_FirstUseEver);
			ImGui::Begin("Diffusion Distance", &boolDiffDist, ImGuiWindowFlags_NoSavedSettings);

			float w2 = ImGui::GetContentRegionAvailWidth();
			float p2 = ImGui::GetStyle().FramePadding.x;

			fastSpectrum.getReducedEigVals(redEigVals);
			int numAllEigs = (int)redEigVals.size();
			ImGui::SliderInt("Eigenpairs to use ##forDiffusionDistance", &diffDistNumEigs, 1, numAllEigs);
			static char tValue[6] = "0.1"; ImGui::InputText("Value of t", tValue, 6);
			diffDistT.at(0) = atof(tValue);
			//diffDistT.push_back(atof(tValue));
			//printf("size of t=%d\n", diffDistT.size());

			Eigen::MatrixXd appEigVecs; 
			Eigen::VectorXd appEigVals; 


			if (ImGui::Button("Compute Diffusion Distance", ImVec2(w2, 30))) {
				fastSpectrum.getApproxEigVects(appEigVecs);
				cout << "Get the eigenvectors: " << appEigVecs.rows() << " by " << appEigVecs.cols() << endl;
				fastSpectrum.getReducedEigVals(appEigVals);
				cout << "Get the eigenvalues: " << appEigVals.size() << endl; 
				constructDiffusionTensor(appEigVecs, appEigVals, V, diffDistNumEigs, diffDistT, DiffusionTensor);
			}

			ImGui::Spacing(); ImGui::Spacing();

			if (ImGui::Button("Color Visualization", ImVec2(w2,0))) {
				cout << "Get the eigenvalues: " << appEigVals.size() << endl;
				//if (appEigVals.size()<1) {
				//	cout << "Error! Please compute the eigenpairs first." << endl;
				//	return;
				//}
				//else 
				if (DiffusionTensor[0].size() < 1) {
					cout << "Error! Please compute the diffusion tensor first." << endl;
					return;
				}
				const int sourceID = 100;
				//visualizeDiffusionDist(V, DiffusionTensor, 0, rand() % V.rows(), diffVector);
				visualizeDiffusionDist(V, DiffusionTensor, 0, sourceID, diffVector);

				Z = diffVector;
				igl::jet(Z, true, vColor);
				viewer.data().set_colors(vColor);
			}

			//if (ImGui::Button("Isolines Visualization", ImVec2(w2, 0))) {
			//	if (DiffusionTensor[0].size() < 1) {
			//		cout << "Error! Please compute the diffusion tensor first." << endl;
			//		return;
			//	}
			//	int vertexID = (int) (rand() % V.rows());
			//	visualizeDiffusionDist(V, DiffusionTensor, 0, vertexID, diffVector);
			//	
			//	printf("Size of diffvector=%d, V=%dx%d, F=%dx%d\n.", diffVector.size(), V.rows(), V.cols(), F.rows(), F.cols());
			//	Eigen::MatrixXd isoV, isoE;
			//	Eigen::VectorXd aaa = diffVector;
			//
			//	viewer.data().clear();
			//	viewer.data().set_mesh(V, F);
			//	viewer.data().add_points(V.row(vertexID), Eigen::RowVector3d(0.7, 0.1, 0.1));
			//	igl::isolines(V, F, diffVector, 10, isoV, isoE);
			//	//igl::isolines(V, F, aaa, 10, isoV, isoE);
			//	printf("Size of isolines=%d\n.", isoE.size());
			//
			//	//for (int i = 0; i < isoE.rows(); i++) {
			//	//	viewer.data().add_edges(isoV.row(isoE.row(i)(0)), isoV.row(isoE.row(i)(1)), Eigen::RowVector3d(0.6, 0.1, 0.1));
			//	//}
			//}
			ImGui::End();
		}

		// WINDOW FOR MESH FILTERING
		if (boolMeshFilter) {	
			int windWidth, windHeight;
			glfwGetWindowSize(viewer.window, &windWidth, &windHeight);
			//ImGui::SetNextWindowPos(ImVec2((5.0f + menuWindowLeft + menuWindowWidth + 5.0f) * menu.menu_scaling(), 5.0f), ImGuiSetCond_FirstUseEver);
			ImGui::SetNextWindowPos(ImVec2(((float)windWidth-5.0f-menuWindowWidth)*menu.menu_scaling(), 5.0f), ImGuiSetCond_FirstUseEver);
			ImGui::SetNextWindowSize(ImVec2(menuWindowWidth, 225), ImGuiSetCond_FirstUseEver);
			ImGui::Begin("Mesh Filtering", &boolMeshFilter, ImGuiWindowFlags_NoSavedSettings);

			float w2 = ImGui::GetContentRegionAvailWidth();
			float p2 = ImGui::GetStyle().FramePadding.x;

			ImGui::RadioButton("Low-Pass Linear", &filterType, Filter_Linear);
			ImGui::RadioButton("Low-Pass Quadratic", &filterType, Filter_LowPass);
			ImGui::RadioButton("Amplify Higher Freq", &filterType, Filter_Polynomial);
			ImGui::RadioButton("Amplify Middle Freq", &filterType, Filter_QuadMiddle);
			ImGui::RadioButton("Suppress Middle Freq", &filterType, Filter_QuadMiddleInverse);

			Eigen::SparseMatrix<double> MassM, BasisM;
			Eigen::MatrixXd Vnew, Vold;

			static int lowPassLimit, highPassLimit, numAllEigs;
			fastSpectrum.getReducedEigVals(redEigVals);
			numAllEigs = (int) redEigVals.size();
			//highPassLimit = numAllEigs;
			//numAllEigs = 1000;

			
			ImGui::SliderInt("Low pass ##forDiffusionDistance", &lowPassLimit, 1, numAllEigs);
			ImGui::SliderInt("High pass ##forDiffusionDistance", &highPassLimit, lowPassLimit, numAllEigs);

			if (ImGui::Button("Compute Mesh Filter", ImVec2(3.0f*(w2-p2)/4.0f,30))) {
				fastSpectrum.getV(Vold);
				fastSpectrum.getMassMatrix(MassM);
				fastSpectrum.getFunctionBasis(BasisM);
				fastSpectrum.getReducedEigVects(redEigVects);
				if (!(Vold.rows()>0) || Vold.rows() != MassM.rows() || Vold.rows() != BasisM.rows() || redEigVects.cols() != BasisM.cols()) {
					cout << "Error! Please compute the eigenvectors of current model before doing Mesh Filtering." << endl;
					return;
				}
				//constructMeshFilter(Vold, MassM, BasisM, redEigVects, (FilterType)filterType, (int) (0.4*redEigVects.cols()), (int)(0.8*redEigVects.cols()), redEigVects.cols(), Vnew);
				constructMeshFilter(Vold, MassM, BasisM, redEigVects, (FilterType)filterType, lowPassLimit, highPassLimit, redEigVects.cols(), Vnew);
				viewer.data().set_mesh(Vnew, F);
				viewer.data().compute_normals();
			}  
			ImGui::SameLine();
			if (ImGui::Button("Reset", ImVec2(1.0f*(w2 - p2) / 4.0f, 30))) {
				fastSpectrum.getV(Vold);
				printf("Vold=%dx%d; F=%dx%d\n", Vold.rows(), Vold.cols(), F.rows(), F.cols());
				viewer.data().set_mesh(Vold, F);
				viewer.data().compute_normals();
			}

			ImGui::End();
		}

		// WINDOW FOR VARIABLE OPERATOR
		if (boolVarOperator) {
			int windWidth, windHeight;
			glfwGetWindowSize(viewer.window, &windWidth, &windHeight);
			ImGui::SetNextWindowPos(ImVec2(((float)windWidth - 5.0f - menuWindowWidth)*menu.menu_scaling(), 5.0f), ImGuiSetCond_FirstUseEver);
			//ImGui::SetNextWindowPos(ImVec2((5.0f+menuWindowLeft+menuWindowWidth+5.0f) * menu.menu_scaling(), 5.0f), ImGuiSetCond_FirstUseEver);
			ImGui::SetNextWindowSize(ImVec2(menuWindowWidth, 150), ImGuiSetCond_FirstUseEver);
			ImGui::Begin("Variable Operator", &boolVarOperator, ImGuiWindowFlags_NoSavedSettings);

			float w2 = ImGui::GetContentRegionAvailWidth();
			float p2 = ImGui::GetStyle().FramePadding.x;
			
			// For Variable-Operator
			Eigen::MatrixXd appEigVects;
			
			if (ImGui::Button("Initiate Variable Operator", ImVec2((w2), 30))) {
				varOperator.constructVariableOperator(V, F, 75, Sample_Poisson_Disk, varOpT, appEigVects);
				varOperator.recomputeVarOpEigVects(varOpT, appEigVects);
				Z = appEigVects.col(eigToShow);
				if (eigToShow>0)
					igl::jet(Z, true, vColor);
				else
					igl::jet(Z, false, vColor);		// for zero-th eigenvector, it doesn't make sense to scale, since they're supposed to be constant
				viewer.data().set_colors(vColor);
			}

			ImGui::Spacing(); ImGui::Spacing();

			ImGui::SliderFloat("\'t\' value", &varOpT, 0.0f, 1.0f, "%.2f");
			if (ImGui::Button("Recompute", ImVec2((w2), 0))) {
				int sampleSize;
				varOperator.getSampleSize(sampleSize);
				if (sampleSize < 1) {
					cout << "Error! You need to 'Initiate variable Operator' first before re-compute the eigenvectors." << endl; 
					return;
				}
				varOperator.recomputeVarOpEigVects(varOpT, appEigVects);

				Z = appEigVects.col(eigToShow);
				if (eigToShow>0)
					igl::jet(Z, true, vColor);
				else
					igl::jet(Z, false, vColor);		// for zero-th eigenvector, it doesn't make sense to scale, since they're supposed to be constant
				viewer.data().set_colors(vColor);
			}

			/* For Eigenvectors */
			ImGui::Text("Show Eigenvect");
			ImGui::SameLine(0, p2);
			if (ImGui::Button("- 1##EigVec2", ImVec2((w2 - 90.0f - 3.0f * p2) / 2.f, 0))) {
				varOperator.getApproxEigVects(appEigVects);
				printf("[%d] Approximated lifted eigenvectors: %dx%d\n", eigToShow, appEigVects.rows(), appEigVects.cols());
				if (eigToShow > 0)
					eigToShow--;

				Z = appEigVects.col(eigToShow);
				if (eigToShow>0)
					igl::jet(Z, true, vColor);
				else
					igl::jet(Z, false, vColor);		// for zero-th eigenvector, it doesn't make sense to scale, since they're supposed to be constant
				viewer.data().set_colors(vColor);
			}
			ImGui::SameLine(0, p2);
			ImGui::Text("[%d]", eigToShow);
			ImGui::SameLine(0, p2);

			if (ImGui::Button("+ 1##EigVec2", ImVec2((w2 - 90.0f - 3.0f * p2) / 2.f, 0))) {	
				varOperator.getApproxEigVects(appEigVects);
				printf("[%d] Approximated lifted eigenvectors: %dx%d\n", eigToShow, appEigVects.rows(), appEigVects.cols());
				if (eigToShow < appEigVects.cols())
					eigToShow++;
				Z = appEigVects.col(eigToShow);
				igl::jet(Z, true, vColor);
				viewer.data().set_colors(vColor);
			}
			ImGui::End();
		}


		ImGui::End();
	}
}

#endif // !FAST_SPECTRUM_GUI

