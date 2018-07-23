#include "FastSpectrum.h"

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

	Eigen::Vector4f bgCol(1.0, 1.0, 1.0, 1.0);
	viewer.core.background_color	= bgCol;
	viewer.data().point_size = 5.0f;
	viewer.data().line_width = 0.5f; 	
	return viewer.launch();
}
