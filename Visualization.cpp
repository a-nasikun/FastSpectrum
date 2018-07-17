#include "Visualization.h"

/* Visualize the samples */
void visualizeSample(igl::viewer::Viewer &viewer, const Eigen::MatrixXd &V, const Eigen::MatrixXi &F, const Eigen::VectorXi &Sample)
{
	viewer.data.clear();
	viewer.data.set_mesh(V, F);

	for (int i = 0; i < Sample.size(); i++)
	{
		viewer.data.add_points(V.row(Sample(i)), Eigen::RowVector3d(1.0, 0.1, 0.0));
	}
}

/* Visualize basis */
void visualizeBasis(igl::viewer::Viewer &viewer, const Eigen::MatrixXd &V, const Eigen::MatrixXi &F, const Eigen::SparseMatrix<double> &U, const int &basisID)
{
	viewer.data.clear();
	viewer.data.set_mesh(V, F);

	// Map the basis value at particular index to the jet coloring scheme
	Eigen::MatrixXd		vColor;
	Eigen::VectorXd		Z;
	Z					= U.col(basisID);
	igl::jet(Z, true, vColor);
	viewer.data.set_colors(vColor);
}

/* Visualize Approximated Eigenvector */
void visualizeEigenvector(igl::viewer::Viewer &viewer, const Eigen::MatrixXd &V, const Eigen::MatrixXi &F, const Eigen::SparseMatrix<double> &U, const Eigen::MatrixXd &LDEigVec, const int &eigenID)
{
	viewer.data.clear();
	viewer.data.set_mesh(V, F);

	// Some validation, ensuring no out-of-bound indices
	if(eigenID<0 || eigenID>LDEigVec.cols()-1) { 
		printf("Error! The index should be between 0 and %d.\n", LDEigVec.cols()-1);
		return;
	}

	// Map the eigenvectors as the color of the mesh (on each of the vertices)
	Eigen::MatrixXd		vColor;
	Eigen::VectorXd		Z;

	Z = U * LDEigVec.col(eigenID);
	if (eigenID == 0) {
		igl::jet(Z, false, vColor);
	}else {
		igl::jet(Z, true, vColor);
	}
	viewer.data.set_colors(vColor);
}