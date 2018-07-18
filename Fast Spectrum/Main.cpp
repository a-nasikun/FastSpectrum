#include "Utility.h"
#include "Sampling.h"
#include "Basis.h"
#include "EigenSolver.h"
#include "Visualization.h"

#define GLEW_STATIC

/* [Definitions of Global Variables] */
Eigen::SparseMatrix<double>		S, M, S_, M_, U, DistanceTableSpM;
Eigen::MatrixXd					V, LDEigVec, EEigVec, REigVec;
Eigen::MatrixXi					F;
Eigen::VectorXd					LDEigVal, REigVal;
Eigen::VectorXi					Sample;
vector<set<int>>				AdM;
igl::viewer::Viewer				viewer;
double							maxNeighDist, avgEdgeLength;
int								sampleSize;
MyMesh							opMeshData;

/* For Visulization/Trivial stuffs */
int				eigenToShow = 0;

/* [FUNCTIONS DECLARATION] */
void readMesh(const string &meshFile, Eigen::MatrixXd &V, Eigen::MatrixXi &F, MyMesh &opMeshData);
void constructLaplacianMatrix(Eigen::MatrixXd &V, Eigen::MatrixXi &F, Eigen::SparseMatrix<double> &S, Eigen::SparseMatrix<double> &M,
	vector<set<int>> &AdM, double &avgEdgeLength, Eigen::SparseMatrix<double> &DistanceTableSpM);
void constructSample(const Eigen::MatrixXd &V, igl::viewer::Viewer &viewer, vector<set<int>> &AdM, int &sampleSize, Eigen::VectorXi &Sample);
void constructBasis(const Eigen::MatrixXd &V, MyMesh &opMeshData, const Eigen::VectorXi &Sample, const vector<set<int>> AdM, const Eigen::SparseMatrix<double> DistanceTableSpM,
	const int sampleSize, const double maxNeighDist, Eigen::SparseMatrix<double> &U);
void formPartitionOfUnity(Eigen::SparseMatrix<double> &U);
void computeEigenPair(Eigen::SparseMatrix<double> &S_, Eigen::SparseMatrix<double> &M_, Eigen::MatrixXd &LDEigVec, Eigen::VectorXd &LDEigVal);

/* [MAIN FUNCTION THAT CALLS EVERYTHING ELSE] */
int main(int argc, char *argv[])
{
	cout << "===============================================" << endl;
	cout << "======== Welcome to Fast Spectrum =============" << endl;
	cout << "========== Visualization tools ================" << endl;
	cout << "Press \"B\" or \"b\" \t: Show 1 (random index) basis." << endl;
	cout << "Press \"S\" or \"s\" \t: Show sample points." << endl;
	cout << "Press \"1\" \t\t: Show a lower eigenvector (i=i-1)" << endl;
	cout << "Press \"2\" \t\t: Show a larger eigenvector (i=i+1)" << endl;
	cout << "===============================================" << endl;

	// [1]		INITIALIZATION
	// [1.1]	Reading the Mesh
	string meshFile = "../Models/AIM894_Chinese Dragon/894_Chinese Dragon.obj";
	readMesh(meshFile, V, F, opMeshData);
	
	// [1.2]	Constructing Laplacian Matrices (Stiffness and Mass-Matrix)
	constructLaplacianMatrix(V, F, S, M, AdM, avgEdgeLength, DistanceTableSpM);

	// [2]		SAMPLING
	constructSample(V, viewer, AdM, sampleSize, Sample);

	// [3]		BASIS CONSTRUCTION
	double	distRatio	= sqrt(pow(0.7, 2) + pow(0.7, 2));
	maxNeighDist		= distRatio * sqrt((double)V.rows() / (double)Sample.size()) * avgEdgeLength;
	constructBasis(V, opMeshData, Sample, AdM, DistanceTableSpM, sampleSize, maxNeighDist, U);
	formPartitionOfUnity(U);

	// [4]		LOW-DIM PROBLEM
	S_ = U.transpose() * S  * U;
	M_ = U.transpose() * M  * U;

	// [5]		SOLVING LOW-DIM EIGENPROBLEM
	computeEigenPair(S_, M_, LDEigVec, LDEigVal); 

	/* User interaction via keyboard */
	const auto &key_down = [](igl::viewer::Viewer &viewer, unsigned char key, int mod)->bool
	{
		int basisID = rand() % Sample.size();

		switch (key)
		{
		case 'S':
		case 's':
			visualizeSample(viewer, V, F, Sample);
			break;
		case 'B':
		case 'b':
			visualizeBasis(viewer, V, F, U, basisID);
			break;
		case '1':
			eigenToShow = max(0, eigenToShow -1);
			visualizeEigenvector(viewer, V, F, U, LDEigVec, eigenToShow);
			break;
		case '2':
			eigenToShow = min((int) (LDEigVec.cols()-1), eigenToShow+1);
			visualizeEigenvector(viewer, V, F, U, LDEigVec, eigenToShow);
			break;

		default:
			return false;
		}
		//viewer.data.set_vertices(V);
		viewer.data.compute_normals();
		viewer.core.align_camera_center(V, F);
		return true;
	};
	viewer.callback_init = [&](igl::viewer::Viewer &viewer){return false;};
	viewer.data.set_mesh(V, F);
	viewer.callback_key_down		= key_down;

	Eigen::Vector4f bgCol(1.0, 1.0, 1.0, 1.0);
	viewer.core.background_color	= bgCol;
	viewer.core.point_size			= 10.0f;
	viewer.core.line_width			= 0.5f;
	
	return viewer.launch();
}

/* [Read mesh model from file] */
void readMesh(const string &meshFile, Eigen::MatrixXd &V, Eigen::MatrixXi &F, MyMesh &opMeshData)
{
	V.resize(0, 0);
	F.resize(0, 0);

	if (meshFile.substr(meshFile.find_last_of(".") + 1) == "off") {
		igl::readOFF(meshFile, V, F);
		OpenMesh::IO::read_mesh(opMeshData, meshFile);
	}
	else if (meshFile.substr(meshFile.find_last_of(".") + 1) == "obj") {
		igl::readOBJ(meshFile, V, F);
		OpenMesh::IO::read_mesh(opMeshData, meshFile);
	}
	else {
		cout << "Error! File type can be either .OFF or .OBJ only." << endl;
		cout << "Program will exit in 2 seconds." << endl;
		Sleep(2000);
		exit(10);
	}
}

/* [Construct the stiffness & mass matrix AND other important tables] */
void constructLaplacianMatrix(Eigen::MatrixXd &V, Eigen::MatrixXi &F, Eigen::SparseMatrix<double> &S, Eigen::SparseMatrix<double> &M,
	vector<set<int>> &AdM, double &avgEdgeLength, Eigen::SparseMatrix<double> &DistanceTableSpM)
{
	/* Construct Stiffness matrix S */
	igl::cotmatrix(V, F, S);
	S = -S;

	/* Construct Mass matrix M*/
	igl::massmatrix(V, F, igl::MASSMATRIX_TYPE_BARYCENTRIC, M);

	/* Constructing Average length of the edges */
	avgEdgeLength = igl::avg_edge_length(V, F);

	/* Constructing Adjacency Matrix (neighborhood) */
	initiateAdM(V, F, AdM);

	/* Constructing table of edge-length using sparse Matrix (faster than vector-of-map and vector-of-unordered map) */
	initiateDistTableSpM(V, AdM, DistanceTableSpM);
}

/* [Construct samples for the subspace] */
void constructSample(const Eigen::MatrixXd &V, igl::viewer::Viewer &viewer, vector<set<int>> &AdM, int &sampleSize, Eigen::VectorXi &Sample)
{
	sampleSize = 1000;
	int nBox = 13; 

	/* Different type of sampling */
	constructPoissonDiskSample(viewer, V, nBox, avgEdgeLength, Sample);
	//constructVoxelSample(viewer, V, nBox, Sample);	
	//constructRandomSample(Sample, V, sampleSize);
	sampleSize = Sample.size();
}

/* [Construct Basis Matrix] */
void constructBasis(const Eigen::MatrixXd &V, MyMesh &opMeshData, const Eigen::VectorXi &Sample, const vector<set<int>> AdM, const Eigen::SparseMatrix<double> DistanceTableSpM,
	const int sampleSize, const double maxNeighDist, Eigen::SparseMatrix<double> &U)
{
	U.resize(V.rows(), Sample.size());
	U.reserve(20 * V.rows());

	vector<vector<Eigen::Triplet<double>>>	UTriplet;
	vector<Eigen::Triplet<double>>			AllTriplet;
	int				eCounter = 0;
	double			aCoef = 2.0/ (maxNeighDist*maxNeighDist*maxNeighDist);
	double			bCoef = -3.0 / (maxNeighDist*maxNeighDist);
	int				i = 0;
	Eigen::VectorXi	Si(1);

	vector<Eigen::Triplet<double>>	emptyTriplet;
	const int						NUM_OF_THREADS = omp_get_num_procs();
	omp_set_num_threads(NUM_OF_THREADS);
	UTriplet.resize(NUM_OF_THREADS);

	int tid, ntids, ipts, istart, iproc;
	chrono::high_resolution_clock::time_point	t4, t5;
	chrono::duration<double>					time_span4;

	vector<int> elementCounter;
	elementCounter.resize(NUM_OF_THREADS);
	int totalElementInsertions = 0;

	/* Manually break the chunck */
#pragma omp parallel private(tid,ntids,ipts,istart,i, t4,t5,time_span4)
	{
		t4			= chrono::high_resolution_clock::now();
		iproc		= omp_get_num_procs();
		tid			= omp_get_thread_num();
		ntids		= omp_get_num_threads();
		ipts		= (int)ceil(1.00*(double)sampleSize / (double)ntids);
		istart		= tid * ipts;
		if (tid == ntids - 1) ipts = sampleSize - istart;
		if (ipts <= 0) ipts = 0;

		// Instantiate D once, used for all Dijkstra
		Eigen::VectorXd				D(V.rows());
		for (int i = 0; i < V.rows(); i++) {
			D(i) = numeric_limits<double>::infinity();
		}

		//printf("Thread %d handles %d samples, starting from %d-th sample.\n", tid, ipts, istart);

		UTriplet[tid].reserve(2.0 * ((double)ipts / (double)sampleSize) * 20.0 * V.rows());

		for (i = istart; i < (istart + ipts); i++) {
			if (i >= sampleSize) break;
			/* Dijkstra basis function */
			ComputeDijkstraCompact(V, Sample(i), AdM, D, aCoef, bCoef, maxNeighDist, i, sampleSize, UTriplet[tid], elementCounter[tid]);			

			/* STVD Basis Function */
			//ComputeSTVDBasis(viewer, opMeshData, V, Sample(i), AdM, DistanceTableSpM, D, aCoef, bCoef, maxNeighDist, i, sampleSize, UTriplet[tid]);
		}
	} /* End of OpenMP parallelization */

	int totalElements = 0;
	for (int j = 0; j < NUM_OF_THREADS; j++) {
		totalElements += UTriplet[j].size();
		totalElementInsertions += elementCounter[j];
	}

	/* Combining the Triplet */
	AllTriplet.resize(totalElements);
#pragma omp parallel for
	for (int j = 0; j < NUM_OF_THREADS; j++) {
		int tripSize = 0;
		for (int k = 0; k < j; k++) {
			tripSize += UTriplet[k].size();
		}
		std::copy(UTriplet[j].begin(), UTriplet[j].end(), AllTriplet.begin() + tripSize);
	}

	U.setFromTriplets(AllTriplet.begin(), AllTriplet.end());
}

/* [Form partition of unity for the basis matrix U] */
void formPartitionOfUnity(Eigen::SparseMatrix<double> &U)
{
	/* Temporary variables */
	vector<Eigen::Triplet<double>>		diagTrip;
	Eigen::SparseMatrix<double>			DiagNormMatrix;	
	Eigen::VectorXd						rowSum = Eigen::VectorXd::Zero(U.rows());
	DiagNormMatrix.resize(U.rows(), U.rows());

	/* Getting the sum of every non-zero elements in a row */
	for (int k = 0; k < U.outerSize(); ++k) {
		for (Eigen::SparseMatrix<double>::InnerIterator it(U, k); it; ++it) {
			rowSum(it.row()) += it.value();
		}
	}

	/* Get the diagonal matrix, reciprocal of the total sum per row */
	for (int i = 0; i < rowSum.size(); i++) {
		diagTrip.push_back(Eigen::Triplet<double>(i, i, 1.0 / rowSum(i)));
	}
	DiagNormMatrix.setFromTriplets(diagTrip.begin(), diagTrip.end());

	/* Form partition of unity */
	U = DiagNormMatrix * U;
}

/* [Compute the eigenpairs of Low-Dim Problem ] */
void computeEigenPair(Eigen::SparseMatrix<double> &S_, Eigen::SparseMatrix<double> &M_, Eigen::MatrixXd &LDEigVec, Eigen::VectorXd &LDEigVal)
{
	computeEigenGPU(S_, M_, LDEigVec, LDEigVal);
	//computeEigenMatlab(S_, M_, LDEigVec, LDEigVal);
}