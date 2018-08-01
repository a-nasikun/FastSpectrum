#include "App_VariableOperator.h"

VariableOperator::VariableOperator()
{

}

VariableOperator::~VariableOperator()
{

}

void VariableOperator::constructVariableOperator(Eigen::MatrixXd &V, Eigen::MatrixXi &F, const int &numOfSamples, SamplingType sampleType, const int &t, Eigen::MatrixXd &approxEigVects, Eigen::VectorXd &approxEigVals)
{
	constructVariableOperator(V, F, numOfSamples, sampleType, t, reducedEigVects);
	approxEigVals = this->reducedEigVals;
}

void VariableOperator::constructVariableOperator(Eigen::MatrixXd &V, Eigen::MatrixXi &F, const int &numOfSamples, SamplingType sampleType, const int &t, Eigen::MatrixXd &approxEigVects)
{
	// [1.2] Constructing Laplacian Matrices (Stiffness, Mass, and Variable Operator Matrix)
	this->V = V;
	this->F = F;
	constructLaplacianMatrix();

	// [1.3] Constructing the Curvature-aware Variable Operator Matrix
	constructMatrixA();

	// [2]	SAMPLING
	this->setSample(numOfSamples, sampleType);
	constructSample();

	// [3]	BASIS CONSTRUCTION
	double	distRatio = sqrt(pow(0.7, 2) + pow(0.7, 2));
	maxNeighDist = distRatio * sqrt((double)V.rows() / (double)Sample.size()) * avgEdgeLength;
	constructBasis();

	// [3.1] Recontruct the VarOp Matrix
	this->t = t;
	AS = t*A + (1.0 - t)*S;

	// [4]	LOW-DIM PROBLEM
	constructRestrictedProblem();

	// [5]	SOLVING LOW-DIM EIGENPROBLEM
	solveRestrictedProblem();

	// OUTPUT
	approxEigVects = Basis * reducedEigVects;
}

void VariableOperator::recomputeVarOpEigVects(const int &t, Eigen::MatrixXd &approxEigVects)
{
	this->t = t;
	AS = t*A + (1.0 - t)*S;

	// [4]	LOW-DIM PROBLEM
	constructRestrictedProblem();

	// [5]	SOLVING LOW-DIM EIGENPROBLEM
	solveRestrictedProblem();

	// OUTPUT
	approxEigVects = Basis * reducedEigVects;
}

void VariableOperator::constructRestrictedProblem()
{
	AS_ = Basis.transpose() * AS  * Basis;
	M_ = Basis.transpose() * M  * Basis;

	printf("A set of reduced VarOp and mass matrices (each %dx%d) is constructed.\n", AS_.rows(), M_.cols());
}

void VariableOperator::solveRestrictedProblem() {
	computeEigenPair();
	normalizeReducedEigVects();
}

void VariableOperator::getMatrixA(Eigen::SparseMatrix<double> &A)
{
	A = this->A;
}

void VariableOperator::getVarOpMatrix(Eigen::SparseMatrix<double> &AS)
{
	AS = this->AS;
}
void VariableOperator::getReducedVarOpMatrix(Eigen::SparseMatrix<double> &ASbar)
{
	ASbar = this->AS_;
}

void VariableOperator::constructMatrixA()
{
	Eigen::MatrixXd					NV;
	Eigen::MatrixXi					edges;
	vector<Eigen::Triplet<double>>	ATriplet;

	igl::per_vertex_normals(V, F, NV);
	A.resize(S.rows(), S.cols());
	A.reserve(20 * V.rows());

	// Construct Matrix A
	for (int k = 0; k < S.outerSize(); ++k) {
		for (Eigen::SparseMatrix<double>::InnerIterator it(S, k); it; ++it) {
			int		row = it.row();
			int		col = it.col();
			double	val = it.value(), ASij;

			double Aij = NV.row(row).dot(NV.row(col));
			if (row == col) {
				ASij = val;
			}
			else {
				ASij = Aij * val;
			}

			ATriplet.push_back(Eigen::Triplet<double>(row, col, ASij));
		}
	}

	A.setFromTriplets(ATriplet.begin(), ATriplet.end());
	printf("A size is %dx%d, with %d non Zeros\n", A.rows(), A.cols(), A.nonZeros());
}

void VariableOperator::computeEigenPair()
{
	computeEigenGPU(AS_, M_, reducedEigVects, reducedEigVals);
}
