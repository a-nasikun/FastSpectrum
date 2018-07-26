#pragma once
#ifndef FAST_SPECTRUM_H
#define FAST_SPECTRUM_H

#include "Basis.h"
#include "Utility.h"
#include "Sampling.h"
#include "EigenSolver.h"
#include "Visualization.h"

/* [LIBIGL HEADERS] */
#include <igl/jet.h>
#include <igl/readOFF.h>
#include <igl/readOBJ.h>
#include <igl/isolines.h>
#include <igl/cotmatrix.h>
#include <igl/doublearea.h>
#include <igl/massmatrix.h>
#include <igl/barycenter.h>
#include <igl/avg_edge_length.h>
#include <igl/opengl/glfw/Viewer.h>

using namespace std;
//#define GLEW_STATIC

class FastSpectrum
{
public:
	/* [MAIN FUNCTIONS CALLED BY USER FOR ABSTRACTION] */

	/* [MAIN FUNCTIONS IN FAST APPROXIMATIONG ALGORITHM] */
	void computeEigenPairs(const string &meshFile, const int &numOfSamples, Eigen::MatrixXd &reducedEigVects, Eigen::VectorXd &reducedEigVals);
	void computeEigenPairs(const string &meshFile, const int &numOfSamples, Eigen::SparseMatrix<double> &Basis, Eigen::MatrixXd &reducedEigVects, Eigen::VectorXd &reducedEigVals);
	void computeEigenPairs(const string &meshFile, const int &numOfSamples, SamplingType sampleType, Eigen::MatrixXd &reducedEigVects, Eigen::VectorXd &reducedEigVals);
	void computeEigenPairs(const string &meshFile, const int &numOfSamples, SamplingType sampleType, Eigen::SparseMatrix<double> &Basis, Eigen::MatrixXd &reducedEigVects, Eigen::VectorXd &reducedEigVals);
	void computeEigenPairs(Eigen::MatrixXd &V, Eigen::MatrixXi &F, const int &numOfSamples, Eigen::MatrixXd &reducedEigVects, Eigen::VectorXd &reducedEigVals);
	void computeEigenPairs(Eigen::MatrixXd &V, Eigen::MatrixXi &F, const int &numOfSamples, Eigen::SparseMatrix<double> &Basis, Eigen::MatrixXd &reducedEigVects, Eigen::VectorXd &reducedEigVals);
	void computeEigenPairs(Eigen::MatrixXd &V, Eigen::MatrixXi &F, const int &numOfSamples, SamplingType sampleType, Eigen::MatrixXd &reducedEigVects, Eigen::VectorXd &reducedEigVals);
	void computeEigenPairs(Eigen::MatrixXd &V, Eigen::MatrixXi &F, const int &numOfSamples, SamplingType sampleType, Eigen::SparseMatrix<double> &Basis, Eigen::MatrixXd &reducedEigVects, Eigen::VectorXd &reducedEigVals);
	
	void constructLaplacianMatrix();
	void constructSample();
	void constructBasis();
	void constructRestrictedProblem();
	void solveRestrictedProblem();

	void getV(Eigen::MatrixXd &V);
	void getF(Eigen::MatrixXi &F);
	void getSamples(Eigen::VectorXi &Sample);
	void getFunctionBasis(Eigen::SparseMatrix<double> &U);
	void getReducedEigVals(Eigen::VectorXd &reducedEigVals);
	void getReducedEigVects(Eigen::MatrixXd &reducedEigVects);
	void getApproxEigVects(Eigen::MatrixXd &approxEigVects);
	void getReducedLaplacian(Eigen::SparseMatrix<double> Sbar, Eigen::SparseMatrix<double> Mbar);
	void setV(const Eigen::MatrixXd &V);
	void setF(const Eigen::MatrixXi &F);
	void setMesh(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F);
	void setMesh(const string &meshFile);
	void setSample(const int &sampleSize, SamplingType sampleType);
	void setSampleNumber(const int &sampleSize);
	void setSampleType(SamplingType sampleType);

	/* [FUNCTIONS DECLARATION] */
	void readMesh(const string &meshFile, Eigen::MatrixXd &V, Eigen::MatrixXi &F);
	void constructBasisFunctions();

private:
	/* [CLAS VARIABLES] */
	Eigen::SparseMatrix<double>		S, M, S_, M_, Basis, DistanceTableSpM;
	Eigen::MatrixXd					V, reducedEigVects;
	Eigen::MatrixXi					F;
	Eigen::VectorXd					reducedEigVals;
	Eigen::VectorXi					Sample;
	vector<set<int>>				AdM;
	double							maxNeighDist, avgEdgeLength;
	int								sampleSize;
	SamplingType					sampleType;

	/* [SUPPLEMENTARY FUNCTIONS] */
	void formPartitionOfUnity(Eigen::SparseMatrix<double> &Basis);
	void computeEigenPair();
	void normalizeReducedEigVects();
};





#endif // !FAST_SPECTRUM_H

