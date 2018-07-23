#pragma once
#ifndef FAST_SPECTRUM_H
#define FAST_SPECTRUM_H

#include "Utility.h"
#include "Sampling.h"
#include "Basis.h"
#include "EigenSolver.h"
#include "Visualization.h"

#define GLEW_STATIC

class FastSpectrum
{
public:
	/* [MAIN FUNCTIONS IN FAST APPROXIMATIONG ALGORITHM] */
	void computeEigenPairs(const string &meshFile, const int &numOfSamples, Eigen::MatrixXd &reducedEigVects, Eigen::VectorXd &reducedEigVals);
	void getFunctionBasis(Eigen::SparseMatrix<double> &U);
	void getReducedLaplacian();
	void FastSpectrum::getV(Eigen::MatrixXd &V);
	void FastSpectrum::getF(Eigen::MatrixXi &F);

private:
	/* [Definitions of Global Variables] */
	Eigen::SparseMatrix<double>		S, M, S_, M_, Basis, DistanceTableSpM;
	Eigen::MatrixXd					V, reducedEigVects;
	Eigen::MatrixXi					F;
	Eigen::VectorXd					reducedEigVals;
	Eigen::VectorXi					Sample;
	vector<set<int>>				AdM;
	double							maxNeighDist, avgEdgeLength;
	int								sampleSize;

	/* [FUNCTIONS DECLARATION] */
	void readMesh(const string &meshFile, Eigen::MatrixXd &V, Eigen::MatrixXi &F);
	void constructLaplacianMatrix(Eigen::MatrixXd &V, Eigen::MatrixXi &F, Eigen::SparseMatrix<double> &S, Eigen::SparseMatrix<double> &M,
		vector<set<int>> &AdM, double &avgEdgeLength, Eigen::SparseMatrix<double> &DistanceTableSpM);
	void constructSample(const Eigen::MatrixXd &V, vector<set<int>> &AdM, int &sampleSize, Eigen::VectorXi &Sample);
	void constructBasis(const Eigen::MatrixXd &V, Eigen::MatrixXi &T, const Eigen::VectorXi &Sample, const vector<set<int>> AdM, const Eigen::SparseMatrix<double> DistanceTableSpM,
		const int sampleSize, const double maxNeighDist, Eigen::SparseMatrix<double> &U);
	void formPartitionOfUnity(Eigen::SparseMatrix<double> &U);
	void computeEigenPair(Eigen::SparseMatrix<double> &S_, Eigen::SparseMatrix<double> &M_, Eigen::MatrixXd &LDEigVec, Eigen::VectorXd &LDEigVal);
};





#endif // !FAST_SPECTRUM_H

