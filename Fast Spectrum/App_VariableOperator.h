#pragma once
#ifndef APP_VARIABLE_OPERATOR_H
#define APP_VARIABLE_OPERATOR_H

#include "FastSpectrum.h"

class VariableOperator : public FastSpectrum
{
public:
	VariableOperator();
	~VariableOperator();

	void constructVariableOperator(Eigen::MatrixXd &V, Eigen::MatrixXi &F, const int &numOfSamples, SamplingType sampleType, const int &t, Eigen::MatrixXd &approxEigVects, Eigen::VectorXd &approxEigVals);
	void constructVariableOperator(Eigen::MatrixXd &V, Eigen::MatrixXi &F, const int &numOfSamples, SamplingType sampleType, const int &t, Eigen::MatrixXd &approxEigVects);
	void recomputeVarOpEigVects(const double &t, Eigen::MatrixXd &approxEigVects);

	void constructRestrictedProblem();
	void solveRestrictedProblem();

	void getMatrixA(Eigen::SparseMatrix<double> &A);
	void getVarOpMatrix(Eigen::SparseMatrix<double> &AS);
	void getReducedVarOpMatrix(Eigen::SparseMatrix<double> &ASbar);

protected:
	Eigen::SparseMatrix<double>		A, AS, AS_;
	double							t;

	void constructMatrixA();
	void computeEigenPair();

};

#endif // !APP_VARIABLE_OPERATOR_H