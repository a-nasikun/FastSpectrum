#pragma once
#ifndef APP_DIFFUSION_DISTANCE_H
#define APP_DIFFUSION_DISTANCE_H

#include <Eigen/Eigen>
#include <Eigen/SparseCore>
#include <vector>

using namespace std;

void constructDiffusionTensor(const Eigen::MatrixXd &approxEigVects, const Eigen::VectorXd &approxEigVals, const Eigen::MatrixXd &V,
	const int &numEigen, const vector<double> &t, vector<vector<Eigen::VectorXd>> &DiffusionTensor);

void constructDiffusionTensor(const Eigen::SparseMatrix<double> &Basis, const Eigen::MatrixXd &reducedEigVects, const Eigen::VectorXd &reducedEigVals, const Eigen::MatrixXd &V,
	const int &numEigen, const vector<double> &t, vector<vector<Eigen::VectorXd>> &DiffusionTensor);

void visualizeDiffusionDist(const Eigen::MatrixXd &V, const vector<vector<Eigen::VectorXd>> &DiffusionTensor, const int &tID, const int &vID, Eigen::VectorXd &diffVector);
#endif // !APP_DIFFUSION_DISTANCE_H
