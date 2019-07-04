#pragma once
#ifndef APP_MESH_FILTER_H
#define APP_MESH_FILTER_H

#include <Eigen/Eigen>
#include <Eigen/SparseCore>

/* Enum for different sampling method */
enum FilterType {
	Filter_LowPass				= 0,
	Filter_Linear				= 1,
	Filter_Polynomial			= 2,
	Filter_QuadMiddle			= 3,
	Filter_QuadMiddleInverse	= 4,
};

// Construct a Mesh Filter from Reduced-Eigenvectors
// INPUT: Vertex matrix & Mass Matrix M & Basis Matrix U & Reduced Eigenvectors & paramters of filter
// OUTPUT: New Vertex matrix
void constructMeshFilter(const Eigen::MatrixXd &V, const Eigen::SparseMatrix<double> &M, const Eigen::SparseMatrix<double> &U,	
	const Eigen::MatrixXd &LDEigVec, FilterType filterType, const int &k, const int &m, const int &n, Eigen::MatrixXd &Vout);

/* Obtain projection of every columns of V (x,y,z) */
void projectSpace(const Eigen::MatrixXd &Vin, const Eigen::SparseMatrix<double> &U, const Eigen::SparseMatrix<double> &M,
	const Eigen::MatrixXd &LDEigVec, const Eigen::VectorXd &filter, Eigen::MatrixXd &Vout);

/* Construct Filters */
// No Filter
void createNoFilter(const int &k, const int &m, Eigen::VectorXd &filterCoef);

// Low Pass Filter
void createLowPassFilter(const int &k, const int &m, Eigen::VectorXd &filterCoef);

// Linear Filter (1 until a certain k, then increase linearly to 2 for m)
void createLinearFilter(const int &k, const int &m, Eigen::VectorXd &filterCoef);

// Polynomial Filter (1 until a certain k, then polynom function from k to m)
void createPolynomialFilter(const int &k, const int &m, Eigen::VectorXd &filterCoef);

// Polynomial Filter (1 until a certain k, then polynom function from k to m)
void createQuadMiddleFilter(const int &k, const int &m, const int &n, Eigen::VectorXd &filterCoef);

// Polynomial Filter (1 until a certain k, then polynom function from k to m)
void createQuadMiddleFilterInverse(const int &k, const int &m, const int &n, Eigen::VectorXd &filterCoef);

#endif // !APP_MESH_FILTER_H

