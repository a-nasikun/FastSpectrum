#pragma once
#ifndef BASIS_H
#define BASIS_H

#include "Utility.h"
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
#include <OpenMesh/Tools/Utils/getopt.h>

typedef OpenMesh::TriMesh_ArrayKernelT<>  MyMesh;

/* [DIJKSTRA-based BASIS] */
// Computing Dijkstra Distance (sing input, D initialized inside, stops when reaching a certain distance, 
// and used for Basis Construction), COMPACT Version
void ComputeDijkstraCompact(const Eigen::MatrixXd &V, const int s1, const vector<set<int>> &AdM,
	Eigen::VectorXd &D, const double a, const double b, const double nDist, const int col, const int sampleSize,
	vector<Eigen::Triplet<double>> &UTriplet, int &elementCounter);


/* [STVD-based BASIS] */
// Unfolding the edges, K times backward, unless reaching the Si (initial point) & distance in form of VECTORS
void UnfoldEdgeSTVD(igl::viewer::Viewer &viewer, MyMesh &mesh, const Eigen::MatrixXd &V, const int Si, const vector<set<int>> &AdM,
	const Eigen::SparseMatrix<double> &DistanceTableSpM, const Eigen::VectorXi &prev, const int &initNode, const int &max_i, vector<Eigen::RowVector2d>	&EVect);

// Return the minimum distance after Unfolding the edges of STVD
double MinDistance4STVD(igl::viewer::Viewer &viewer, MyMesh &mesh, const Eigen::MatrixXd &V, const int Si, const vector<set<int>> &AdM,
	const Eigen::SparseMatrix<double> &DistanceTableSpM, Eigen::VectorXd &D, Eigen::VectorXi &prev, const int &u, const int &w, const int &max_pred);

// STVD for Basis Construction
void ComputeSTVDBasis(igl::viewer::Viewer &viewer, MyMesh &mesh, const Eigen::MatrixXd &V, const int s1, const vector<set<int>> &AdM, const Eigen::SparseMatrix<double> &DistanceTableSpM,
	Eigen::VectorXd &D, const double a, const double b, const double nDist, const int col, const int sampleSize,
	vector<Eigen::Triplet<double>> &UTriplet);

#endif // !BASIS_H

