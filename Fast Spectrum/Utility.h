#pragma once
#ifndef UTILITY_H
#define UTILITY_H

/* [LIBIGL HEADERS] */
#include <igl/barycenter.h>
#include <igl/cotmatrix.h>
#include <igl/doublearea.h>
#include <igl/grad.h>
#include <igl/jet.h>
#include <igl/massmatrix.h>
#include <igl/per_vertex_normals.h>
#include <igl/readDMAT.h>
#include <igl/readOFF.h>
#include <igl/readOBJ.h>
#include <igl/repdiag.h>
#include <igl/invert_diag.h>
#include <igl/avg_edge_length.h>
#include <igl/eigs.h>
#include <igl/isolines.h>
#include <igl/viewer/Viewer.h>

/* [STANDARD LIBRARIES] */
#include <iostream>
#include <fstream>
#include <set>
#include <vector>
#include <map>
#include <unordered_map>
#include <queue>
#include <time.h>
#include <random>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <algorithm>

/* [ADDITIONAL LIBRARIES] */
#include <omp.h>
#include <chrono>

/* [MATLAB] */
#include <MatlabDataArray.hpp>
#include <MatlabEngine.hpp>
#include <tuple>
#include "engine.h"
#include "mex.h"

/* [CUDA] */
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include "cusparse.h"
#include "cusolverDn.h"
#include "cusolverSp.h"

using namespace std;

/* Data structure for Priority Queue */
struct VertexPair
{
	int		vId;
	double	distance;
	bool	operator> (const VertexPair &ref) const { return distance > ref.distance; }
	bool	operator< (const VertexPair &ref) const { return distance < ref.distance; }
};

/* Create a list of neighbors of every vertex */
void initiateAdM(Eigen::MatrixXd &V, Eigen::MatrixXi &F, vector<set<int>> &AM);
/* Precompute the distance between a vertex and its neighbor(s), for a bit of speed up */
void initiateDistTableSpM(Eigen::MatrixXd &V, const vector<set<int>> &AM, Eigen::SparseMatrix<double> &DistanceTable);


/* Compute distance between two points */
void VtoVDist(const Eigen::VectorXd &V1, const Eigen::VectorXd &V2, double &dist);
/* Rotate the edge*/
void rotateVector2D(const Eigen::RowVector2d &e, const double &angle, Eigen::RowVector2d &e2);
/* Get the angle between two edges */
double VtoVangle(const Eigen::Vector3d &V1, const Eigen::Vector3d &V2, const Eigen::Vector3d &P);

#endif // !UTILITY_H

