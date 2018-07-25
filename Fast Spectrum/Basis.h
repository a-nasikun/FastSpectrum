#pragma once
#ifndef BASIS_H
#define BASIS_H

#include <set>
#include <queue>
#include <vector>
#include <iostream>

#include <omp.h>			// OpenMP, parallelizing basis construction

#include <Eigen/Eigen>

#include "Utility.h"

using namespace std;

/* [DIJKSTRA-based BASIS] */
/* Computing Dijkstra Distance (single input, distance vector D initialized inside, stops when reaching a certain distance. 
** COMPACT Version, meaning that distance vector D is only initialized when necessary **/
void ComputeDijkstraCompact(const Eigen::MatrixXd &V, const int s1, const vector<set<int>> &AdM,
	Eigen::VectorXd &D, const double a, const double b, const double nDist, const int col, const int sampleSize,
	vector<Eigen::Triplet<double>> &UTriplet, int &elementCounter);	

#endif // !BASIS_H

