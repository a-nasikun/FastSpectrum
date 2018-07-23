#pragma once
#ifndef BASIS_H
#define BASIS_H

#include "Utility.h"

/* [DIJKSTRA-based BASIS] */
// Computing Dijkstra Distance (sing input, D initialized inside, stops when reaching a certain distance, 
// and used for Basis Construction), COMPACT Version
void ComputeDijkstraCompact(const Eigen::MatrixXd &V, const int s1, const vector<set<int>> &AdM,
	Eigen::VectorXd &D, const double a, const double b, const double nDist, const int col, const int sampleSize,
	vector<Eigen::Triplet<double>> &UTriplet, int &elementCounter);	

#endif // !BASIS_H

