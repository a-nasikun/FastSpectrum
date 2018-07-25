#pragma once
#ifndef SAMPLING_H
#define SAMPLING_H

#include <set>
#include <queue>
#include <time.h>
#include <random>
#include <chrono>
#include <iostream>

#include <igl/opengl/glfw/Viewer.h>

#include "Utility.h"

using namespace std;

//#include "Utility.h"

/* Data structure for box numbering in Poisson-Disk */
struct BoxStruct{
	int		id;
	double	distance;
	bool	operator> (const BoxStruct &ref) const { return distance > ref.distance; }
	bool	operator< (const BoxStruct &ref) const { return distance < ref.distance; }
};

void constructRandomSample(Eigen::VectorXi &Sample, const Eigen::MatrixXd &V, int n);
void constructVoxelSample(igl::opengl::glfw::Viewer &viewer, const Eigen::MatrixXd &V, const int &n, Eigen::VectorXi &Sample);
void constructPoissonDiskSample(const Eigen::MatrixXd &V, const int &n, const double &avgEdgeLength, Eigen::VectorXi &Sample);

#endif // !SAMPLING_H

