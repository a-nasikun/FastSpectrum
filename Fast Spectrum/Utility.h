#pragma once
#ifndef UTILITY_H
#define UTILITY_H

#include <set>
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>

#include <chrono>				// For quite precise timing, needed (perhaps) in most sub-task

#include <Eigen/Eigen>

using namespace std;

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

