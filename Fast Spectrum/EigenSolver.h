#pragma once
#ifndef EIGEN_SOLVER_H
#define EIGEN_SOLVER_H

//#include "Utility.h"

#include <Eigen/Eigen>

// For CUDA
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include "cusolverDn.h"
#include "cublas_v2.h"

#include <assert.h> 
#include <iostream>

using namespace std;

/* Computing Eigenstructure in GPU */
void computeEigenGPU(Eigen::SparseMatrix<double> &S_, Eigen::SparseMatrix<double> &M_, Eigen::MatrixXd &LDEigVec, Eigen::VectorXd &LDEigVal);

#endif // !EIGEN_SOLVER_H

