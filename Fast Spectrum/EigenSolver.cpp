#include "EigenSolver.h"


/* Computing Eigenstructure in GPU */
/* This code is mostly taken from CUDA's cuSOLVER's example codes */
void computeEigenGPU(Eigen::SparseMatrix<double> &S_, Eigen::SparseMatrix<double> &M_, Eigen::MatrixXd &LDEigVec, Eigen::VectorXd &LDEigVal)
{
	cusolverDnHandle_t	cusolverH = NULL;
	cusolverStatus_t	cusolver_status = CUSOLVER_STATUS_SUCCESS;
	cudaError_t			cudaStat1 = cudaSuccess;
	cudaError_t			cudaStat2 = cudaSuccess;
	cudaError_t			cudaStat3 = cudaSuccess;
	cudaError_t			cudaStat4 = cudaSuccess;

	const int m = S_.rows();
	const int lda = m;

	double *A = (double*)std::malloc(m*lda * sizeof(double));
	double *B = (double*)std::malloc(m*lda * sizeof(double));
	Eigen::MatrixXd Snew(S_);
	Eigen::MatrixXd Mnew(M_);
	A = Snew.data();
	B = Mnew.data();

	double	*V = (double*)std::malloc(m*lda * sizeof(double)); // [lda*m];		// eigenvectors 
	double	*W = (double*)std::malloc(m * sizeof(double)); 			// eigenvalues 
	double	*d_A = NULL;
	double	*d_B = NULL;
	double	*d_W = NULL;
	int		*devInfo = NULL;
	double	*d_work = NULL;
	int		lwork = 0;
	int		info_gpu = 0;

	// step 1: create cusolver/cublas handle 
	cusolver_status = cusolverDnCreate(&cusolverH);
	assert(CUSOLVER_STATUS_SUCCESS == cusolver_status);

	// step 2: copy A and B to device 
	cudaStat1 = cudaMalloc((void**)&d_A, sizeof(double) * lda * m);
	cudaStat2 = cudaMalloc((void**)&d_B, sizeof(double) * lda * m);
	cudaStat3 = cudaMalloc((void**)&d_W, sizeof(double) * m);
	cudaStat4 = cudaMalloc((void**)&devInfo, sizeof(int));
	if (cudaStat1 != 0 || cudaStat2 != 0 || cudaStat3 != 0 || cudaStat4 !=0) {
		cout << "Error! Cannot allocate memory in GPU." << endl; 
	}	
	assert(cudaSuccess == cudaStat1);
	assert(cudaSuccess == cudaStat2);
	assert(cudaSuccess == cudaStat3);
	assert(cudaSuccess == cudaStat4);
	
	cudaStat1			= cudaMemcpy(d_A, A, sizeof(double) * lda * m, cudaMemcpyHostToDevice);
	cudaStat2			= cudaMemcpy(d_B, B, sizeof(double) * lda * m, cudaMemcpyHostToDevice);
	assert(cudaSuccess == cudaStat1);
	assert(cudaSuccess == cudaStat2);
	if (cudaStat1 != 0 || cudaStat2 != 0) {
		cout << "Error! Cannot copy data to GPU." << endl;
	}

	// step 3: query working space of sygvd 
	cusolverEigType_t	itype = CUSOLVER_EIG_TYPE_1;		// A*x = (lambda)*B*x 
	cusolverEigMode_t	jobz = CUSOLVER_EIG_MODE_VECTOR;	// compute eigenvalues and eigenvectors. 
	cublasFillMode_t	uplo = CUBLAS_FILL_MODE_UPPER;
	cusolver_status = cusolverDnDsygvd_bufferSize(cusolverH, itype, jobz, uplo, m, d_A, lda, d_B, lda, d_W, &lwork);
	assert(cusolver_status == CUSOLVER_STATUS_SUCCESS);
	cudaStat1 = cudaMalloc((void**)&d_work, sizeof(double)*lwork);
	assert(cudaSuccess == cudaStat1);
	if (cudaStat1 != 0) {
		cout << "Error! Cannot request space to solve the eigenproblem." << endl;
	}

	// step 4: compute spectrum of (A,B) 
	cusolver_status = cusolverDnDsygvd(cusolverH, itype, jobz, uplo, m, d_A, lda, d_B, lda, d_W, d_work, lwork, devInfo);
	cudaStat1 = cudaDeviceSynchronize();
	assert(CUSOLVER_STATUS_SUCCESS == cusolver_status);
	assert(cudaSuccess == cudaStat1);
	if (cudaStat1 != 0 || cudaStat2 != 0) {
		cout << "Error! CUDA cuSolver's eigensolver FAILED to solve the eigenproblem." << endl;
	}

	cudaStat1 = cudaMemcpy(W, d_W, sizeof(double)*m, cudaMemcpyDeviceToHost);
	cudaStat2 = cudaMemcpy(V, d_A, sizeof(double)*lda*m, cudaMemcpyDeviceToHost);
	cudaStat3 = cudaMemcpy(&info_gpu, devInfo, sizeof(int), cudaMemcpyDeviceToHost);
	if (cudaStat1 != 0 || cudaStat2 != 0 || cudaStat3 != 0) {
		cout << "Error! Cannot copy data from GPU to CPU." << endl;
	}
	assert(cudaSuccess == cudaStat1);
	assert(cudaSuccess == cudaStat2);
	assert(cudaSuccess == cudaStat3);

	if (info_gpu != 0) {
		cout << "Error! There's some error with sygvd." << endl;
	}
	assert(0 == info_gpu);

	// Copying back to EigenFormat (dense) for Eigenvalues and Eigenvectors
	LDEigVec = Eigen::Map<Eigen::MatrixXd, 0, Eigen::OuterStride<>>(V, m, lda, Eigen::OuterStride<>(m));
	LDEigVal = Eigen::Map<Eigen::VectorXd>(W, m);

	// free resources 
	if (d_A) cudaFree(d_A);
	if (d_B) cudaFree(d_B);
	if (d_W) cudaFree(d_W);
	if (devInfo) cudaFree(devInfo);
	if (d_work) cudaFree(d_work);
	if (cusolverH) cusolverDnDestroy(cusolverH);
	cudaDeviceReset();
}