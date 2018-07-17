#include "EigenSolver.h"


/* Computing Eigenstructure in GPU */
void computeEigenGPU(Eigen::SparseMatrix<double> &S_, Eigen::SparseMatrix<double> &M_, Eigen::MatrixXd &LDEigVec, Eigen::VectorXd &LDEigVal)
{
	cusolverDnHandle_t	cusolverH = NULL;
	cusolverStatus_t	cusolver_status = CUSOLVER_STATUS_SUCCESS;
	cudaError_t			cudaStat1 = cudaSuccess;
	cudaError_t			cudaStat2 = cudaSuccess;
	cudaError_t			cudaStat3 = cudaSuccess;
	cudaError_t			cudaStat4 = cudaSuccess;

	chrono::high_resolution_clock::time_point	t1, t2;
	chrono::duration<double>					time_span;

	t1 = chrono::high_resolution_clock::now();

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
	t2 = chrono::high_resolution_clock::now();
	time_span = chrono::duration_cast<chrono::duration<double>>(t2 - t1);
	printf("Time to create handle in GPU is %.5f. \n", time_span);

	// step 2: copy A and B to device 
	cudaStat1 = cudaMalloc((void**)&d_A, sizeof(double) * lda * m);
	cudaStat2 = cudaMalloc((void**)&d_B, sizeof(double) * lda * m);
	cudaStat3 = cudaMalloc((void**)&d_W, sizeof(double) * m);
	cudaStat4 = cudaMalloc((void**)&devInfo, sizeof(int));
	cout << "Status of Allocating Memory" << endl;
	cout << cudaStat1 << ", " << cudaStat2 << ", " << cudaStat3 << endl;
	
	assert(cudaSuccess == cudaStat1);
	assert(cudaSuccess == cudaStat2);
	assert(cudaSuccess == cudaStat3);
	assert(cudaSuccess == cudaStat4);
	cudaStat1 = cudaMemcpy(d_A, A, sizeof(double) * lda * m, cudaMemcpyHostToDevice);
	cudaStat2 = cudaMemcpy(d_B, B, sizeof(double) * lda * m, cudaMemcpyHostToDevice);
	assert(cudaSuccess == cudaStat1);
	assert(cudaSuccess == cudaStat2);
	cout << "Status of Copying CPU => GPU" << endl;
	cout << cudaStat1 << ", " << cudaStat2 << endl;

	// step 3: query working space of sygvd 
	cusolverEigType_t	itype = CUSOLVER_EIG_TYPE_1;		// A*x = (lambda)*B*x 
	cusolverEigMode_t	jobz = CUSOLVER_EIG_MODE_VECTOR;	// compute eigenvalues and eigenvectors. 
	cublasFillMode_t	uplo = CUBLAS_FILL_MODE_UPPER;
	cusolver_status = cusolverDnDsygvd_bufferSize(cusolverH, itype, jobz, uplo, m, d_A, lda, d_B, lda, d_W, &lwork);
	assert(cusolver_status == CUSOLVER_STATUS_SUCCESS);
	cudaStat1 = cudaMalloc((void**)&d_work, sizeof(double)*lwork);
	assert(cudaSuccess == cudaStat1);
	cout << "Status of asking working space" << endl;
	cout << cudaStat1 << endl;
	
	// step 4: compute spectrum of (A,B) 
	cusolver_status = cusolverDnDsygvd(cusolverH, itype, jobz, uplo, m, d_A, lda, d_B, lda, d_W, d_work, lwork, devInfo);
	cudaStat1 = cudaDeviceSynchronize();
	assert(CUSOLVER_STATUS_SUCCESS == cusolver_status);
	assert(cudaSuccess == cudaStat1);

	cout << "Status of eigensolver" << endl;
	cout << cusolver_status << ", " << cudaStat1 << endl;

	cudaStat1 = cudaMemcpy(W, d_W, sizeof(double)*m, cudaMemcpyDeviceToHost);
	cudaStat2 = cudaMemcpy(V, d_A, sizeof(double)*lda*m, cudaMemcpyDeviceToHost);
	cudaStat3 = cudaMemcpy(&info_gpu, devInfo, sizeof(int), cudaMemcpyDeviceToHost);
	
	cout << "Status of Copying" << endl;
	cout << cudaStat1 << ", " << cudaStat2 << ", " << cudaStat3 << endl;

	assert(cudaSuccess == cudaStat1);
	assert(cudaSuccess == cudaStat2);
	assert(cudaSuccess == cudaStat3);
	printf("after sygvd: info_gpu = %d\n", info_gpu);
	assert(0 == info_gpu);

	// Copying back to EigenFormat (dense) for Eigenvalues and Eigenvectors
	LDEigVec = Eigen::Map<Eigen::MatrixXd, 0, Eigen::OuterStride<>>(V, m, lda, Eigen::OuterStride<>(m));
	LDEigVal = Eigen::Map<Eigen::VectorXd>(W, m);
	//Eigen::MatrixXd newMatA(a.rows(), a.cols());
	//newMatA = Eigen::Map<Eigen::MatrixXd, 0, Eigen::OuterStride<>>(stdFormat, a.rows(), a.cols(), Eigen::OuterStride<>(a.rows()));

	t2 = chrono::high_resolution_clock::now();
	time_span = chrono::duration_cast<chrono::duration<double>>(t2 - t1);
	printf("Time to copy data to convert to EigenFormat is %.5f. \n", time_span);

	// free resources 
	if (d_A) cudaFree(d_A);
	if (d_B) cudaFree(d_B);
	if (d_W) cudaFree(d_W);
	if (devInfo) cudaFree(devInfo);
	if (d_work) cudaFree(d_work);
	if (cusolverH) cusolverDnDestroy(cusolverH);
	cudaDeviceReset();
}

/* Computing Eigenstructure in Matlab */
void computeEigenMatlab(Eigen::SparseMatrix<double> &S, Eigen::SparseMatrix<double> &M, Eigen::MatrixXd &EigVec, Eigen::VectorXd &EigVal)
{
	printf("Size of S = %dx%d\n", S.rows(), S.cols());
	using namespace matlab::engine;
	Engine *ep;
	mxArray *MM = NULL, *MS = NULL, *result = NULL, *eigVecResult, *nEigs;

	chrono::high_resolution_clock::time_point	t1, t2, t3, t4;
	chrono::duration<double>					time_span, ts2;

	const int NNZ_S = S.nonZeros();
	const int NNZ_M = M.nonZeros();
	double *eigVal, *eigVec;

	// Allocate memory for S and M (sparse representation)
	double	*srs = (double*)malloc(NNZ_S * sizeof(double));
	mwIndex *irs = (mwIndex*)malloc(NNZ_S * sizeof(mwIndex));
	mwIndex *jcs = (mwIndex*)malloc((S.cols() + 1) * sizeof(mwIndex));

	double	*srm = (double*)malloc(NNZ_M * sizeof(double));
	mwIndex *irm = (mwIndex*)malloc(NNZ_M * sizeof(mwIndex));
	mwIndex *jcm = (mwIndex*)malloc((M.cols() + 1) * sizeof(mwIndex));

	// Bind MM with M, and MS with S
	MS = mxCreateSparse(S.rows(), S.cols(), NNZ_S, mxREAL);
	srs = mxGetPr(MS);
	irs = mxGetIr(MS);
	jcs = mxGetJc(MS);

	MM = mxCreateSparse(M.rows(), M.cols(), NNZ_M, mxREAL);
	srm = mxGetPr(MM);
	irm = mxGetIr(MM);
	jcm = mxGetJc(MM);

	// Setting initial variable value
	int nnzSCounter = 0;
	int nnzMCounter = 0;

	t1 = chrono::high_resolution_clock::now();

	// Getting matrix S
	jcs[0] = nnzSCounter;
	for (int i = 0; i < S.outerSize(); i++) {
		for (Eigen::SparseMatrix<double>::InnerIterator it(S, i); it; ++it) {
			srs[nnzSCounter] = it.value();
			irs[nnzSCounter] = it.row();
			nnzSCounter++;
		}
		jcs[i + 1] = nnzSCounter;
	}

	// Getting matrix M
	jcm[0] = nnzMCounter;
	for (int i = 0; i < M.outerSize(); i++) {
		for (Eigen::SparseMatrix<double>::InnerIterator it(M, i); it; ++it) {
			srm[nnzMCounter] = it.value();
			irm[nnzMCounter] = it.row();
			nnzMCounter++;
		}
		jcm[i + 1] = nnzMCounter;
	}

	// Start Matlab Engine
	ep = engOpen(NULL);
	if (!(ep = engOpen(""))) {
		fprintf(stderr, "\nCan't start MATLAB engine\n");
		cout << "CANNOT START MATLAB " << endl;
	}
	else {
		cout << "MATLAB STARTS. OH YEAH!!!" << endl;
	}

	// Compute Eigenvalue in Matlab
	int NUM_EIGEN = 50;

	engPutVariable(ep, "MS", MS);
	engPutVariable(ep, "MM", MM);

	t3 = chrono::high_resolution_clock::now();
	engEvalString(ep, "[LDEigVec,LDEigVal]=eigs(MS,MM,50,'smallestreal');");
	engEvalString(ep, "LDEigVal=diag(LDEigVal);");
	engEvalString(ep, "hold on; plot(1:50,LDEigVal(1:50),'LineWidth',1.5);"); // has to do it this way for "correct" plot
	t4 = chrono::high_resolution_clock::now();

	//engEvalString(ep, "save('D:/4_SCHOOL/TU Delft/Research/Projects/EigenTrial/MATLAB Implementation/Data/Dino_LDEigVec_15_5000','LDEigVec');");
	
	result = engGetVariable(ep, "LDEigVal");
	eigVal = (double*)malloc(NUM_EIGEN * sizeof(double));
	memcpy((void *)eigVal, (void *)mxGetPr(result), NUM_EIGEN * sizeof(double));

	eigVecResult = engGetVariable(ep, "LDEigVec");
	eigVec = (double*)malloc(M.rows() * NUM_EIGEN * sizeof(double));
	memcpy((void *)eigVec, (void *)mxGetPr(eigVecResult), M.rows() * NUM_EIGEN * sizeof(double));

	EigVal = Eigen::Map<Eigen::VectorXd>(eigVal, NUM_EIGEN);
	EigVec = Eigen::Map<Eigen::MatrixXd, 0, Eigen::OuterStride<>>(eigVec, M.rows(), NUM_EIGEN, Eigen::OuterStride<>(M.rows()));

	t2 = chrono::high_resolution_clock::now();
	time_span = chrono::duration_cast<chrono::duration<double>>(t2 - t1);
	ts2 = chrono::duration_cast<chrono::duration<double>>(t4 - t3);
	printf("Time to compute %d first eigenstructure is %.5f. \n", NUM_EIGEN, time_span);
	printf("Time to compute %d first eigenstructure (in MATLAB ONLY) is %.5f. \n", NUM_EIGEN, ts2);

	// Testing out
	//printf("Check orthogonality: \n Equal=%.5f, Not-Equal=%.5f\n", EigVec.col(0).transpose()*M*EigVec.col(0), EigVec.col(0).transpose()*M*EigVec.col(50));
	//cout << "Eigen Vectors" << endl << EigVec.block(0, 0, 2, 10) << endl;
	printf("EigenVector dimension: %dx%d\n", EigVec.rows(), EigVec.cols());
}
