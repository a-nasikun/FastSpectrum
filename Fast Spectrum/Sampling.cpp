#include "Sampling.h"

/* Construct random n number of samples */
void constructRandomSample(Eigen::VectorXi &Sample, const Eigen::MatrixXd &V, int n)
{
	std::set<int>						sample;
	std::random_device					rd;
	std::default_random_engine			generator(rd());
	std::uniform_int_distribution<long> distribution(0, V.rows());

	do {
		int newSample = distribution(generator);
		sample.insert(newSample);
	} while (sample.size() < n);

	Sample.resize(n);

	int i = 0;
	for (auto s : sample) {
		Sample(i++) = s;
	}
}

/* Construct Voxel-based sampling */
void constructVoxelSample(igl::opengl::glfw::Viewer &viewer, const Eigen::MatrixXd &V, const int &n, Eigen::VectorXi &Sample)
{
	// Initialization of computing boundary of each box	
	Eigen::Vector3d minV(V.col(0).minCoeff(), V.col(1).minCoeff(), V.col(2).minCoeff());
	Eigen::Vector3d maxV(V.col(0).maxCoeff(), V.col(1).maxCoeff(), V.col(2).maxCoeff());
	Eigen::Vector3d range, boxDist;

	range = maxV - minV;
	boxDist = range / (double)n;
	set<int>		locSample;


	// Constructing tensor (n x n x n x (box-members)) for sample candidates
	vector<int>							v1;
	vector<vector<int>>					v2;
	vector<vector<vector<int>>>			v3;
	vector<vector<vector<vector<int>>>> v4;

	// Creating empty vectors
	for (int i = 0; i < n; i++) { v2.push_back(v1); }
	for (int i = 0; i < n; i++) { v3.push_back(v2); }
	for (int i = 0; i < n; i++) { v4.push_back(v3); }


	//Put each vertex to its corresponding box
	for (int i = 0; i < V.rows(); i++) {
		int xId, yId, zId;
		xId = (int)abs(floor((V(i, 0) - minV(0)) / boxDist(0)));
		yId = (int)abs(floor((V(i, 1) - minV(1)) / boxDist(1)));
		zId = (int)abs(floor((V(i, 2) - minV(2)) / boxDist(2)));
		if (xId >(n - 1)) xId = n - 1;
		if (yId >(n - 1)) yId = n - 1;
		if (zId > (n - 1)) zId = n - 1;
		v4[xId][yId][zId].push_back(i);
	}

	// Populate the boxes of samples
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			for (int k = 0; k < n; k++) {
				// Only cares for non-empty boxes
				if (v4[i][j][k].size() > 0) {
					Eigen::Vector3d boxMin(minV(0) + i*boxDist(0), minV(1) + j*boxDist(1), minV(2) + k*boxDist(2));

					double	lambda = 3.0;
					double	minVal = numeric_limits<double>::infinity();
					int		minDistId;

					// Selecting vertex that closest to the center of each box/container;
					for (int l = 0; l < v4[i][j][k].size(); l++) {
						double	dist;
						VtoVDist(V.row(v4[i][j][k][l]), (boxMin + boxDist / 2.0), dist);
						//double pDist = lambda * exp(-lambda*dist);
						if (dist < minVal) {
							minVal = dist;
							minDistId = v4[i][j][k][l];
						}
					}
					locSample.insert(minDistId);
				}
			}
		}
	}
	Sample.resize(locSample.size());
	int i = 0;
	for (int s : locSample) {
		Sample(i++) = s;
		//viewer.data.add_points(V.row(s), Eigen::Vector3d(1.0, 0.0, 0.0).transpose());
	}
	printf("%d-box produces %d samples.\n", n, locSample.size());
}

/* Construct Poisson-Disk Sample */
void constructPoissonDiskSample(const Eigen::MatrixXd &V, const int &numSamples, const int &n, const double &avgEdgeLength, Eigen::VectorXi &Sample)
{
	/* Initialization of computing boundary of each box	 */
	Eigen::Vector3d minV(V.col(0).minCoeff(), V.col(1).minCoeff(), V.col(2).minCoeff());
	Eigen::Vector3d maxV(V.col(0).maxCoeff(), V.col(1).maxCoeff(), V.col(2).maxCoeff());
	Eigen::Vector3d range, boxDist;
	Eigen::Vector3i nn;

	range = maxV - minV;
	double minRange = range.minCoeff();

	for (int i = 0; i < 3; i++) {
		boxDist(i) = minRange / (double) n;
		nn(i) = (int)ceil(range(i) / boxDist(i));
		printf("nn(%d)=%d || boxDis(%d)=%.3f. \n", i, nn(i), i, boxDist(i));
	}

	//boxDist = range / (double)n;

	set<int>		locSample;
	const int		SamplePerBox = 10;
	//double			radius = 0.50 * 1.0 / 3.0 * (boxDist(0) + boxDist(1) + boxDist(2));
	double			radius = sqrt((double)V.rows() / (2.5*(double)numSamples)) * avgEdgeLength;

	/* Container for boxes */
	priority_queue<BoxStruct, std::vector<BoxStruct>, std::greater<BoxStruct>>	UnvisitedBoxQueue;
	list<int>																	VisitedEmptyBoxList;

	/* Constructing tensor (n x n x n x (box-members)) for sample candidates */
	list<int>							v1;
	vector<list<int>>					v2;
	vector<vector<list<int>>>			v3;
	vector<vector<vector<list<int>>>>	v4;

	/* Instatiating the tensor with empty elements */
	v2.reserve(nn(2));	v3.reserve(nn(1));	v4.reserve(nn(0));
	for (int i = 0; i < nn(2); i++) { v2.push_back(v1); }
	for (int i = 0; i < nn(1); i++) { v3.push_back(v2); }
	for (int i = 0; i < nn(0); i++) { v4.push_back(v3); }
	
	printf("v2 size=%d, v3=%d, v4=%d \n", v2.size(), v3.size(), v4.size());
	
	/* Put each vertex to its corresponding box */
	for (int vId = 0; vId < V.rows(); vId++) {
		int xId, yId, zId;
		xId = (int)abs(floor((V(vId, 0) - minV(0)) / boxDist(0)));
		yId = (int)abs(floor((V(vId, 1) - minV(1)) / boxDist(1)));
		zId = (int)abs(floor((V(vId, 2) - minV(2)) / boxDist(2)));
		if (xId >(nn(0) - 1)) xId = nn(0) - 1;
		if (yId >(nn(1) - 1)) yId = nn(1) - 1;
		if (zId >(nn(2) - 1)) zId = nn(2) - 1;		

		v4[xId][yId][zId].push_back(vId);
	}

	// Randomizer for large number (LONG)
	std::random_device					rd;
	std::default_random_engine			generator(rd());
	int									trimCounter = 0;

	/* FOR TESTING ONLY*/
	int poolCounter = 0;
	/* END OF TESTING */

	// Populating Boxes with 10 random vertices
	for (int iIdx = 0; iIdx < nn(0); iIdx++) {
		for (int jIdx = 0; jIdx < nn(1); jIdx++) {
			for (int kIdx = 0; kIdx < nn(2); kIdx++) {
				Eigen::Vector3d		boxMin(minV(0) + iIdx*boxDist(0), minV(1) + jIdx*boxDist(1), minV(2) + kIdx*boxDist(2));

				/* IF there are more samples on a box then required, select randomly 10 of them; */
				int boxSize = v4[iIdx][jIdx][kIdx].size();
				if (boxSize == 0) {
					//VisitedEmptyBoxList.push_back(iIdx*nn(0)*nn(0) + jIdx*nn(1) + kIdx);
					VisitedEmptyBoxList.push_back(iIdx*nn(2)*nn(1) + jIdx*nn(2) + kIdx);
				}
				else
				{	// Work on NON-EMPTY boxes only
					double randDist = (double)(rand() % 100) / (double)10;
					//BoxStruct curBox{ iIdx*nn(0)*nn(0) + jIdx*nn(1) + kIdx, randDist };
					BoxStruct curBox{ iIdx*nn(2)*nn(1) + jIdx*nn(2) + kIdx, randDist };
					UnvisitedBoxQueue.push(curBox);

					/* If there are more elements than required, randomly discard the rest */
					if (boxSize > SamplePerBox) {
						set<int> vertexToDel;
						srand(time(NULL));
						do {
							std::uniform_int_distribution<long> distribution(0, boxSize);
							int vSel = distribution(generator);
							vertexToDel.insert(vertexToDel.end(), vSel);
						} while (vertexToDel.size() <= (boxSize - SamplePerBox));

						int l = 0;
						for (list<int>::iterator it = v4[iIdx][jIdx][kIdx].begin(); it != v4[iIdx][jIdx][kIdx].end(); ++it) {
							if (vertexToDel.find(l) != vertexToDel.end()) {
								v4[iIdx][jIdx][kIdx].remove(*it);
							}
							l++;
						}
						trimCounter++;
					} // End of vertex(-ices) deletion
				}
				/* FOR TESTING ONLY*/
				poolCounter += v4[iIdx][jIdx][kIdx].size();
				/* END OF TESTING */
			}
		}
	} // End of Populating the boxes with samples POOL
	
	printf("Has set up %d samplepool. \n", poolCounter);


	/* Getting samples from POOL */
	while (VisitedEmptyBoxList.size()<nn(0)*nn(1)*nn(2)) {
	//while(!UnvisitedBoxQueue.empty()){
		/* Randomly select a certain box */
		BoxStruct box1 = UnvisitedBoxQueue.top();
		UnvisitedBoxQueue.pop();

		/* Determining in which box is this point located */
		int b_xId, b_yId, b_zId;
		b_xId = round(box1.id / (nn(1)*nn(2)));
		b_yId = round((box1.id - b_xId*nn(1)*nn(2)) / nn(2));
		b_zId = box1.id % nn(2);
		//if (box1.id % 100 == 70) 
			//printf("%d => (%d, %d, %d): SIZE=%d.\n", box1.id, b_xId, b_yId, b_zId, VisitedEmptyBoxList.size());
		if (v4[b_xId][b_yId][b_zId].size() < 1) continue;
		
		
		/* Picking a random sample from a box */
		srand(time(NULL));
		int randV1 = rand() % v4[b_xId][b_yId][b_zId].size();
		int boxSample;
		int locCounter = 0;
		for (list<int>::iterator it = v4[b_xId][b_yId][b_zId].begin(); it != v4[b_xId][b_yId][b_zId].end(); ++it) {
			if (randV1 == locCounter) boxSample = *it;
			locCounter++;
		}
		locSample.insert(locSample.end(), boxSample);
		
		/* Iterating through neigboring boxes */
		for (int i = b_xId - 1; i <= b_xId + 1; i++) {
			if (i < 0 || i >(nn(0) - 1)) continue;
			for (int j = b_yId - 1; j <= b_yId + 1; j++) {
				if (j < 0 || j >(nn(1) - 1)) continue;
				for (int k = b_zId - 1; k <= b_zId + 1; k++) {
					if (k < 0 || k >(nn(2) - 1)) continue;
					if (v4[i][j][k].size() > 0) {
						/* Removing samples closer than a certain radius */
						for (list<int>::iterator it = v4[i][j][k].begin(); it != v4[i][j][k].end(); ++it) {
							double distAB;
							VtoVDist(V.row(boxSample), V.row(*it), distAB);
		
							if (distAB <= radius) {
								v4[i][j][k].remove(*it);
							}
						} /* End of Removing samples closer than a certain radius */
		
						  /* Put the empty boxes to its list */
						if (v4[i][j][k].size() < 1) {
							VisitedEmptyBoxList.push_back(i*nn(1)*nn(2) + j*nn(2) + k);
						}
						else
						{
							/* If the original box is not empty, put it back to queue with new priority */
							if (i == b_xId && j == b_yId && k == b_zId) {
								double nRandDist = (double)(rand() % 100) / (double)10;
								// ID = a*(y*z) + b*(z) + c.
								BoxStruct curBox{ i*nn(1)*nn(2) + j*nn(2) + k, nRandDist };
								UnvisitedBoxQueue.push(curBox);
							}
						}
					}
				} /* End of iteration on each box */
			}
		} /* End of Iterating through neigboring boxes */
	}
	/* END of Iterating to get All samples */

	/* Passing sample to caller */
	Sample.resize(locSample.size());
	int sId = 0;
	for (int samp : locSample) {
		Sample(sId++) = samp;
	}

	printf("Sample: %d from %d pool (%.1f percent deletion) (boxes=%d).\n", Sample.size(), poolCounter, 100.0*(double)(poolCounter - Sample.size()) / (double)poolCounter, nn(0)*nn(1)*nn(2));
}

void ComputeDijkstra(const Eigen::MatrixXd &V, const int Si, const vector<set<int>> &AdM, Eigen::VectorXd &D)
{
	// Define containers for (1) visited nodes and (2) priority queue
	priority_queue<VertexPair, std::vector<VertexPair>, std::greater<VertexPair>> DistPQueue;
	int insertCounter = 0;

	// Computing distance for initial sample points S
	Eigen::VectorXd vRow1 = V.block<1, 3>(Si, 0);
	D(Si) = 0.0f;
	VertexPair vp{ Si,D(Si) };
	DistPQueue.push(vp);
	insertCounter++;

	double distFromCenter = numeric_limits<double>::infinity();
	// For other vertices in mesh
	while (DistPQueue.size()>0)
	{
		VertexPair vp1 = DistPQueue.top();
		Eigen::VectorXd vRow1 = V.block<1, 3>(vp1.vId, 0);
		DistPQueue.pop();

		for (set<int>::iterator it = AdM[vp1.vId].begin(); it != AdM[vp1.vId].end(); ++it) {
			double dist;
			Eigen::VectorXd vRow2 = V.block<1, 3>(*it, 0);
			VtoVDist(vRow1, vRow2, dist);
			double tempDist = vp1.distance + dist;
			if (tempDist < D(*it)) {
				D(*it) = tempDist;
				distFromCenter = tempDist;
				VertexPair vp2{ *it,tempDist };
				DistPQueue.push(vp2);
				insertCounter++;
			}
		}
	}
}

void constructFarthestPointSample(const Eigen::MatrixXd &V, const std::vector<std::set<int>> &AdM, const int &numSamples, Eigen::VectorXi &Samples)
{
	Samples.resize(numSamples);
	Eigen::VectorXd D;
	Eigen::VectorXi Si;
	D.resize(V.rows());

	// Initialize the value of D
	for (int i = 0; i < V.rows(); i++) {
		D(i) = numeric_limits<double>::infinity();
	}

	srand(time(NULL));
	Samples(0) = rand() % V.rows();
	Si = Samples.head(1);

	ComputeDijkstra(V, Samples(0), AdM, D);
	Eigen::VectorXi::Index maxIndex1;
	D.maxCoeff(&maxIndex1);
	Samples(0) = maxIndex1;

	for (int i = 1; i < numSamples; i++) {
		Eigen::VectorXi::Index maxIndex;
		ComputeDijkstra(V, Samples(i - 1), AdM, D);
		D.maxCoeff(&maxIndex);
		Samples(i) = maxIndex;
	}
}