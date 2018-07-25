#include "Basis.h"

/* [DIJKSTRA-based BASIS] */
// Computing Dijkstra Distance (sing input, D initialized inside, stops when reaching a certain distance, 
// and used for Basis Construction), COMPACT Version
void ComputeDijkstraCompact(const Eigen::MatrixXd &V, const int s1, const vector<set<int>> &AdM,
	Eigen::VectorXd &D, const double a, const double b, const double nDist, const int col, const int sampleSize,
	vector<Eigen::Triplet<double>> &UTriplet, int &elementCounter)
{
	vector<VertexPair> container;
	container.reserve(V.rows());
	priority_queue<VertexPair, std::vector<VertexPair>, std::greater<VertexPair>> DistPQueue;
	
	int elCounter = 1;

	// Computing distance for initial sample points S
	D(s1) = 0.0f;
	VertexPair vp{ s1,D(s1) };
	DistPQueue.push(vp);

	double distFromCenter = numeric_limits<double>::infinity();

	// For other vertices in mesh
	do {
		if (DistPQueue.size() == 0) break;
		VertexPair vp1 = DistPQueue.top();
		distFromCenter = vp1.distance;
		DistPQueue.pop();

		// Updating the distance for neighbors of vertex of lowest distance in priority queue
		auto const& elem = AdM[vp1.vId];
		for (auto it = elem.begin(); it != elem.end(); ++it) {
			double dist;
			/* Regular Dikjstra */
			//VtoVDist(V.row(vp1.vId), V.row(*it), dist);
			//double tempDist = vp1.distance + dist;

			/* Correct distance using Euclidean Distance */
			double tempDist = (V.row(s1) - V.row(*it)).norm();

			/* updating the distance */
			if (tempDist < D(*it)) {
				D(*it) = tempDist;
				VertexPair vp2{ *it,tempDist };
				DistPQueue.push(vp2);
			}
		}
	} while (distFromCenter <= nDist);

	// Constructing a sparse matrix based on the triplet

	/*Parameters for 5th-Order Polynomial */
	double a5 = -6.0 / (nDist*nDist*nDist*nDist*nDist);
	double b5 = 15.0 / (nDist*nDist*nDist*nDist);
	double c5 = -10.0 / (nDist*nDist*nDist);

	/* Gaussian Mapping */
	double tHalf = 0.45; tHalf = tHalf*tHalf;
	double yHalf = 0.5; yHalf = log(yHalf);
	double aCoeff = yHalf / (nDist*nDist*tHalf);

	for (int i = 0; i < D.size(); i++) {
		double aVal = D(i);
		if (aVal < nDist) {
			/* Linear Mapping */
			//double val = 1.0 - aVal/nDist;
			//UTriplet.push_back(Eigen::Triplet<double>(i, col, val));

			/* Polynomial mapping (3rd Order) */
			double val = a*aVal*aVal*aVal + b*aVal*aVal + 1.0;
			UTriplet.push_back(Eigen::Triplet<double>(i, col, val));

			/* Polynomial mapping (5th Order) */
			//double val = a5*aVal*aVal*aVal*aVal*aVal + b5*aVal*aVal*aVal*aVal + c5*aVal*aVal*aVal + 1.0;
			//UTriplet.push_back(Eigen::Triplet<double>(i, col, val));

			/* Gaussian Mapping */
			//double val = exp(aCoeff*aVal*aVal);
			//UTriplet.push_back(Eigen::Triplet<double>(i, col, val));
		}

		// For the next samples, reset to infinity
		D(i) = numeric_limits<double>::infinity();
	}
}