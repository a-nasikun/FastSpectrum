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
		//for (set<int>::iterator it = AdM[vp1.vId].begin(); it != AdM[vp1.vId].end(); ++it) {
		for (auto it = elem.begin(); it != elem.end(); ++it) {
			double dist;
			/* Dikjstra: Compute on the fly/online */
			//VtoVDist(V.row(vp1.vId), V.row(*it), dist);
			//double tempDist = vp1.distance + dist;

			/* Euclidean Distance */
			double tempDist = (V.row(s1) - V.row(*it)).norm();
			if (tempDist < D(*it)) {
				D(*it) = tempDist;
				VertexPair vp2{ *it,tempDist };
				{
					DistPQueue.push(vp2);
				}
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

/* [STVD-based BASIS] */
// Get angle1 and angle2 (total not 360) between p1 and p2, w.r.t vertex v1
void ComputeAngle4STVD(igl::viewer::Viewer &viewer, MyMesh &mesh, const Eigen::MatrixXd &V, int &p1, int &p2, const int &v1, double &totAngle1, double &totAngle2)
{
	MyMesh::VertexIter          v_it, v_end(mesh.vertices_end());
	MyMesh::VertexVertexIter    vv_it;
	MyMesh::Point               cog;
	int							hit1 = 0;
	vector<int>					neighbors;
	double						tempTot = 0.0;

	totAngle1			= 0.0, totAngle2 = 0.0;
	bool flagAngle1		= false;
	int hitID			= 0;

	// Iterate through v1
	for (vv_it = mesh.vv_iter(mesh.vertex_handle(v1)); vv_it.is_valid(); ++vv_it){
		neighbors.push_back((*vv_it).idx());
	}
	
	/* Getting the angles between vertices */
	for (int i = 0; i < neighbors.size(); i++) {
		if (neighbors[i] == p1) { hitID = 1; }
		else if (neighbors[i] == p2) {hitID = 2;}

		double locAngle;
		if (i == (neighbors.size() - 1)) {
			locAngle = VtoVangle(V.row(neighbors[i]), V.row(neighbors[0]), V.row(v1));
		} else {
			locAngle = VtoVangle(V.row(neighbors[i]), V.row(neighbors[i + 1]), V.row(v1));
		}
		
		/* Summing the angles */
		if (hitID == 0) {
			tempTot += locAngle;
		} else if (hitID == 1) {
			totAngle1 += locAngle;
		} else if (hitID == 2) {
			totAngle2 += locAngle;
		}
	}

	if (hitID == 1) totAngle1 += tempTot;
	else if (hitID == 2) totAngle2 += tempTot;
}

// Unfolding the edges, K times backward, unless reaching the Si (initial point) & distance in form of VECTORS
void UnfoldEdgeSTVD(igl::viewer::Viewer &viewer, MyMesh &mesh, const Eigen::MatrixXd &V, const int Si, const vector<set<int>> &AdM,
	const Eigen::SparseMatrix<double> &DistanceTableSpM, const Eigen::VectorXi &prev, const int &initNode, const int &max_i, vector<Eigen::RowVector2d>	&EVect)
{
	/* Initiating variables */
	int							curNode = initNode;
	int							p1 = curNode, p2, v1;
	int							pCounter = 0;
	double						angle1, angle2;
	//vector<Eigen::RowVector2d>	EVect;
	Eigen::RowVector2d			e, E(0.0, 0.0), ENorm;
	Eigen::Vector3d				VectorDijkstra(0, 0, 0);
	double						eTot_distance = 0.0;

	double totDistance = 0.0;

	do {
		//viewer.data.add_edges(V.row(curNode), V.row(prev(curNode)), Eigen::RowVector3d(0.1, 0.1, 0.9));
		VectorDijkstra += (V.row(curNode) - V.row(prev(curNode)));

		if (pCounter == 0) {	// For the first edge, just make it into (e, 0.0)
			e(0) = DistanceTableSpM.coeff(curNode, prev(curNode));
			e(1) = 0.0;
			EVect.push_back(e);
			Eigen::RowVector3d		e_(e(0), e(1), 0.0);
		}
		else {
			v1 = curNode;
			p2 = prev(curNode);
			ComputeAngle4STVD(viewer, mesh, V, p1, p2, v1, angle1, angle2);
			//Eigen::RowVector2d		e1(-1.0*DistanceTableSpM.coeff(v1, p2), 0.0);
			Eigen::RowVector2d		e1 = EVect[pCounter - 1];
			double rotAngle = 2.0*M_PI*(angle1 / (angle1 + angle2));
			rotateVector2D(-e1, rotAngle, e);
			e = DistanceTableSpM.coeff(v1, p2) / e1.norm() * e;
			Eigen::RowVector3d		e_(e(0), e(1), 0.0);
			EVect.push_back(e);
			p1 = curNode;			// NEXT iteration: p1<=curNode, 
		}

		curNode = prev(curNode);
		pCounter++;
	} while (curNode != Si && EVect.size()<max_i);

	/* Normalizing E */
	for (int i = 0; i < EVect.size(); i++) {
		E += EVect[i];
	}

	ENorm = E / E.norm();
}


// Return the minimum distance after Unfolding the edges of STVD
double MinDistance4STVD(igl::viewer::Viewer &viewer, MyMesh &mesh, const Eigen::MatrixXd &V, const int Si, const vector<set<int>> &AdM,
	const Eigen::SparseMatrix<double> &DistanceTableSpM, Eigen::VectorXd &D, Eigen::VectorXi &prev, const int &u, const int &w, const int &max_pred)
{	
	int tempNode = prev(w);
	prev(w) = u;
	double minDist = numeric_limits<double>::infinity();
	int					max_iter = max_pred;

	vector<int>			nodePrev;
	int					curNode = w;

	/* Save the predecessor until (k+1) or when it reaches the origin */
	do {
		nodePrev.push_back(prev(curNode));
		curNode = prev(curNode);
	} while (curNode != -2 && nodePrev.size() <= (max_pred));

	/* If max_pred is larger than possible values (nodePrev size), truncate it */
	if (nodePrev.size() < max_pred)
		max_iter = nodePrev.size();

	/* Run the unfolding only ONCE (on longest trackback) and use it over and over again */
	vector<Eigen::RowVector2d>	EVect;
	EVect.reserve(max_iter);
	UnfoldEdgeSTVD(viewer, mesh, V, Si, AdM, DistanceTableSpM, prev, w, max_iter, EVect);

	/* Computing the Unfolded distance */
	for (int i = 1; i <= max_iter; i++) {
		if (nodePrev[i-1] == -2) break;
		double tempDist;
		//printf("[%d] prevNode %d has distance %.3f\n", i, nodePrev[i], tempDist);
		if (i == 1) {		// For pred=1, it's similar to regular Dijkstra distance
			tempDist = D(u) + DistanceTableSpM.coeff(u, w);
		}
		else {				// For pred>1, unfold the edges			
			int prevNode = nodePrev[i - 1];
			
			/* Unfolding ONCE and called several times */
			double normDist;
			{
				Eigen::RowVector2d E(0, 0), ENorm;
				for (int j = 0; j < i; j++) {
					E += EVect[j];
				}
				ENorm = E / E.norm();
				double totDistance = 0.0;
				for (int j = 0; j < i; j++) {
					double locDist = EVect[j].dot(ENorm);
					totDistance += abs(locDist);
				}
				normDist = totDistance;
			}
			tempDist = D(prevNode) + normDist;
			//printf("prevNode=%d, normDist=%.3f, tempDist=%.5f\n", prevNode, normDist, tempDist);
		}

		/* Get the minimum distance */
		if (tempDist < minDist) {
			minDist = tempDist;
		}

	}
	prev(w) = tempNode;
	return minDist;
}


// STVD for Basis Construction
void ComputeSTVDBasis(igl::viewer::Viewer &viewer, MyMesh &mesh, const Eigen::MatrixXd &V, const int Si, const vector<set<int>> &AdM, const Eigen::SparseMatrix<double> &DistanceTableSpM,
	Eigen::VectorXd &D, const double a, const double b, const double nDist, const int col, const int sampleSize,
	vector<Eigen::Triplet<double>> &UTriplet)
{
	// Define containers for (1) visited nodes and (2) priority queue
	set<int>					vVisited;
	//Eigen::VectorXd				D;
	Eigen::SparseVector<double> spV;
	Eigen::VectorXi				prev;
	int							insertCounter = 0;
	const int					MAX_K = 3;
	priority_queue<VertexPair, std::vector<VertexPair>, std::greater<VertexPair>> DistPQueue;

	spV.resize(V.rows());
	spV.reserve(V.rows() / sampleSize * 30);

	// Initiate value of Prev(v) with -1
	prev.resize(V.rows());
	for (int i = 0; i < V.rows(); i++) {
		prev(i) = -1;
	}

	// Computing distance for initial sample points S		
	D(Si)			= 0.0f;
	prev(Si)		= -2;
	VertexPair		vp{ Si,D(Si) };
	DistPQueue.push(vp);
	insertCounter++;

	double distFromCenter = numeric_limits<double>::infinity();

	// For other vertices in mesh
	do
	{
		VertexPair			vp1 = DistPQueue.top();
		Eigen::VectorXd		vRow1 = V.block<1, 3>(vp1.vId, 0);
		distFromCenter		= vp1.distance;
		DistPQueue.pop();

		for (set<int>::iterator it = AdM[vp1.vId].begin(); it != AdM[vp1.vId].end(); ++it) {
			double dist;

			/* The real STVD */
			double tempDist = MinDistance4STVD(viewer, mesh, V, Si, AdM, DistanceTableSpM, D, prev, vp1.vId, *it, MAX_K);
			if (tempDist < D(*it)) {
				D(*it)		= tempDist;
				prev(*it)	= vp1.vId;
				VertexPair	vp2{ *it,tempDist };
				DistPQueue.push(vp2);
				insertCounter++;
			}
		}
		vVisited.insert(vp1.vId);
	} while (distFromCenter <= nDist && DistPQueue.size()>0);

	// Constructing a sparse matrix based on the triplet
	//int counter = 0;		/* Gaussian Mapping */
	double tHalf = 0.3; tHalf = tHalf*tHalf;
	double yHalf = 0.5; yHalf = log(yHalf);
	double aCoeff = yHalf / (nDist*nDist*tHalf);

	//printf("Max neigh=%.3f; a=%.3f; tHalf=%.3f; yHalf=%.3f\n", nDist, aCoeff, tHalf, yHalf);
	//for (Eigen::SparseVector<double>::InnerIterator it(spV); it; ++it) {
	for (int i = 0; i<D.size(); i++) {
		double aVal = D(i);

		if (aVal <= nDist) {
			/* Polynomial mapping*/
			double val = a*aVal*aVal*aVal + b*aVal*aVal + 1;
			UTriplet.push_back(Eigen::Triplet<double>(i, col, val));	

			D(i) = numeric_limits<double>::infinity();
		}
	}
}
