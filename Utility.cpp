#include "Utility.h"

/* Create a list of neighbors of every vertex *
** Input:
**		V: matrix of vertices
**		F: matrix of faces (e.g. triangles)
** Output:
**		AM: Adjacency matrix */
void initiateAdM(Eigen::MatrixXd &V, Eigen::MatrixXi &F, vector<set<int>> &AM) {
	for (int i = 0; i < V.rows(); i++) {
		set<int>		emptySet;
		AM.insert(AM.end(), emptySet);
	}

	for (int i = 0; i < F.rows(); i++)
	{
		int a1 = F(i, 0);
		int a2 = F(i, 1);
		int a3 = F(i, 2);

		AM[a1].insert(a2);
		AM[a1].insert(a3);
		AM[a2].insert(a1);
		AM[a2].insert(a3);
		AM[a3].insert(a1);
		AM[a3].insert(a2);
	}
}

/* Precompute the distance between a vertex and its neighbor(s), for a bit of speed up */
void initiateDistTableSpM(Eigen::MatrixXd &V, const vector<set<int>> &AM, Eigen::SparseMatrix<double> &DistanceTable)
{
	std::vector<Eigen::Triplet<double>> locTrip;
	locTrip.reserve(V.rows());
	DistanceTable.resize(V.rows(), V.rows());
	DistanceTable.reserve(V.rows() * 6 * 2);

	for (int i = 0; i < V.rows(); i++) {
		double distAB;
		for (std::set<int, double>::iterator it = AM[i].begin(); it != AM[i].end(); ++it) {
			distAB = (V.row(i) - V.row(*it)).norm();
			locTrip.push_back(Eigen::Triplet<double>(i, *it, distAB));
		}
	}

	DistanceTable.setFromTriplets(locTrip.begin(), locTrip.end());
}

/* Compute distance between two points */
void VtoVDist(const Eigen::VectorXd &V1, const Eigen::VectorXd &V2, double &dist)
{
	Eigen::VectorXd	  v = V2 - V1;
	dist				= v.norm();
}

/* Rotate the edge*/
void rotateVector2D(const Eigen::RowVector2d &e, const double &angle, Eigen::RowVector2d &e2)
{
	double a = cos(angle);
	double b = sin(angle);
	e2(0) = e(0)*a - e(1)*b;
	e2(1) = e(0)*b + e(1)*a;
}

/* Get the angle between two edges */
double VtoVangle(const Eigen::Vector3d &V1, const Eigen::Vector3d &V2, const Eigen::Vector3d &P)
{
	Eigen::Vector3d v1 = V1 - P;
	Eigen::Vector3d v2 = V2 - P;
	double angle_ = acos((v1.dot(v2)) / (v1.norm()*v2.norm()));
	return angle_;
}