#include "2DFluidSimulator/ODESolver.h"
#include <iostream>


void setNeighbourCells(double& Ve,
	double& Vw,
	double& Vn,
	double& Vs,
	const Eigen::VectorXd& V,
	int n,
	int k,
	int N)
{
	double leftBoundary = 100;
	double rightBoundary = 0;
	double topBoundary = 0;
	double bottomBoundary = 100;

	if ((n + 1) % k != 0) // Check right boundary
	{
		Ve = V[n + 1];
	}
	else // Apply right boundary condition
	{
		Ve = 2 * rightBoundary - V[n - 1];
	}

	if (n % k != 0) // Check left boundary
	{
		Vw = V[n - 1];
	}
	else // Apply left boundary condition
	{
		Vw = 2 * leftBoundary - V[n + 1];
	}

	if (n + k < N) // Check top boundary
	{
		Vn = V[n + k];
	}
	else // Apply top boundary condition
	{
		Vn = 2 * topBoundary - V[n - k];
	}

	if (n - k >= 0) // Check bottom boundary
	{
		Vs = V[n - k];
	}
	else // Apply bottom boundary condition
	{
		Vs = 2 * bottomBoundary - V[n + k];
	}
}

void ODESolver::initialise(int numberOfNodes)
{
	k = numberOfNodes;
}

void ODESolver::solve(Eigen::VectorXd& V)
{
	// Gauss-seidel method
	double error = getResudule(V);
	int itterations = 0;

	while (error > m_tolerance)
	{
		update(V);
		error = getResudule(V);
		itterations++;
	}

	std::cout << "Equation Solved Successfully. \nError : " << error << "\tItterations : " << itterations << std::endl;
}

void ODESolver::update(Eigen::VectorXd& V)
{
	auto N = (int)V.size();

	double Ve = 0.0;
	double Vw = 0.0;
	double Vn = 0.0;
	double Vs = 0.0;

	for (int n = 0; n < N; ++n)
	{
		setNeighbourCells(Ve, Vw, Vn, Vs, V, n, k, N);
		V[n] = (-1 / Ao[n]) * (Ae[n] * Ve + Aw[n] * Vw + An[n] * Vn + As[n] * Vs - S[n]);
	}
}

double ODESolver::getResudule(const Eigen::VectorXd& V)
{
	auto res = V; // Create a copy to write the error to
	auto N = (int)V.size();

	double Ve = 0.0;
	double Vw = 0.0;
	double Vn = 0.0;
	double Vs = 0.0;

	for (int n = 0; n < N; ++n)
	{
		setNeighbourCells(Ve, Vw, Vn, Vs, V, n, k, N);
		res[n] = V[n] * Ao[n] + Ae[n] * Ve + Aw[n] * Vw + An[n] * Vn + As[n] * Vs - S[n];
	}

	return res.norm();
}

