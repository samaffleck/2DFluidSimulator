#include "2DFluidSimulator/ODESolver.h"
#include <iostream>


void ODESolver::initialise(int numberOfXNodes, int numberOfYNodes)
{
	nx = numberOfXNodes;
	ny = numberOfYNodes;
}


void ODESolver::solve(Eigen::MatrixXd& var,
	const Eigen::MatrixXd& Ao,
	const Eigen::MatrixXd& Ae,
	const Eigen::MatrixXd& Aw,
	const Eigen::MatrixXd& An,
	const Eigen::MatrixXd& As,
	const Eigen::MatrixXd& S)
{
	// Gauss-seidel method
	double error = getResidual(var, Ao, Ae, Aw, An, As, S);
	int itterations = 0;

	while (error > m_tolerance)
	{
		update(var, Ao, Ae, Aw, An, As, S);
		error = getResidual(var, Ao, Ae, Aw, An, As, S);
		itterations++;
	}

	std::cout << "Equation Solved Successfully. \nError : " << error << "\tInner Itterations : " << itterations << std::endl;
}


void ODESolver::update(Eigen::MatrixXd& var,
	const Eigen::MatrixXd& Ao,
	const Eigen::MatrixXd& Ae,
	const Eigen::MatrixXd& Aw,
	const Eigen::MatrixXd& An,
	const Eigen::MatrixXd& As,
	const Eigen::MatrixXd& S)
{
	for (int x = 0; x < nx; ++x)
	{
		for (int y = 0; y < ny; ++y)
		{
			setNeighbourCells(var, x, y);
			var(x, y) = (-1 / Ao(x, y)) * (Ae(x, y) * ve + Aw(x, y) * vw + An(x, y) * vn + As(x, y) * vs - S(x, y));
		}
	}
}


double ODESolver::getResidual(const Eigen::MatrixXd& var,
	const Eigen::MatrixXd& Ao,
	const Eigen::MatrixXd& Ae,
	const Eigen::MatrixXd& Aw,
	const Eigen::MatrixXd& An,
	const Eigen::MatrixXd& As,
	const Eigen::MatrixXd& S)
{
	auto res = var; // Create a copy to write the error to

	for (int x = 0; x < nx; ++x)
	{
		for (int y = 0; y < ny; ++y)
		{
			setNeighbourCells(var, x, y);
			res(x, y) = Ao(x, y) * var(x, y) + Ae(x, y) * ve + Aw(x, y) * vw + An(x, y) * vn + As(x, y) * vs - S(x, y);
		}
	}

	return res.norm();
}


void ODESolver::setNeighbourCells(const Eigen::MatrixXd& var, int x, int y)
{
	ve = 0.0;
	vw = 0.0;
	vs = 0.0;
	vn = 0.0;

	if (x < nx - 1) // not at the right boundary
	{
		ve = var(x + 1, y);
	}

	if (x > 0) // not at the left boundary
	{
		vw = var(x - 1, y);
	}

	if (y < ny - 1) // not at the top boundary
	{
		ve = var(x, y + 1);
	}

	if (y > 0) // not at the bottom boundary
	{
		vw = var(x, y - 1);
	}
}
