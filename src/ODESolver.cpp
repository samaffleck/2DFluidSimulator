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
	const Eigen::MatrixXd& S,
	double& normalisedResidual,
	double tolerance,
	double damping_factor)
{
	// Gauss-seidel method
	auto res = var;
	double absoluteResidual = tolerance * 10; //getResidualVector(var, Ao, Ae, Aw, An, As, S, res);
	//double maximumResidual = absoluteResidual;
	//normalisedResidual = absoluteResidual / (maximumResidual + 1e-15); // = 1 so it gets us to enter the while loop unless tol = 1
	m_damping_factor = damping_factor;

	int itterations = 0;
	int maxItterations = 10000;

	while (absoluteResidual > tolerance && itterations < maxItterations)
	{
		update(var, Ao, Ae, Aw, An, As, S);
		absoluteResidual = getResidualVector(var, Ao, Ae, Aw, An, As, S, res);
		//if (absoluteResidual > maximumResidual)
		//{
		//	maximumResidual = absoluteResidual;
		//}
		//normalisedResidual = absoluteResidual / (maximumResidual + 1e-15);
		itterations++;
	}

	std::cout << "Equation Solved Successfully. \nNormalised Residual : " << absoluteResidual << "\tInner Itterations : " << itterations << std::endl;
}

void ODESolver::solveInCorrectionForm(Eigen::MatrixXd& var, const Eigen::MatrixXd& Ao, const Eigen::MatrixXd& Ae, const Eigen::MatrixXd& Aw, const Eigen::MatrixXd& An, const Eigen::MatrixXd& As, const Eigen::MatrixXd& S, double& normalisedResidual, double tolerance, double damping_factor)
{
	// Gauss-seidel method
	auto res = var;
	double absoluteResidual = getResidualVector(var, Ao, Ae, Aw, An, As, S, res);
	double maximumResidual = absoluteResidual;
	normalisedResidual = absoluteResidual / (maximumResidual + 1e-15); // = 1 so it gets us to enter the while loop unless tol = 1
	m_damping_factor = damping_factor;

	int itterations = 0;
	int maxItterations = 10000;

	auto var_c = var;
	var_c.setConstant(0.0);

	while (normalisedResidual > tolerance && itterations < maxItterations)
	{
		update(var_c, Ao, Ae, Aw, An, As, res);
		var += var_c;
		absoluteResidual = getResidualVector(var, Ao, Ae, Aw, An, As, S, res);
		if (absoluteResidual > maximumResidual)
		{
			maximumResidual = absoluteResidual;
		}
		assert(maximumResidual != 0);
		normalisedResidual = absoluteResidual / (maximumResidual + 1e-15);
		itterations++;
	}

	std::cout << "Equation Solved Successfully. \nNormalised Residual : " << normalisedResidual << "\tInner Itterations : " << itterations << std::endl;
}


void ODESolver::update(Eigen::MatrixXd& var,
	const Eigen::MatrixXd& Ao,
	const Eigen::MatrixXd& Ae,
	const Eigen::MatrixXd& Aw,
	const Eigen::MatrixXd& An,
	const Eigen::MatrixXd& As,
	const Eigen::MatrixXd& S)
{
	auto var_old = var;

	// SWEEP: LEFT TO RIGHT, TOP TO BOTTOM
	for (int y = ny - 1; y >= 0; y--)
	{
		for (int x = 0; x < nx; ++x)
		{
			setNeighbourCells(var, x, y);
			var(x, y) = (-1 / ((1 + m_damping_factor) * Ao(x, y))) * (Ae(x, y) * ve + Aw(x, y) * vw + An(x, y) * vn + As(x, y) * vs - S(x, y) - m_damping_factor * Ao(x, y) * var_old(x, y));
		}
	}
}


double ODESolver::getResidualVector(const Eigen::MatrixXd& var,
	const Eigen::MatrixXd& Ao,
	const Eigen::MatrixXd& Ae,
	const Eigen::MatrixXd& Aw,
	const Eigen::MatrixXd& An,
	const Eigen::MatrixXd& As,
	const Eigen::MatrixXd& S,
	Eigen::MatrixXd& res)
{
	for (int x = 0; x < nx; ++x)
	{
		for (int y = 0; y < ny; ++y)
		{
			setNeighbourCells(var, x, y);
			//res(x, y) = Ao(x, y) * var(x, y) + Ae(x, y) * ve + Aw(x, y) * vw + An(x, y) * vn + As(x, y) * vs - S(x, y);
			res(x, y) = S(x, y) - Ao(x, y) * var(x, y) - Ae(x, y) * ve - Aw(x, y) * vw - An(x, y) * vn - As(x, y) * vs;
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
		vn = var(x, y + 1);
	}

	if (y > 0) // not at the bottom boundary
	{
		vs = var(x, y - 1);
	}
}
