#include "2DFluidSimulator/Equation_NavierStokes.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <filesystem>


void Equation_NavierStokes::initialiseEquation(int numberOfXCells, int numberOfYCells)
{
	nx = numberOfXCells;
	ny = numberOfYCells;
	
	double viscosity = 1e-3;
	double density = 1.0;
	double initialPressure = 0.0;
	double initialXVelocity = 0.0;
	double initialYVelocity = 0.0;

	u.setConstant(nx, ny, initialXVelocity);
	u_c.setConstant(nx, ny, 0.0);
	v.setConstant(nx, ny, initialYVelocity);
	v_c.setConstant(nx, ny, 0.0);
	p.setConstant(nx, ny, initialPressure);
	p_c.setConstant(nx, ny, 0.0);
	vis.setConstant(nx, ny, viscosity);
	rho.setConstant(nx, ny, density);
	Sp_x.setConstant(nx, ny, 0.0);
	Sp_y.setConstant(nx, ny, 0.0);

	Ao.setConstant(nx, ny, 0.0);
	Ae.setConstant(nx, ny, 0.0);
	Aw.setConstant(nx, ny, 0.0);
	An.setConstant(nx, ny, 0.0);
	As.setConstant(nx, ny, 0.0);
	
	Ap_o.setConstant(nx, ny, 0.0);
	Ap_e.setConstant(nx, ny, 0.0);
	Ap_w.setConstant(nx, ny, 0.0);
	Ap_n.setConstant(nx, ny, 0.0);
	Ap_s.setConstant(nx, ny, 0.0);
	Sp.setConstant(nx, ny, 0.0);

	vel_face.resize(nx, ny);
	vel_c_face.resize(nx, ny);
	rho_face.resize(nx, ny);
	vis_face.resize(nx, ny);
	mflux_face.resize(nx, ny);

	for (int x = 0; x < nx; ++x)
	{
		for (int y = 0; y < ny; ++y)
		{
			vis_face(x, y).east = viscosity;
			vis_face(x, y).west = viscosity;
			vis_face(x, y).north = viscosity;
			vis_face(x, y).south = viscosity;

			rho_face(x, y).east = density;
			rho_face(x, y).west = density;
			rho_face(x, y).north = density;
			rho_face(x, y).south = density;
		}
	}
}


void Equation_NavierStokes::update()
{
	int itt = 0;
	int maxItterations = 200;
	bool hasConverged = false;

	while (!hasConverged && itt < maxItterations)
	{
		// Update momentum links and sources for x and y-direction
		updateLinkCoefficient();

		// solve x-momentum
		m_solver.solve(u, Ao, Ae, Aw, An, As, Sp_x);

		// solve y-momentum
		m_solver.solve(v, Ao, Ae, Aw, An, As, Sp_y);

		// Calculate the face velocities using PWIM
		updateFaceVelocities();

		// Calculate pressure links and mass imbalance
		updatePressureLinks();

		// Solve pressure correction equation
		m_solver.solve(p_c, Ap_o, Ap_e, Ap_w, Ap_n, Ap_s, Sp);

		// Correct velocities and pressure
		updateVelocityCorrection();
		correctVelocityAndPressure();

		// Check for overall convergance
		hasConverged = isConverged();

		itt++;
		//logVelocity();
	}

	if (hasConverged)
	{
		std::cout << "Solution has converged in: \t" << itt << " itteration(s).\n";
		std::cout << "U_x residual:\t" << res_u << "\n";
		std::cout << "U_y residual:\t" << res_v << "\n";
		std::cout << "P residual:\t" << res_p << "\n";
	}
	else
	{
		std::cout << "Solution could not converge after: \t" << itt << " itteration(s).\n";
		std::cout << "Consider increasing the maximum number of itterations, or reducing your tolerance.\n";
	}
}


void Equation_NavierStokes::logVelocity()
{	
	std::cout << "x-veloctiy\n";
	for (int y = nx - 1; y >= 0; --y)
	{
		for (int x = 0; x < nx; ++x)
		{
			std::cout << u(x, y) << "\t";
		}
		std::cout << "\n";
	}
	std::cout << "\n";

	std::cout << "y-veloctiy\n";
	for (int y = nx - 1; y >= 0; --y)
	{
		for (int x = 0; x < nx; ++x)
		{
			std::cout << v(x, y) << "\t";
		}
		std::cout << "\n";
	}
	std::cout << "\n";

	std::cout << "pressure\n";
	for (int y = nx - 1; y >= 0; --y)
	{
		for (int x = 0; x < nx; ++x)
		{
			std::cout << p(x, y) << "\t";
		}
		std::cout << "\n";
	}
	std::cout << "\n";

	std::cout << "pressure correction\n";
	for (int y = nx - 1; y >= 0; --y)
	{
		for (int x = 0; x < nx; ++x)
		{
			std::cout << p_c(x, y) << "\t";
		}
		std::cout << "\n";
	}
	std::cout << "\n";

	std::cout << "Co-ordinates\n";
	for (int y = nx - 1; y >= 0; --y)
	{
		for (int x = 0; x < nx; ++x)
		{
			std::cout << "(" << x << ", " << y << ")\t";
		}
		std::cout << "\n";
	}
	std::cout << "\n";

	std::cout << "face velocity (east)\n";
	for (int y = nx - 1; y >= 0; --y)
	{
		for (int x = 0; x < nx; ++x)
		{
			std::cout << vel_face(x, y).east << "\t";
		}
		std::cout << "\n";
	}
	std::cout << "\n";

	std::cout << "face velocity (west)\n";
	for (int y = nx - 1; y >= 0; --y)
	{
		for (int x = 0; x < nx; ++x)
		{
			std::cout << vel_face(x, y).west << "\t";
		}
		std::cout << "\n";
	}
	std::cout << "\n";

	std::cout << "face velocity (north)\n";
	for (int y = nx - 1; y >= 0; --y)
	{
		for (int x = 0; x < nx; ++x)
		{
			std::cout << vel_face(x, y).north << "\t";
		}
		std::cout << "\n";
	}
	std::cout << "\n";

	std::cout << "face velocity (south)\n";
	for (int y = nx - 1; y >= 0; --y)
	{
		for (int x = 0; x < nx; ++x)
		{
			std::cout << vel_face(x, y).south << "\t";
		}
		std::cout << "\n";
	}
	std::cout << "\n";

	std::cout << "x-face-error\n";
	for (int y = nx - 2; y >= 1; --y)
	{
		for (int x = 0; x < nx - 1; ++x)
		{
			std::cout << vel_face(x, y).east - vel_face(x + 1, y).west << "\t";
		}
		std::cout << "\n";
	}
	std::cout << "\n";

	std::cout << "x-face-correction-error\n";
	for (int y = nx - 2; y >= 1; --y)
	{
		for (int x = 0; x < nx - 1; ++x)
		{
			std::cout << vel_c_face(x, y).east - vel_c_face(x + 1, y).west << "\t";
		}
		std::cout << "\n";
	}
	std::cout << "\n";

	std::cout << "y-face-error\n";
	for (int y = nx - 2; y >= 1; --y)
	{
		for (int x = 1; x < nx - 1; ++x)
		{
			std::cout << vel_face(x, y).north - vel_face(x, y + 1).south << "\t";
		}
		std::cout << "\n";
	}
	std::cout << "\n";

	std::cout << "Mass imbalance\n";
	for (int y = nx - 2; y >= 1; --y)
	{
		for (int x = 1; x < nx - 1; ++x)
		{
			std::cout << (vel_face(x, y).east - vel_face(x, y).west) + (vel_face(x, y).north - vel_face(x, y).south) << "\t";
		}
		std::cout << "\n";
	}
	std::cout << "\n";

}


void Equation_NavierStokes::logData(std::string resultsDirectory)
{
	std::filesystem::path xVelocityFilePath = std::filesystem::path(resultsDirectory) / "x-velocity.csv";
	std::filesystem::path yVelocityFilePath = std::filesystem::path(resultsDirectory) / "y-velocity.csv";
	std::filesystem::path pressureFilePath = std::filesystem::path(resultsDirectory) / "Pressure.csv";

	std::ofstream xVelocityFile(xVelocityFilePath);
	std::ofstream yVelocityFile(yVelocityFilePath);
	std::ofstream pressureFile(pressureFilePath);

	if (!xVelocityFile.is_open()) 
	{
		std::cout << "Failed to open file: " << xVelocityFilePath << std::endl;
		return;
	}

	if (!yVelocityFile.is_open())
	{
		std::cout << "Failed to open file: " << yVelocityFilePath << std::endl;
		return;
	}

	if (!pressureFile.is_open())
	{
		std::cout << "Failed to open file: " << pressureFilePath << std::endl;
		return;
	}

	for (int y = ny - 1; y >= 0; y--)
	{
		for (int x = 0; x < nx; ++x)
		{
			xVelocityFile << u(x, y) << ",";
			yVelocityFile << v(x, y) << ",";
			pressureFile << p(x, y) << ",";
		}
		xVelocityFile << std::endl;
		yVelocityFile << std::endl;
		pressureFile << std::endl;
	}

	xVelocityFile.close();
	yVelocityFile.close();
	pressureFile.close();
}


bool Equation_NavierStokes::isConverged()
{
	res_u = m_solver.getResidual(u, Ao, Ae, Aw, An, As, Sp_x);
	res_v = m_solver.getResidual(v, Ao, Ae, Aw, An, As, Sp_y);
	res_p = m_solver.getResidual(p_c, Ap_o, Ap_e, Ap_w, Ap_n, Ap_s, Sp);

	if (res_u > tolerance || res_v > tolerance || res_p > tolerance)
	{
		return false;
	}

	return true;
}


void Equation_NavierStokes::updateFaceVelocities()
{
	// Using PWIM Equations

	// Get constants
	double dx = p_mesh->getCells()[0][0].dx;
	double dy = p_mesh->getCells()[0][0].dy;
	double ao_o = 0.0;
	double ao_e = 0.0;
	double ao_w = 0.0;
	double ao_n = 0.0;
	double ao_s = 0.0;

	// 1st LAYER
	// BOTTOM LEFT CORNER (1st LAYER)
	int x = 0;
	int y = 0;

	ao_o = 1 / Ao(x, y);
	ao_e = 1 / Ao(x + 1, y);
	ao_w = 0.0;
	ao_n = 1 / Ao(x, y + 1);
	ao_s = 0.0;

	vel_face(x, y).east = 0.5 * (u(x, y) + u(x + 1, y)) +
		0.25 * dy * ao_o * (p(x + 1, y) - p(x, y)) +
		0.25 * dy * ao_e * (p(x + 2, y) - p(x, y)) -
		0.5 * dy * (ao_o + ao_e) * (p(x + 1, y) - p(x, y));

	vel_face(x, y).west = 0.0;

	vel_face(x, y).north = 0.5 * (v(x, y) + v(x, y + 1)) +
		0.25 * dx * ao_o * (p(x, y + 1) - p(x, y)) +
		0.25 * dx * ao_n * (p(x, y + 2) - p(x, y)) -
		0.5 * dx * (ao_o + ao_n) * (p(x, y + 1) - p(x, y));

	vel_face(x, y).south = 0.0;

	// x = 1, y = 0
	x = 1;
	y = 0;

	ao_o = 1 / Ao(x, y);
	ao_e = 1 / Ao(x + 1, y);
	ao_w = 1 / Ao(x - 1, y);
	ao_n = 1 / Ao(x, y + 1);
	ao_s = 0.0;

	vel_face(x, y).east = 0.5 * (u(x, y) + u(x + 1, y)) +
		0.25 * dy * ao_o * (p(x + 1, y) - p(x, y)) +
		0.25 * dy * ao_e * (p(x + 2, y) - p(x, y)) -
		0.5 * dy * (ao_o + ao_e) * (p(x + 1, y) - p(x, y));

	vel_face(x, y).west = 0.5 * (u(x, y) + u(x - 1, y)) +
		0.25 * dy * ao_o * (p(x + 1, y) - p(x - 1, y)) +
		0.25 * dy * ao_w * (p(x, y) - p(x - 1, y)) -
		0.5 * dy * (ao_o + ao_w) * (p(x, y) - p(x - 1, y));

	vel_face(x, y).north = 0.5 * (v(x, y) + v(x, y + 1)) +
		0.25 * dx * ao_o * (p(x, y + 1) - p(x, y)) +
		0.25 * dx * ao_n * (p(x, y + 2) - p(x, y)) -
		0.5 * dx * (ao_o + ao_n) * (p(x, y + 1) - p(x, y));

	vel_face(x, y).south = 0.0;


	// x = 0, y = 1
	x = 0;
	y = 1;

	ao_o = 1 / Ao(x, y);
	ao_e = 1 / Ao(x + 1, y);
	ao_w = 0.0;
	ao_n = 1 / Ao(x, y + 1);
	ao_s = 1 / Ao(x, y - 1);

	vel_face(x, y).east = 0.5 * (u(x, y) + u(x + 1, y)) +
		0.25 * dy * ao_o * (p(x + 1, y) - p(x, y)) +
		0.25 * dy * ao_e * (p(x + 2, y) - p(x, y)) -
		0.5 * dy * (ao_o + ao_e) * (p(x + 1, y) - p(x, y));

	vel_face(x, y).west = 0.0;

	vel_face(x, y).north = 0.5 * (v(x, y) + v(x, y + 1)) +
		0.25 * dx * ao_o * (p(x, y + 1) - p(x, y)) +
		0.25 * dx * ao_n * (p(x, y + 2) - p(x, y)) -
		0.5 * dx * (ao_o + ao_n) * (p(x, y + 1) - p(x, y));

	vel_face(x, y).south = 0.5 * (v(x, y) + v(x, y - 1)) +
		0.25 * dx * ao_o * (p(x, y + 1) - p(x, y - 1)) +
		0.25 * dx * ao_s * (p(x, y) - p(x, y - 1)) -
		0.5 * dx * (ao_o + ao_s) * (p(x, y) - p(x, y - 1));

	// BOTTOM RIGHT CORNER (1st LAYER)
	x = nx - 1;
	y = 0;

	ao_o = 1 / Ao(x, y);
	ao_e = 0.0;
	ao_w = 1 / Ao(x - 1, y);
	ao_n = 1 / Ao(x, y + 1);
	ao_s = 0.0;

	vel_face(x, y).east = 0.0;

	vel_face(x, y).west = 0.5 * (u(x, y) + u(x - 1, y)) +
		0.25 * dy * ao_o * (p(x, y) - p(x - 1, y)) +
		0.25 * dy * ao_w * (p(x, y) - p(x - 2, y)) -
		0.5 * dy * (ao_o + ao_w) * (p(x, y) - p(x - 1, y));

	vel_face(x, y).north = 0.5 * (v(x, y) + v(x, y + 1)) +
		0.25 * dx * ao_o * (p(x, y + 1) - p(x, y)) +
		0.25 * dx * ao_n * (p(x, y + 2) - p(x, y)) -
		0.5 * dx * (ao_o + ao_n) * (p(x, y + 1) - p(x, y));

	vel_face(x, y).south = 0.0;

	// x = nx - 2, y = 0
	x = nx - 2;
	y = 0;

	ao_o = 1 / Ao(x, y);
	ao_e = 1 / Ao(x + 1, y);
	ao_w = 1 / Ao(x - 1, y);
	ao_n = 1 / Ao(x, y + 1);
	ao_s = 0.0;

	vel_face(x, y).east = 0.5 * (u(x, y) + u(x + 1, y)) +
		0.25 * dy * ao_o * (p(x + 1, y) - p(x - 1, y)) +
		0.25 * dy * ao_e * (p(x + 1, y) - p(x, y)) -
		0.5 * dy * (ao_o + ao_e) * (p(x + 1, y) - p(x, y));

	vel_face(x, y).west = 0.5 * (u(x, y) + u(x - 1, y)) +
		0.25 * dy * ao_o * (p(x, y) - p(x - 1, y)) +
		0.25 * dy * ao_w * (p(x, y) - p(x - 2, y)) -
		0.5 * dy * (ao_o + ao_w) * (p(x, y) - p(x - 1, y));

	vel_face(x, y).north = 0.5 * (v(x, y) + v(x, y + 1)) +
		0.25 * dx * ao_o * (p(x, y + 1) - p(x, y)) +
		0.25 * dx * ao_n * (p(x, y + 2) - p(x, y)) -
		0.5 * dx * (ao_o + ao_n) * (p(x, y + 1) - p(x, y));

	vel_face(x, y).south = 0.0;

	// x = nx - 1, y = 1
	x = nx - 1;
	y = 1;

	ao_o = 1 / Ao(x, y);
	ao_e = 0.0;
	ao_w = 1 / Ao(x - 1, y);
	ao_n = 1 / Ao(x, y + 1);
	ao_s = 1 / Ao(x, y - 1);

	vel_face(x, y).east = 0.0;

	vel_face(x, y).west = 0.5 * (u(x, y) + u(x - 1, y)) +
		0.25 * dy * ao_o * (p(x, y) - p(x - 1, y)) +
		0.25 * dy * ao_w * (p(x, y) - p(x - 2, y)) -
		0.5 * dy * (ao_o + ao_w) * (p(x, y) - p(x - 1, y));

	vel_face(x, y).north = 0.5 * (v(x, y) + v(x, y + 1)) +
		0.25 * dx * ao_o * (p(x, y + 1) - p(x, y)) +
		0.25 * dx * ao_n * (p(x, y + 2) - p(x, y)) -
		0.5 * dx * (ao_o + ao_n) * (p(x, y + 1) - p(x, y));

	vel_face(x, y).south = 0.5 * (v(x, y) + v(x, y - 1)) +
		0.25 * dx * ao_o * (p(x, y + 1) - p(x, y - 1)) +
		0.25 * dx * ao_s * (p(x, y) - p(x, y - 1)) -
		0.5 * dx * (ao_o + ao_s) * (p(x, y) - p(x, y - 1));

	// TOP LEFT CORNER (1st LAYER)
	x = 0;
	y = ny - 1;

	ao_o = 1 / Ao(x, y);
	ao_e = 1 / Ao(x + 1, y);
	ao_w = 0.0;
	ao_n = 0.0;
	ao_s = 1 / Ao(x, y - 1);

	vel_face(x, y).east = 0.5 * (u(x, y) + u(x + 1, y)) +
		0.25 * dy * ao_o * (p(x + 1, y) - p(x, y)) +
		0.25 * dy * ao_e * (p(x + 2, y) - p(x, y)) -
		0.5 * dy * (ao_o + ao_e) * (p(x + 1, y) - p(x, y));

	vel_face(x, y).west = 0.0;

	vel_face(x, y).north = 0.0;

	vel_face(x, y).south = 0.5 * (v(x, y) + v(x, y - 1)) +
		0.25 * dx * ao_o * (p(x, y) - p(x, y - 1)) +
		0.25 * dx * ao_s * (p(x, y) - p(x, y - 2)) -
		0.5 * dx * (ao_o + ao_s) * (p(x, y) - p(x, y - 1));

	// x = 0, y = ny - 2
	x = 0;
	y = ny - 2;

	ao_o = 1 / Ao(x, y);
	ao_e = 1 / Ao(x + 1, y);
	ao_w = 0.0;
	ao_n = 1 / Ao(x, y + 1);
	ao_s = 1 / Ao(x, y - 1);

	vel_face(x, y).east = 0.5 * (u(x, y) + u(x + 1, y)) +
		0.25 * dy * ao_o * (p(x + 1, y) - p(x, y)) +
		0.25 * dy * ao_e * (p(x + 2, y) - p(x, y)) -
		0.5 * dy * (ao_o + ao_e) * (p(x + 1, y) - p(x, y));

	vel_face(x, y).west = 0.0;

	vel_face(x, y).north = 0.5 * (v(x, y) + v(x, y + 1)) +
		0.25 * dx * ao_o * (p(x, y + 1) - p(x, y - 1)) +
		0.25 * dx * ao_n * (p(x, y + 1) - p(x, y)) -
		0.5 * dx * (ao_o + ao_n) * (p(x, y + 1) - p(x, y));

	vel_face(x, y).south = 0.5 * (v(x, y) + v(x, y - 1)) +
		0.25 * dx * ao_o * (p(x, y) - p(x, y - 1)) +
		0.25 * dx * ao_s * (p(x, y) - p(x, y - 2)) -
		0.5 * dx * (ao_o + ao_s) * (p(x, y) - p(x, y - 1));

	// x = 1, y = ny - 1
	x = 1;
	y = ny - 1;

	ao_o = 1 / Ao(x, y);
	ao_e = 1 / Ao(x + 1, y);
	ao_w = 1 / Ao(x - 1, y);
	ao_n = 0.0;
	ao_s = 1 / Ao(x, y - 1);

	vel_face(x, y).east = 0.5 * (u(x, y) + u(x + 1, y)) +
		0.25 * dy * ao_o * (p(x + 1, y) - p(x, y)) +
		0.25 * dy * ao_e * (p(x + 2, y) - p(x, y)) -
		0.5 * dy * (ao_o + ao_e) * (p(x + 1, y) - p(x, y));

	vel_face(x, y).west = 0.5 * (u(x, y) + u(x - 1, y)) +
		0.25 * dy * ao_o * (p(x + 1, y) - p(x - 1, y)) +
		0.25 * dy * ao_w * (p(x, y) - p(x - 1, y)) -
		0.5 * dy * (ao_o + ao_w) * (p(x, y) - p(x - 1, y));

	vel_face(x, y).north = 0.0;

	vel_face(x, y).south = 0.5 * (v(x, y) + v(x, y - 1)) +
		0.25 * dx * ao_o * (p(x, y) - p(x, y - 1)) +
		0.25 * dx * ao_s * (p(x, y) - p(x, y - 2)) -
		0.5 * dx * (ao_o + ao_s) * (p(x, y) - p(x, y - 1));

	// TOP RIGHT CORNER (1st LAYER)
	x = nx - 1;
	y = ny - 1;

	ao_o = 1 / Ao(x, y);
	ao_e = 0.0;
	ao_w = 1 / Ao(x - 1, y);
	ao_n = 0.0;
	ao_s = 1 / Ao(x, y - 1);

	vel_face(x, y).east = 0.0;

	vel_face(x, y).west = 0.5 * (u(x, y) + u(x - 1, y)) +
		0.25 * dy * ao_o * (p(x, y) - p(x - 1, y)) +
		0.25 * dy * ao_w * (p(x, y) - p(x - 2, y)) -
		0.5 * dy * (ao_o + ao_w) * (p(x, y) - p(x - 1, y));

	vel_face(x, y).north = 0.0;

	vel_face(x, y).south = 0.5 * (v(x, y) + v(x, y - 1)) +
		0.25 * dx * ao_o * (p(x, y) - p(x, y - 1)) +
		0.25 * dx * ao_s * (p(x, y) - p(x, y - 2)) -
		0.5 * dx * (ao_o + ao_s) * (p(x, y) - p(x, y - 1));

	// x = nx - 2, y = ny - 1
	x = nx - 2;
	y = ny - 1;

	ao_o = 1 / Ao(x, y);
	ao_e = 1 / Ao(x + 1, y);
	ao_w = 1 / Ao(x - 1, y);
	ao_n = 0.0;
	ao_s = 1 / Ao(x, y - 1);

	vel_face(x, y).east = 0.5 * (u(x, y) + u(x + 1, y)) +
		0.25 * dy * ao_o * (p(x + 1, y) - p(x - 1, y)) +
		0.25 * dy * ao_e * (p(x + 1, y) - p(x, y)) -
		0.5 * dy * (ao_o + ao_e) * (p(x + 1, y) - p(x, y));

	vel_face(x, y).west = 0.5 * (u(x, y) + u(x - 1, y)) +
		0.25 * dy * ao_o * (p(x, y) - p(x - 1, y)) +
		0.25 * dy * ao_w * (p(x, y) - p(x - 2, y)) -
		0.5 * dy * (ao_o + ao_w) * (p(x, y) - p(x - 1, y));

	vel_face(x, y).north = 0.0;

	vel_face(x, y).south = 0.5 * (v(x, y) + v(x, y - 1)) +
		0.25 * dx * ao_o * (p(x, y) - p(x, y - 1)) +
		0.25 * dx * ao_s * (p(x, y) - p(x, y - 2)) -
		0.5 * dx * (ao_o + ao_s) * (p(x, y) - p(x, y - 1));

	// x = nx - 1, y = ny - 2
	x = nx - 1;
	y = ny - 2;

	ao_o = 1 / Ao(x, y);
	ao_e = 0.0;
	ao_w = 1 / Ao(x - 1, y);
	ao_n = 1 / Ao(x, y + 1);
	ao_s = 1 / Ao(x, y - 1);

	vel_face(x, y).east = 0.0;

	vel_face(x, y).west = 0.5 * (u(x, y) + u(x - 1, y)) +
		0.25 * dy * ao_o * (p(x, y) - p(x - 1, y)) +
		0.25 * dy * ao_w * (p(x, y) - p(x - 2, y)) -
		0.5 * dy * (ao_o + ao_w) * (p(x, y) - p(x - 1, y));

	vel_face(x, y).north = 0.5 * (v(x, y) + v(x, y + 1)) +
		0.25 * dx * ao_o * (p(x, y + 1) - p(x, y - 1)) +
		0.25 * dx * ao_n * (p(x, y + 1) - p(x, y)) -
		0.5 * dx * (ao_o + ao_n) * (p(x, y + 1) - p(x, y));

	vel_face(x, y).south = 0.5 * (v(x, y) + v(x, y - 1)) +
		0.25 * dx * ao_o * (p(x, y) - p(x, y - 1)) +
		0.25 * dx * ao_s * (p(x, y) - p(x, y - 2)) -
		0.5 * dx * (ao_o + ao_s) * (p(x, y) - p(x, y - 1));

	// BOTTOM FACE (1st LAYER)
	y = 0;
	for (x = 2; x < nx - 2; x++)
	{
		ao_o = 1 / Ao(x, y);
		ao_e = 1 / Ao(x + 1, y);
		ao_w = 1 / Ao(x - 1, y);
		ao_n = 1 / Ao(x, y + 1);
		ao_s = 0.0;

		vel_face(x, y).east = 0.5 * (u(x, y) + u(x + 1, y)) +
			0.25 * dy * ao_o * (p(x + 1, y) - p(x - 1, y)) +
			0.25 * dy * ao_e * (p(x + 2, y) - p(x, y)) -
			0.5 * dy * (ao_o + ao_e) * (p(x + 1, y) - p(x, y));

		vel_face(x, y).west = 0.5 * (u(x, y) + u(x - 1, y)) +
			0.25 * dy * ao_o * (p(x + 1, y) - p(x - 1, y)) +
			0.25 * dy * ao_w * (p(x, y) - p(x - 2, y)) -
			0.5 * dy * (ao_o + ao_w) * (p(x, y) - p(x - 1, y));

		vel_face(x, y).north = 0.5 * (v(x, y) + v(x, y + 1)) +
			0.25 * dx * ao_o * (p(x, y + 1) - p(x, y)) +
			0.25 * dx * ao_n * (p(x, y + 2) - p(x, y)) -
			0.5 * dx * (ao_o + ao_n) * (p(x, y + 1) - p(x, y));

		vel_face(x, y).south = 0.0;
	}

	// TOP FACE (1st LAYER)
	y = ny - 1;
	for (x = 2; x < nx - 2; x++)
	{
		ao_o = 1 / Ao(x, y);
		ao_e = 1 / Ao(x + 1, y);
		ao_w = 1 / Ao(x - 1, y);
		ao_n = 0.0;
		ao_s = 1 / Ao(x, y - 1);

		vel_face(x, y).east = 0.5 * (u(x, y) + u(x + 1, y)) +
			0.25 * dy * ao_o * (p(x + 1, y) - p(x - 1, y)) +
			0.25 * dy * ao_e * (p(x + 2, y) - p(x, y)) -
			0.5 * dy * (ao_o + ao_e) * (p(x + 1, y) - p(x, y));

		vel_face(x, y).west = 0.5 * (u(x, y) + u(x - 1, y)) +
			0.25 * dy * ao_o * (p(x + 1, y) - p(x - 1, y)) +
			0.25 * dy * ao_w * (p(x, y) - p(x - 2, y)) -
			0.5 * dy * (ao_o + ao_w) * (p(x, y) - p(x - 1, y));

		vel_face(x, y).north = 0.0;

		vel_face(x, y).south = 0.5 * (v(x, y) + v(x, y - 1)) +
			0.25 * dx * ao_o * (p(x, y) - p(x, y - 1)) +
			0.25 * dx * ao_s * (p(x, y) - p(x, y - 2)) -
			0.5 * dx * (ao_o + ao_s) * (p(x, y) - p(x, y - 1));
	}

	// LEFT FACE (1st LAYER)
	x = 0;
	for (y = 2; y < ny - 2; y++)
	{
		ao_o = 1 / Ao(x, y);
		ao_e = 1 / Ao(x + 1, y);
		ao_w = 0.0;
		ao_n = 1 / Ao(x, y + 1);
		ao_s = 1 / Ao(x, y - 1);

		vel_face(x, y).east = 0.5 * (u(x, y) + u(x + 1, y)) +
			0.25 * dy * ao_o * (p(x + 1, y) - p(x, y)) +
			0.25 * dy * ao_e * (p(x + 2, y) - p(x, y)) -
			0.5 * dy * (ao_o + ao_e) * (p(x + 1, y) - p(x, y));
		
		vel_face(x, y).west = 0.0;

		vel_face(x, y).north = 0.5 * (v(x, y) + v(x, y + 1)) +
			0.25 * dx * ao_o * (p(x, y + 1) - p(x, y - 1)) +
			0.25 * dx * ao_n * (p(x, y + 2) - p(x, y)) -
			0.5 * dx * (ao_o + ao_n) * (p(x, y + 1) - p(x, y));

		vel_face(x, y).south = 0.5 * (v(x, y) + v(x, y - 1)) +
			0.25 * dx * ao_o * (p(x, y + 1) - p(x, y - 1)) +
			0.25 * dx * ao_s * (p(x, y) - p(x, y - 2)) -
			0.5 * dx * (ao_o + ao_s) * (p(x, y) - p(x, y - 1));
	}

	// RIGHT FACE (1st LAYER)
	x = nx - 1;
	for (y = 2; y < ny - 2; y++)
	{
		ao_o = 1 / Ao(x, y);
		ao_e = 0.0;
		ao_w = 1 / Ao(x - 1, y);
		ao_n = 1 / Ao(x, y + 1);
		ao_s = 1 / Ao(x, y - 1);

		vel_face(x, y).east = 0.0;

		vel_face(x, y).west = 0.5 * (u(x, y) + u(x - 1, y)) +
			0.25 * dy * ao_o * (p(x, y) - p(x - 1, y)) +
			0.25 * dy * ao_w * (p(x, y) - p(x - 2, y)) -
			0.5 * dy * (ao_o + ao_w) * (p(x, y) - p(x - 1, y));

		vel_face(x, y).north = 0.5 * (v(x, y) + v(x, y + 1)) +
			0.25 * dx * ao_o * (p(x, y + 1) - p(x, y - 1)) +
			0.25 * dx * ao_n * (p(x, y + 2) - p(x, y)) -
			0.5 * dx * (ao_o + ao_n) * (p(x, y + 1) - p(x, y));

		vel_face(x, y).south = 0.5 * (v(x, y) + v(x, y - 1)) +
			0.25 * dx * ao_o * (p(x, y + 1) - p(x, y - 1)) +
			0.25 * dx * ao_s * (p(x, y) - p(x, y - 2)) -
			0.5 * dx * (ao_o + ao_s) * (p(x, y) - p(x, y - 1));
	}

	// 2nd LAYER
	// BOTTOM LEFT CORNER (2nd LAYER)
	x = 1;
	y = 1;

	ao_o = 1 / Ao(x, y);
	ao_e = 1 / Ao(x + 1, y);
	ao_w = 1 / Ao(x - 1, y);
	ao_n = 1 / Ao(x, y + 1);
	ao_s = 1 / Ao(x, y - 1);

	vel_face(x, y).east = 0.5 * (u(x, y) + u(x + 1, y)) +
		0.25 * dy * ao_o * (p(x + 1, y) - p(x - 1, y)) +
		0.25 * dy * ao_e * (p(x + 2, y) - p(x, y)) -
		0.5 * dy * (ao_o + ao_e) * (p(x + 1, y) - p(x, y));

	vel_face(x, y).west = vel_face(x - 1, y).east;
		//0.5 * (u(x, y) + u(x - 1, y)) +
		//0.25 * dy * ao_o * (p(x + 1, y) - p(x - 1, y)) +
		//0.25 * dy * ao_w * (p(x, y) - p(x - 1, y)) -
		//0.5 * dy * (ao_o + ao_w) * (p(x, y) - p(x - 1, y));

	vel_face(x, y).north = 0.5 * (v(x, y) + v(x, y + 1)) +
		0.25 * dx * ao_o * (p(x, y + 1) - p(x, y - 1)) +
		0.25 * dx * ao_n * (p(x, y + 2) - p(x, y)) -
		0.5 * dx * (ao_o + ao_n) * (p(x, y + 1) - p(x, y));

	vel_face(x, y).south = vel_face(x, y - 1).north;
		//0.5 * (v(x, y) + v(x, y - 1)) +
		//0.25 * dx * ao_o * (p(x, y + 1) - p(x, y - 1)) +
		//0.25 * dx * ao_s * (p(x, y) - p(x, y - 1)) -
		//0.5 * dx * (ao_o + ao_s) * (p(x, y) - p(x, y - 1));

	// BOTTOM RIGHT CORNER (2nd LAYER)
	x = nx - 2;
	y = 1;

	ao_o = 1 / Ao(x, y);
	ao_e = 1 / Ao(x + 1, y);
	ao_w = 1 / Ao(x - 1, y);
	ao_n = 1 / Ao(x, y + 1);
	ao_s = 1 / Ao(x, y - 1);

	vel_face(x, y).east = vel_face(x + 1, y).west;
		//0.5 * (u(x, y) + u(x + 1, y)) +
		//0.25 * dy * ao_o * (p(x + 1, y) - p(x - 1, y)) +
		//0.25 * dy * ao_e * (p(x + 1, y) - p(x, y)) -
		//0.5 * dy * (ao_o + ao_e) * (p(x + 1, y) - p(x, y));

	vel_face(x, y).west = 0.5 * (u(x, y) + u(x - 1, y)) +
		0.25 * dy * ao_o * (p(x + 1, y) - p(x - 1, y)) +
		0.25 * dy * ao_w * (p(x, y) - p(x - 2, y)) -
		0.5 * dy * (ao_o + ao_w) * (p(x, y) - p(x - 1, y));

	vel_face(x, y).north = 0.5 * (v(x, y) + v(x, y + 1)) +
		0.25 * dx * ao_o * (p(x, y + 1) - p(x, y - 1)) +
		0.25 * dx * ao_n * (p(x, y + 2) - p(x, y)) -
		0.5 * dx * (ao_o + ao_n) * (p(x, y + 1) - p(x, y));

	vel_face(x, y).south = vel_face(x, y - 1).north;
		//0.5 * (v(x, y) + v(x, y - 1)) +
		//0.25 * dx * ao_o * (p(x, y + 1) - p(x, y - 1)) +
		//0.25 * dx * ao_s * (p(x, y) - p(x, y - 1)) -
		//0.5 * dx * (ao_o + ao_s) * (p(x, y) - p(x, y - 1));

	// TOP LEFT CORNER (2nd LAYER)
	x = 1;
	y = ny - 2;

	ao_o = 1 / Ao(x, y);
	ao_e = 1 / Ao(x + 1, y);
	ao_w = 1 / Ao(x - 1, y);
	ao_n = 1 / Ao(x, y + 1);
	ao_s = 1 / Ao(x, y - 1);

	vel_face(x, y).east = 0.5 * (u(x, y) + u(x + 1, y)) +
		0.25 * dy * ao_o * (p(x + 1, y) - p(x - 1, y)) +
		0.25 * dy * ao_e * (p(x + 2, y) - p(x, y)) -
		0.5 * dy * (ao_o + ao_e) * (p(x + 1, y) - p(x, y));

	vel_face(x, y).west = vel_face(x - 1, y).east;
		//0.5 * (u(x, y) + u(x - 1, y)) +
		//0.25 * dy * ao_o * (p(x + 1, y) - p(x - 1, y)) +
		//0.25 * dy * ao_w * (p(x, y) - p(x - 1, y)) -
		//0.5 * dy * (ao_o + ao_w) * (p(x, y) - p(x - 1, y));

	vel_face(x, y).north = vel_face(x, y + 1).south;
		//0.5 * (v(x, y) + v(x, y + 1)) +
		//0.25 * dx * ao_o * (p(x, y + 1) - p(x, y - 1)) +
		//0.25 * dx * ao_n * (p(x, y + 1) - p(x, y)) -
		//0.5 * dx * (ao_o + ao_n) * (p(x, y + 1) - p(x, y));

	vel_face(x, y).south = 0.5 * (v(x, y) + v(x, y - 1)) +
		0.25 * dx * ao_o * (p(x, y + 1) - p(x, y - 1)) +
		0.25 * dx * ao_s * (p(x, y) - p(x, y - 2)) -
		0.5 * dx * (ao_o + ao_s) * (p(x, y) - p(x, y - 1));

	// TOP RIGHT CORNER (2nd LAYER)
	x = nx - 2;
	y = ny - 2;

	ao_o = 1 / Ao(x, y);
	ao_e = 1 / Ao(x + 1, y);
	ao_w = 1 / Ao(x - 1, y);
	ao_n = 1 / Ao(x, y + 1);
	ao_s = 1 / Ao(x, y - 1);

	vel_face(x, y).east = vel_face(x + 1, y).west;
		//0.5 * (u(x, y) + u(x + 1, y)) +
		//0.25 * dy * ao_o * (p(x + 1, y) - p(x - 1, y)) +
		//0.25 * dy * ao_e * (p(x + 1, y) - p(x, y)) -
		//0.5 * dy * (ao_o + ao_e) * (p(x + 1, y) - p(x, y));

	vel_face(x, y).west = 0.5 * (u(x, y) + u(x - 1, y)) +
		0.25 * dy * ao_o * (p(x + 1, y) - p(x - 1, y)) +
		0.25 * dy * ao_w * (p(x, y) - p(x - 2, y)) -
		0.5 * dy * (ao_o + ao_w) * (p(x, y) - p(x - 1, y));

	vel_face(x, y).north = vel_face(x, y + 1).south;
		//0.5 * (v(x, y) + v(x, y + 1)) +
		//0.25 * dx * ao_o * (p(x, y + 1) - p(x, y - 1)) +
		//0.25 * dx * ao_n * (p(x, y + 1) - p(x, y)) -
		//0.5 * dx * (ao_o + ao_n) * (p(x, y + 1) - p(x, y));

	vel_face(x, y).south = 0.5 * (v(x, y) + v(x, y - 1)) +
		0.25 * dx * ao_o * (p(x, y + 1) - p(x, y - 1)) +
		0.25 * dx * ao_s * (p(x, y) - p(x, y - 2)) -
		0.5 * dx * (ao_o + ao_s) * (p(x, y) - p(x, y - 1));

	// BOTTOM FACE (2nd LAYER)
	y = 1;
	for (x = 2; x < nx - 2; x++)
	{
		ao_o = 1 / Ao(x, y);
		ao_e = 1 / Ao(x + 1, y);
		ao_w = 1 / Ao(x - 1, y);
		ao_n = 1 / Ao(x, y + 1);
		ao_s = 1 / Ao(x, y - 1);

		vel_face(x, y).east = 0.5 * (u(x, y) + u(x + 1, y)) +
			0.25 * dy * ao_o * (p(x + 1, y) - p(x - 1, y)) +
			0.25 * dy * ao_e * (p(x + 2, y) - p(x, y)) -
			0.5 * dy * (ao_o + ao_e) * (p(x + 1, y) - p(x, y));

		vel_face(x, y).west = 0.5 * (u(x, y) + u(x - 1, y)) +
			0.25 * dy * ao_o * (p(x + 1, y) - p(x - 1, y)) +
			0.25 * dy * ao_w * (p(x, y) - p(x - 2, y)) -
			0.5 * dy * (ao_o + ao_w) * (p(x, y) - p(x - 1, y));

		vel_face(x, y).north = 0.5 * (v(x, y) + v(x, y + 1)) +
			0.25 * dx * ao_o * (p(x, y + 1) - p(x, y - 1)) +
			0.25 * dx * ao_n * (p(x, y + 2) - p(x, y)) -
			0.5 * dx * (ao_o + ao_n) * (p(x, y + 1) - p(x, y));

		vel_face(x, y).south = vel_face(x, y - 1).north;
			//0.5 * (v(x, y) + v(x, y - 1)) +
			//0.25 * dx * ao_o * (p(x, y + 1) - p(x, y - 1)) +
			//0.25 * dx * ao_s * (p(x, y) - p(x, y - 1)) -
			//0.5 * dx * (ao_o + ao_s) * (p(x, y) - p(x, y - 1));
	}

	// TOP FACE (2nd LAYER)
	y = ny - 2;
	for (x = 2; x < nx - 2; x++)
	{
		ao_o = 1 / Ao(x, y);
		ao_e = 1 / Ao(x + 1, y);
		ao_w = 1 / Ao(x - 1, y);
		ao_n = 1 / Ao(x, y + 1);
		ao_s = 1 / Ao(x, y - 1);

		vel_face(x, y).east = 0.5 * (u(x, y) + u(x + 1, y)) +
			0.25 * dy * ao_o * (p(x + 1, y) - p(x - 1, y)) +
			0.25 * dy * ao_e * (p(x + 2, y) - p(x, y)) -
			0.5 * dy * (ao_o + ao_e) * (p(x + 1, y) - p(x, y));

		vel_face(x, y).west = 0.5 * (u(x, y) + u(x - 1, y)) +
			0.25 * dy * ao_o * (p(x + 1, y) - p(x - 1, y)) +
			0.25 * dy * ao_w * (p(x, y) - p(x - 2, y)) -
			0.5 * dy * (ao_o + ao_w) * (p(x, y) - p(x - 1, y));

		vel_face(x, y).north = vel_face(x, y + 1).south;
			//0.5 * (v(x, y) + v(x, y + 1)) +
			//0.25 * dx * ao_o * (p(x, y + 1) - p(x, y - 1)) +
			//0.25 * dx * ao_n * (p(x, y + 1) - p(x, y)) -
			//0.5 * dx * (ao_o + ao_n) * (p(x, y + 1) - p(x, y));

		vel_face(x, y).south = 0.5 * (v(x, y) + v(x, y - 1)) +
			0.25 * dx * ao_o * (p(x, y + 1) - p(x, y - 1)) +
			0.25 * dx * ao_s * (p(x, y) - p(x, y - 2)) -
			0.5 * dx * (ao_o + ao_s) * (p(x, y) - p(x, y - 1));
	}

	// LEFT FACE (2nd LAYER)
	x = 1;
	for (y = 2; y < ny - 2; y++)
	{
		ao_o = 1 / Ao(x, y);
		ao_e = 1 / Ao(x + 1, y);
		ao_w = 1 / Ao(x - 1, y);
		ao_n = 1 / Ao(x, y + 1);
		ao_s = 1 / Ao(x, y - 1);

		vel_face(x, y).east = 0.5 * (u(x, y) + u(x + 1, y)) +
			0.25 * dy * ao_o * (p(x + 1, y) - p(x - 1, y)) +
			0.25 * dy * ao_e * (p(x + 2, y) - p(x, y)) -
			0.5 * dy * (ao_o + ao_e) * (p(x + 1, y) - p(x, y));

		vel_face(x, y).west = vel_face(x - 1, y).east;
			//0.5 * (u(x, y) + u(x - 1, y)) +
			//0.25 * dy * ao_o * (p(x + 1, y) - p(x - 1, y)) +
			//0.25 * dy * ao_w * (p(x, y) - p(x - 1, y)) -
			//0.5 * dy * (ao_o + ao_w) * (p(x, y) - p(x - 1, y));

		vel_face(x, y).north = 0.5 * (v(x, y) + v(x, y + 1)) +
			0.25 * dx * ao_o * (p(x, y + 1) - p(x, y - 1)) +
			0.25 * dx * ao_n * (p(x, y + 2) - p(x, y)) -
			0.5 * dx * (ao_o + ao_n) * (p(x, y + 1) - p(x, y));

		vel_face(x, y).south = 0.5 * (v(x, y) + v(x, y - 1)) +
			0.25 * dx * ao_o * (p(x, y + 1) - p(x, y - 1)) +
			0.25 * dx * ao_s * (p(x, y) - p(x, y - 2)) -
			0.5 * dx * (ao_o + ao_s) * (p(x, y) - p(x, y - 1));

	}

	// RIGHT FACE (2nd LAYER)
	x = nx - 2;
	for (y = 2; y < ny - 2; y++)
	{
		ao_o = 1 / Ao(x, y);
		ao_e = 1 / Ao(x + 1, y);
		ao_w = 1 / Ao(x - 1, y);
		ao_n = 1 / Ao(x, y + 1);
		ao_s = 1 / Ao(x, y - 1);

		vel_face(x, y).east = vel_face(x + 1, y).west;
			//0.5 * (u(x, y) + u(x + 1, y)) +
			//0.25 * dy * ao_o * (p(x + 1, y) - p(x - 1, y)) +
			//0.25 * dy * ao_e * (p(x + 1, y) - p(x, y)) -
			//0.5 * dy * (ao_o + ao_e) * (p(x + 1, y) - p(x, y));

		vel_face(x, y).west = 0.5 * (u(x, y) + u(x - 1, y)) +
			0.25 * dy * ao_o * (p(x + 1, y) - p(x - 1, y)) +
			0.25 * dy * ao_w * (p(x, y) - p(x - 2, y)) -
			0.5 * dy * (ao_o + ao_w) * (p(x, y) - p(x - 1, y));

		vel_face(x, y).north = 0.5 * (v(x, y) + v(x, y + 1)) +
			0.25 * dx * ao_o * (p(x, y + 1) - p(x, y - 1)) +
			0.25 * dx * ao_n * (p(x, y + 2) - p(x, y)) -
			0.5 * dx * (ao_o + ao_n) * (p(x, y + 1) - p(x, y));

		vel_face(x, y).south = 0.5 * (v(x, y) + v(x, y - 1)) +
			0.25 * dx * ao_o * (p(x, y + 1) - p(x, y - 1)) +
			0.25 * dx * ao_s * (p(x, y) - p(x, y - 2)) -
			0.5 * dx * (ao_o + ao_s) * (p(x, y) - p(x, y - 1));
	}

	// INTERIOR CELLS
	for (x = 2; x < nx - 2; ++x)
	{
		for (y = 2; y < ny - 2; ++y)
		{
			ao_o = 1 / Ao(x, y);
			ao_e = 1 / Ao(x + 1, y);
			ao_w = 1 / Ao(x - 1, y);
			ao_n = 1 / Ao(x, y + 1);
			ao_s = 1 / Ao(x, y - 1);

			vel_face(x, y).east = 0.5 * (u(x, y) + u(x + 1, y)) +
				0.25 * dy * ao_o * (p(x + 1, y) - p(x - 1, y)) +
				0.25 * dy * ao_e * (p(x + 2, y) - p(x, y)) -
				0.5 * dy * (ao_o + ao_e) * (p(x + 1, y) - p(x, y));

			vel_face(x, y).west = 0.5 * (u(x, y) + u(x - 1, y)) +
				0.25 * dy * ao_o * (p(x + 1, y) - p(x - 1, y)) +
				0.25 * dy * ao_w * (p(x, y) - p(x - 2, y)) -
				0.5 * dy * (ao_o + ao_w) * (p(x, y) - p(x - 1, y));

			vel_face(x, y).north = 0.5 * (v(x, y) + v(x, y + 1)) +
				0.25 * dx * ao_o * (p(x, y + 1) - p(x, y - 1)) +
				0.25 * dx * ao_n * (p(x, y + 2) - p(x, y)) -
				0.5 * dx * (ao_o + ao_n) * (p(x, y + 1) - p(x, y));

			vel_face(x, y).south = 0.5 * (v(x, y) + v(x, y - 1)) +
				0.25 * dx * ao_o * (p(x, y + 1) - p(x, y - 1)) +
				0.25 * dx * ao_s * (p(x, y) - p(x, y - 2)) -
				0.5 * dx * (ao_o + ao_s) * (p(x, y) - p(x, y - 1));
		}
	}

	// END OF MESH
}


void Equation_NavierStokes::updatePressureLinks()
{
	// Get constants
	double dx = p_mesh->getCells()[0][0].dx;
	double dy = p_mesh->getCells()[0][0].dy;
	double dyy = dy * dy;
	double dxx = dx * dx;
	double ao_o = 0.0;
	double ao_e = 0.0;
	double ao_w = 0.0;
	double ao_n = 0.0;
	double ao_s = 0.0;

	// BOTTOM LEFT CORNER
	int x = 0;
	int y = 0;

	ao_o = 1 / Ao(x, y);
	ao_e = 1 / Ao(x + 1, y);
	ao_w = 0.0;
	ao_n = 1 / Ao(x, y + 1);
	ao_s = 0.0;

	double rho_e = rho_face(x, y).east;
	double rho_w = rho_face(x, y).west;
	double rho_n = rho_face(x, y).north;
	double rho_s = rho_face(x, y).south;

	double u_e = vel_face(x, y).east;
	double u_w = vel_face(x, y).west;
	double v_n = vel_face(x, y).north;
	double v_s = vel_face(x, y).south;

	Ap_e(x, y) = -0.5 * rho_e * dyy * (ao_o + ao_e);
	Ap_w(x, y) = 0.0;
	Ap_n(x, y) = -0.5 * rho_n * dxx * (ao_o + ao_n);
	Ap_s(x, y) = 0.0;
	Ap_o(x, y) = -Ap_e(x, y) - Ap_n(x, y);
	Sp(x, y) = -(dy * (rho_e * u_e) + dx * (rho_n * v_n));

	// BOTTOM RIGHT CORNER
	x = nx - 1;
	y = 0;

	ao_o = 1 / Ao(x, y);
	ao_e = 0.0;
	ao_w = 1 / Ao(x - 1, y);
	ao_n = 1 / Ao(x, y + 1);
	ao_s = 0.0;

	rho_e = rho_face(x, y).east;
	rho_w = rho_face(x, y).west;
	rho_n = rho_face(x, y).north;
	rho_s = rho_face(x, y).south;

	u_e = vel_face(x, y).east;
	u_w = vel_face(x, y).west;
	v_n = vel_face(x, y).north;
	v_s = vel_face(x, y).south;

	Ap_e(x, y) = 0.0;
	Ap_w(x, y) = -0.5 * rho_w * dyy * (ao_o + ao_w);
	Ap_n(x, y) = -0.5 * rho_n * dxx * (ao_o + ao_n);
	Ap_s(x, y) = 0.0;
	Ap_o(x, y) = -Ap_w(x, y) - Ap_n(x, y);
	Sp(x, y) = -(dy * (- rho_w * u_w) + dx * (rho_n * v_n));

	// TOP LEFT
	x = 0;
	y = ny - 1;

	ao_o = 1 / Ao(x, y);
	ao_e = 1 / Ao(x + 1, y);
	ao_w = 0.0;
	ao_n = 0.0;
	ao_s = 1 / Ao(x, y - 1);

	rho_e = rho_face(x, y).east;
	rho_w = rho_face(x, y).west;
	rho_n = rho_face(x, y).north;
	rho_s = rho_face(x, y).south;

	u_e = vel_face(x, y).east;
	u_w = vel_face(x, y).west;
	v_n = vel_face(x, y).north;
	v_s = vel_face(x, y).south;

	Ap_e(x, y) = -0.5 * rho_e * dyy * (ao_o + ao_e);
	Ap_w(x, y) = 0.0;
	Ap_n(x, y) = 0.0;
	Ap_s(x, y) = -0.5 * rho_s * dxx * (ao_o + ao_s);
	Ap_o(x, y) = -Ap_e(x, y) - Ap_s(x, y);
	Sp(x, y) = -(dy * (rho_e * u_e) + dx * (- rho_s * v_s));

	// TOP RIGHT
	x = nx - 1;
	y = ny - 1;

	ao_o = 1 / Ao(x, y);
	ao_e = 0.0;
	ao_w = 1 / Ao(x - 1, y);
	ao_n = 0.0;
	ao_s = 1 / Ao(x, y - 1);

	rho_e = rho_face(x, y).east;
	rho_w = rho_face(x, y).west;
	rho_n = rho_face(x, y).north;
	rho_s = rho_face(x, y).south;

	u_e = vel_face(x, y).east;
	u_w = vel_face(x, y).west;
	v_n = vel_face(x, y).north;
	v_s = vel_face(x, y).south;

	Ap_e(x, y) = 0.0;
	Ap_w(x, y) = -0.5 * rho_w * dyy * (ao_o + ao_w);
	Ap_n(x, y) = 0.0;
	Ap_s(x, y) = -0.5 * rho_s * dxx * (ao_o + ao_s);
	Ap_o(x, y) = -Ap_w(x, y) - Ap_s(x, y);
	Sp(x, y) = -(dy * (- rho_w * u_w) + dx * (- rho_s * v_s));

	// BOTTOM FACE
	y = 0;
	for (x = 1; x < nx - 1; ++x)
	{
		ao_o = 1 / Ao(x, y);
		ao_e = 1 / Ao(x + 1, y);
		ao_w = 1 / Ao(x - 1, y);
		ao_n = 1 / Ao(x, y + 1);
		ao_s = 0.0;

		rho_e = rho_face(x, y).east;
		rho_w = rho_face(x, y).west;
		rho_n = rho_face(x, y).north;
		rho_s = rho_face(x, y).south;

		u_e = vel_face(x, y).east;
		u_w = vel_face(x, y).west;
		v_n = vel_face(x, y).north;
		v_s = vel_face(x, y).south;

		Ap_e(x, y) = -0.5 * rho_e * dyy * (ao_o + ao_e);
		Ap_w(x, y) = -0.5 * rho_w * dyy * (ao_o + ao_w);
		Ap_n(x, y) = -0.5 * rho_n * dxx * (ao_o + ao_n);
		Ap_s(x, y) = 0.0;
		Ap_o(x, y) = -Ap_e(x, y) - Ap_w(x, y) - Ap_n(x, y);
		Sp(x, y) = -(dy * (rho_e * u_e - rho_w * u_w) + dx * (rho_n * v_n));
	}

	// TOP FACE
	y = ny - 1;
	for (x = 1; x < nx - 1; ++x)
	{
		ao_o = 1 / Ao(x, y);
		ao_e = 1 / Ao(x + 1, y);
		ao_w = 1 / Ao(x - 1, y);
		ao_n = 0.0;
		ao_s = 1 / Ao(x, y - 1);

		rho_e = rho_face(x, y).east;
		rho_w = rho_face(x, y).west;
		rho_n = rho_face(x, y).north;
		rho_s = rho_face(x, y).south;

		u_e = vel_face(x, y).east;
		u_w = vel_face(x, y).west;
		v_n = vel_face(x, y).north;
		v_s = vel_face(x, y).south;

		Ap_e(x, y) = -0.5 * rho_e * dyy * (ao_o + ao_e);
		Ap_w(x, y) = -0.5 * rho_w * dyy * (ao_o + ao_w);
		Ap_n(x, y) = 0.0;
		Ap_s(x, y) = -0.5 * rho_s * dxx * (ao_o + ao_s);
		Ap_o(x, y) = -Ap_e(x, y) - Ap_w(x, y) - Ap_n(x, y) - Ap_s(x, y);
		Sp(x, y) = -(dy * (rho_e * u_e - rho_w * u_w) + dx * (- rho_s * v_s));
	}

	// LEFT FACE
	x = 0;
	for (y = 1; y < ny - 1; ++y)
	{
		ao_o = 1 / Ao(x, y);
		ao_e = 1 / Ao(x + 1, y);
		ao_w = 0.0;
		ao_n = 1 / Ao(x, y + 1);
		ao_s = 1 / Ao(x, y - 1);

		rho_e = rho_face(x, y).east;
		rho_w = rho_face(x, y).west;
		rho_n = rho_face(x, y).north;
		rho_s = rho_face(x, y).south;

		u_e = vel_face(x, y).east;
		u_w = vel_face(x, y).west;
		v_n = vel_face(x, y).north;
		v_s = vel_face(x, y).south;

		Ap_e(x, y) = -0.5 * rho_e * dyy * (ao_o + ao_e);
		Ap_w(x, y) = 0.0;
		Ap_n(x, y) = -0.5 * rho_n * dxx * (ao_o + ao_n);
		Ap_s(x, y) = -0.5 * rho_s * dxx * (ao_o + ao_s);
		Ap_o(x, y) = -Ap_e(x, y) - Ap_n(x, y) - Ap_s(x, y);
		Sp(x, y) = -(dy * (rho_e * u_e) + dx * (rho_n * v_n - rho_s * v_s));
	}

	// RIGHT FACE
	x = nx - 1;
	for (y = 1; y < ny - 1; ++y)
	{
		ao_o = 1 / Ao(x, y);
		ao_e = 0.0;
		ao_w = 1 / Ao(x - 1, y);
		ao_n = 1 / Ao(x, y + 1);
		ao_s = 1 / Ao(x, y - 1);

		rho_e = rho_face(x, y).east;
		rho_w = rho_face(x, y).west;
		rho_n = rho_face(x, y).north;
		rho_s = rho_face(x, y).south;

		u_e = vel_face(x, y).east;
		u_w = vel_face(x, y).west;
		v_n = vel_face(x, y).north;
		v_s = vel_face(x, y).south;

		Ap_e(x, y) = 0.0;
		Ap_w(x, y) = -0.5 * rho_w * dyy * (ao_o + ao_w);
		Ap_n(x, y) = -0.5 * rho_n * dxx * (ao_o + ao_n);
		Ap_s(x, y) = -0.5 * rho_s * dxx * (ao_o + ao_s);
		Ap_o(x, y) = - Ap_w(x, y) - Ap_n(x, y) - Ap_s(x, y);
		Sp(x, y) = -(dy * (- rho_w * u_w) + dx * (rho_n * v_n - rho_s * v_s));
	}

	// INTERIOR CELLS
	for (x = 1; x < nx - 1; ++x)
	{
		for (y = 1; y < ny - 1; ++y)
		{
			ao_o = 1 / Ao(x, y);
			ao_e = 1 / Ao(x + 1, y);
			ao_w = 1 / Ao(x - 1, y);
			ao_n = 1 / Ao(x, y + 1);
			ao_s = 1 / Ao(x, y - 1);

			rho_e = rho_face(x, y).east;
			rho_w = rho_face(x, y).west;
			rho_n = rho_face(x, y).north;
			rho_s = rho_face(x, y).south;

			u_e = vel_face(x, y).east;
			u_w = vel_face(x, y).west;
			v_n = vel_face(x, y).north;
			v_s = vel_face(x, y).south;

			Ap_e(x, y) = -0.5 * rho_e * dyy * (ao_o + ao_e);
			Ap_w(x, y) = -0.5 * rho_w * dyy * (ao_o + ao_w);
			Ap_n(x, y) = -0.5 * rho_n * dxx * (ao_o + ao_n);
			Ap_s(x, y) = -0.5 * rho_s * dxx * (ao_o + ao_s);
			Ap_o(x, y) = -Ap_e(x, y) - Ap_w(x, y) - Ap_n(x, y) - Ap_s(x, y);
			Sp(x, y) = -(dy * (rho_e * u_e - rho_w * u_w) + dx * (rho_n * v_n - rho_s * v_s));
		}
	}
}


void Equation_NavierStokes::updateVelocityCorrection()
{
	// Get constants
	double dx = p_mesh->getCells()[0][0].dx;
	double dy = p_mesh->getCells()[0][0].dy;
	double ao_o = 0.0;
	double ao_e = 0.0;
	double ao_w = 0.0;
	double ao_n = 0.0;
	double ao_s = 0.0;

	// BOTTOM LEFT CORNER
	int x = 0;
	int y = 0;

	ao_o = 1 / Ao(x, y);
	ao_e = 1 / Ao(x + 1, y);
	ao_w = 0.0;
	ao_n = 1 / Ao(x, y + 1);
	ao_s = 0.0;

	//u_c(x, y) = 0.25 * dy * ao_o * (p_c(x, y) - p_c(x + 1, y));
	u_c(x, y) = 0.5 * dy * ao_o * (p_c(x, y) - p_c(x + 1, y));
	//u_c(x, y) = 0.5 * dy * ao_o * (-p_c(x + 1, y));
	//v_c(x, y) = 0.25 * dx * ao_o * (p_c(x, y) - p_c(x, y + 1));
	v_c(x, y) = 0.5 * dx * ao_o * (p_c(x, y) - p_c(x, y + 1));
	//v_c(x, y) = 0.5 * dx * ao_o * (-p_c(x, y + 1));

	vel_c_face(x, y).east = 0.5 * dy * (ao_o + ao_e) * (p_c(x, y) - p_c(x + 1, y));
	vel_c_face(x, y).west = 0.0;
	vel_c_face(x, y).north = 0.5 * dx * (ao_o + ao_n) * (p_c(x, y) - p_c(x, y + 1));
	vel_c_face(x, y).south = 0.0;

	// BOTTOM RIGHT CORNER
	x = nx - 1;
	y = 0;

	ao_o = 1 / Ao(x, y);
	ao_e = 0.0;
	ao_w = 1 / Ao(x - 1, y);
	ao_n = 1 / Ao(x, y + 1);
	ao_s = 0.0;

	//u_c(x, y) = 0.25 * dy * ao_o * (p_c(x - 1, y) - p_c(x, y));
	u_c(x, y) = 0.5 * dy * ao_o * (p_c(x - 1, y) - p_c(x, y));
	//u_c(x, y) = 0.5 * dy * ao_o * (p_c(x - 1, y));
	//v_c(x, y) = 0.25 * dx * ao_o * (p_c(x, y) - p_c(x, y + 1));
	v_c(x, y) = 0.5 * dx * ao_o * (p_c(x, y) - p_c(x, y + 1));
	//v_c(x, y) = 0.5 * dx * ao_o * (-p_c(x, y + 1));

	vel_c_face(x, y).east = 0.0;
	vel_c_face(x, y).west = 0.5 * dy * (ao_o + ao_w) * (p_c(x - 1, y) - p_c(x, y));
	vel_c_face(x, y).north = 0.5 * dx * (ao_o + ao_n) * (p_c(x, y) - p_c(x, y + 1));
	vel_c_face(x, y).south = 0.0;

	// TOP LEFT CORNER
	x = 0;
	y = ny - 1;

	ao_o = 1 / Ao(x, y);
	ao_e = 1 / Ao(x + 1, y);
	ao_w = 0.0;
	ao_n = 0.0;
	ao_s = 1 / Ao(x, y - 1);

	//u_c(x, y) = 0.25 * dy * ao_o * (p_c(x, y) - p_c(x + 1, y));
	u_c(x, y) = 0.5 * dy * ao_o * (p_c(x, y) - p_c(x + 1, y));
	//u_c(x, y) = 0.5 * dy * ao_o * (-p_c(x + 1, y));
	//v_c(x, y) = 0.25 * dx * ao_o * (p_c(x, y - 1) - p_c(x, y));
	v_c(x, y) = 0.5 * dx * ao_o * (p_c(x, y - 1) - p_c(x, y));
	//v_c(x, y) = 0.5 * dx * ao_o * (p_c(x, y - 1));

	vel_c_face(x, y).east = 0.5 * dy * (ao_o + ao_e) * (p_c(x, y) - p_c(x + 1, y));
	vel_c_face(x, y).west = 0.0;
	vel_c_face(x, y).north = 0.0;
	vel_c_face(x, y).south = 0.5 * dx * (ao_o + ao_s) * (p_c(x, y - 1) - p_c(x, y));

	// TOP RIGHT CORNER
	x = nx - 1;
	y = ny - 1;

	ao_o = 1 / Ao(x, y);
	ao_e = 0.0;
	ao_w = 1 / Ao(x - 1, y);
	ao_n = 0.0;
	ao_s = 1 / Ao(x, y - 1);

	//u_c(x, y) = 0.25 * dy * ao_o * (p_c(x - 1, y) - p_c(x, y));
	u_c(x, y) = 0.5 * dy * ao_o * (p_c(x - 1, y) - p_c(x, y));
	//u_c(x, y) = 0.5 * dy * ao_o * (p_c(x - 1, y));
	//v_c(x, y) = 0.25 * dx * ao_o * (p_c(x, y - 1) - p_c(x, y));
	v_c(x, y) = 0.5 * dx * ao_o * (p_c(x, y - 1) - p_c(x, y));
	//v_c(x, y) = 0.5 * dx * ao_o * (p_c(x, y - 1));

	vel_c_face(x, y).east = 0.0;
	vel_c_face(x, y).west = 0.5 * dy * (ao_o + ao_w) * (p_c(x - 1, y) - p_c(x, y));
	vel_c_face(x, y).north = 0.0;
	vel_c_face(x, y).south = 0.5 * dx * (ao_o + ao_s) * (p_c(x, y - 1) - p_c(x, y));

	// BOTTOM FACE
	y = 0;
	for (x = 1; x < nx - 1; ++x)
	{
		ao_o = 1 / Ao(x, y);
		ao_e = 1 / Ao(x + 1, y);
		ao_w = 1 / Ao(x - 1, y);
		ao_n = 1 / Ao(x, y + 1);
		ao_s = 0.0;

		u_c(x, y) = 0.5 * dy * ao_o * (p_c(x - 1, y) - p_c(x + 1, y));
		//v_c(x, y) = 0.25 * dx * ao_o * (p_c(x, y) - p_c(x, y + 1));
		v_c(x, y) = 0.5 * dx * ao_o * (p_c(x, y) - p_c(x, y + 1));
		//v_c(x, y) = 0.5 * dx * ao_o * (-p_c(x, y + 1));

		vel_c_face(x, y).east = 0.5 * dy * (ao_o + ao_e) * (p_c(x, y) - p_c(x + 1, y));
		vel_c_face(x, y).west = 0.5 * dy * (ao_o + ao_w) * (p_c(x - 1, y) - p_c(x, y));
		vel_c_face(x, y).north = 0.5 * dx * (ao_o + ao_n) * (p_c(x, y) - p_c(x, y + 1));
		vel_c_face(x, y).south = 0.0;
	}

	// TOP FACE
	y = ny - 1;
	for (x = 1; x < nx - 1; ++x)
	{
		ao_o = 1 / Ao(x, y);
		ao_e = 1 / Ao(x + 1, y);
		ao_w = 1 / Ao(x - 1, y);
		ao_n = 0.0;
		ao_s = 1 / Ao(x, y - 1);

		u_c(x, y) = 0.5 * dy * ao_o * (p_c(x - 1, y) - p_c(x + 1, y));
		//v_c(x, y) = 0.25 * dx * ao_o * (p_c(x, y - 1) - p_c(x, y));
		v_c(x, y) = 0.5 * dx * ao_o * (p_c(x, y - 1) - p_c(x, y));
		//v_c(x, y) = 0.5 * dx * ao_o * (p_c(x, y - 1));

		vel_c_face(x, y).east = 0.5 * dy * (ao_o + ao_e) * (p_c(x, y) - p_c(x + 1, y));
		vel_c_face(x, y).west = 0.5 * dy * (ao_o + ao_w) * (p_c(x - 1, y) - p_c(x, y));
		vel_c_face(x, y).north = 0.0;
		vel_c_face(x, y).south = 0.5 * dx * (ao_o + ao_s) * (p_c(x, y - 1) - p_c(x, y));
	}

	// LEFT FACE
	x = 0;
	for (y = 1; y < ny - 1; ++y)
	{
		ao_o = 1 / Ao(x, y);
		ao_e = 1 / Ao(x + 1, y);
		ao_w = 0.0;
		ao_n = 1 / Ao(x, y + 1);
		ao_s = 1 / Ao(x, y - 1);

		//u_c(x, y) = 0.25 * dy * ao_o * (p_c(x, y) - p_c(x + 1, y));
		u_c(x, y) = 0.5 * dy * ao_o * (p_c(x, y) - p_c(x + 1, y));
		//u_c(x, y) = 0.5 * dy * ao_o * (-p_c(x + 1, y));
		v_c(x, y) = 0.5 * dx * ao_o * (p_c(x, y - 1) - p_c(x, y + 1));

		vel_c_face(x, y).east = 0.5 * dy * (ao_o + ao_e) * (p_c(x, y) - p_c(x + 1, y));
		vel_c_face(x, y).west = 0.0;
		vel_c_face(x, y).north = 0.5 * dx * (ao_o + ao_n) * (p_c(x, y) - p_c(x, y + 1));
		vel_c_face(x, y).south = 0.5 * dx * (ao_o + ao_s) * (p_c(x, y - 1) - p_c(x, y));
	}

	// RIGHT FACE
	x = nx - 1;
	for (y = 1; y < ny - 1; ++y)
	{
		ao_o = 1 / Ao(x, y);
		ao_e = 0.0;
		ao_w = 1 / Ao(x - 1, y);
		ao_n = 1 / Ao(x, y + 1);
		ao_s = 1 / Ao(x, y - 1);

		//u_c(x, y) = 0.25 * dy * ao_o * (p_c(x - 1, y) - p_c(x, y));
		u_c(x, y) = 0.5 * dy * ao_o * (p_c(x - 1, y) - p_c(x, y));
		//u_c(x, y) = 0.5 * dy * ao_o * (p_c(x - 1, y));
		v_c(x, y) = 0.5 * dx * ao_o * (p_c(x, y - 1) - p_c(x, y + 1));

		vel_c_face(x, y).east = 0.0;
		vel_c_face(x, y).west = 0.5 * dy * (ao_o + ao_w) * (p_c(x - 1, y) - p_c(x, y));
		vel_c_face(x, y).north = 0.5 * dx * (ao_o + ao_n) * (p_c(x, y) - p_c(x, y + 1));
		vel_c_face(x, y).south = 0.5 * dx * (ao_o + ao_s) * (p_c(x, y - 1) - p_c(x, y));
	}

	// INTERIOR CELLS
	for (x = 1; x < nx - 1; ++x)
	{
		for (y = 1; y < ny - 1; ++y)
		{
			ao_o = 1 / Ao(x, y);
			ao_e = 1 / Ao(x + 1, y);
			ao_w = 1 / Ao(x - 1, y);
			ao_n = 1 / Ao(x, y + 1);
			ao_s = 1 / Ao(x, y - 1);

			u_c(x, y) = 0.5 * dy * ao_o * (p_c(x - 1, y) - p_c(x + 1, y));
			v_c(x, y) = 0.5 * dx * ao_o * (p_c(x, y - 1) - p_c(x, y + 1));
			
			vel_c_face(x, y).east = 0.5 * dy * (ao_o + ao_e) * (p_c(x, y) - p_c(x + 1, y));
			vel_c_face(x, y).west = 0.5 * dy * (ao_o + ao_w) * (p_c(x - 1, y) - p_c(x, y));
			vel_c_face(x, y).north = 0.5 * dx * (ao_o + ao_n) * (p_c(x, y) - p_c(x, y + 1));
			vel_c_face(x, y).south = 0.5 * dx * (ao_o + ao_s) * (p_c(x, y - 1) - p_c(x, y));
		}
	}

	// END OF MESH
}


void Equation_NavierStokes::correctVelocityAndPressure()
{
	// ALL CELLS
	for (int x = 0; x < nx; ++x)
	{
		for (int y = 0; y < ny; ++y)
		{
			p(x, y) = p(x, y) + p_relax * p_c(x, y);
			u(x, y) = u(x, y) + vel_relax * u_c(x, y);
			v(x, y) = v(x, y) + vel_relax * v_c(x, y);
			vel_face(x, y).east = vel_face(x, y).east + vel_relax * vel_c_face(x, y).east;
			vel_face(x, y).west = vel_face(x, y).west + vel_relax * vel_c_face(x, y).west;
			vel_face(x, y).north = vel_face(x, y).north + vel_relax * vel_c_face(x, y).north;
			vel_face(x, y).south = vel_face(x, y).south + vel_relax * vel_c_face(x, y).south;
		}
	}
}


void Equation_NavierStokes::updateLinkCoefficient()
{
	// Get constants
	double dx = p_mesh->getCells()[0][0].dx;
	double dy = p_mesh->getCells()[0][0].dy;
	double dyx = dy / dx;
	double dxy = dx / dy;

	double mflux = 0.0;

	// BOTTOM LEFT CORNER
	int x = 0;
	int y = 0;
	
	mflux_face(x, y).east = rho_face(x, y).east * vel_face(x, y).east;
	Ae(x, y) = -0.5 * dy * (abs(mflux_face(x, y).east) - mflux_face(x, y).east) -
		vis_face(x, y).east * dyx;

	mflux_face(x, y).west = 0.0;
	Aw(x, y) = 0.0;

	mflux_face(x, y).north = rho_face(x, y).north * vel_face(x, y).north;
	An(x, y) = -0.5 * dx * (abs(mflux_face(x, y).north) - mflux_face(x, y).north) -
		vis_face(x, y).north * dxy;

	mflux_face(x, y).south = 0.0;
	As(x, y) = 0.0;

	Ao(x, y) = 0.5 * dy * (abs(mflux_face(x, y).east) + mflux_face(x, y).east) +
		0.5 * dx * (abs(mflux_face(x, y).north) + mflux_face(x, y).north) +
		vis_face(x, y).east * dyx +
		2 * vis_face(x, y).west * dxy +
		vis_face(x, y).north * dxy + 
		2 * vis_face(x, y).south * dyx;

	Sp_x(x, y) = 0.5 * dy * (p(x, y) - p(x + 1, y));
	Sp_y(x, y) = 0.5 * dx * (p(x, y) - p(x, y + 1));

	// BOTTOM RIGHT CORNER
	x = nx - 1;
	y = 0;

	mflux_face(x, y).east = 0.0;
	Ae(x, y) = 0.0;

	mflux_face(x, y).west = rho_face(x, y).west * vel_face(x, y).west;
	Aw(x, y) = -0.5 * dy * (abs(mflux_face(x, y).west) + mflux_face(x, y).west) -
		vis_face(x, y).west * dyx;

	mflux_face(x, y).north = rho_face(x, y).north * vel_face(x, y).north;
	An(x, y) = -0.5 * dx * (abs(mflux_face(x, y).north) - mflux_face(x, y).north) -
		vis_face(x, y).north * dxy;

	mflux_face(x, y).south = 0.0;
	As(x, y) = 0.0;

	Ao(x, y) = 0.5 * dy * (abs(mflux_face(x, y).west) - mflux_face(x, y).west) +
		0.5 * dx * (abs(mflux_face(x, y).north) + mflux_face(x, y).north) +
		2 * vis_face(x, y).east * dxy + 
		vis_face(x, y).west * dyx +
		vis_face(x, y).north * dxy +
		2 * vis_face(x, y).south * dyx;
	
	Sp_x(x, y) = 0.5 * dy * (p(x - 1, y) - p(x, y));
	Sp_y(x, y) = 0.5 * dx * (p(x, y) - p(x, y + 1));

	// TOP LEFT CORNER
	x = 0;
	y = ny - 1;

	mflux_face(x, y).east = rho_face(x, y).east * vel_face(x, y).east;
	Ae(x, y) = -0.5 * dy * (abs(mflux_face(x, y).east) - mflux_face(x, y).east) -
		vis_face(x, y).east * dyx;

	mflux_face(x, y).west = 0.0;
	Aw(x, y) = 0.0;

	mflux_face(x, y).north = 0.0;
	An(x, y) = 0.0;

	mflux_face(x, y).south = rho_face(x, y).south * vel_face(x, y).south;
	As(x, y) = -0.5 * dx * (abs(mflux_face(x, y).south) + mflux_face(x, y).south) -
		vis_face(x, y).south * dxy;

	Ao(x, y) = 0.5 * dy * (abs(mflux_face(x, y).east) + mflux_face(x, y).east) +
		0.5 * dx * (abs(mflux_face(x, y).south) - mflux_face(x, y).south) +
		vis_face(x, y).east * dyx +
		2 * vis_face(x, y).west * dyx +
		2 * vis_face(x, y).north * dxy +
		vis_face(x, y).south * dxy;

	Sp_x(x, y) = 0.5 * dy * (p(x, y) - p(x + 1, y)) + 2 * u_lid * vis_face(x, y).north * dxy;
	Sp_y(x, y) = 0.5 * dx * (p(x, y - 1) - p(x, y));

	// TOP RIGHT CORNER
	x = nx - 1;
	y = ny - 1;

	mflux_face(x, y).east = 0.0;
	Ae(x, y) = 0.0;

	mflux_face(x, y).west = rho_face(x, y).west * vel_face(x, y).west;
	Aw(x, y) = -0.5 * dy * (abs(mflux_face(x, y).west) + mflux_face(x, y).west) -
		vis_face(x, y).west * dyx;

	mflux_face(x, y).north = 0.0;
	An(x, y) = 0.0;

	mflux_face(x, y).south = rho_face(x, y).south * vel_face(x, y).south;
	As(x, y) = -0.5 * dx * (abs(mflux_face(x, y).south) + mflux_face(x, y).south) -
		vis_face(x, y).south * dxy;

	Ao(x, y) = 0.5 * dy * (abs(mflux_face(x, y).west) - mflux_face(x, y).west) +
		0.5 * dx * (abs(mflux_face(x, y).south) - mflux_face(x, y).south) +
		2 * vis_face(x, y).east * dyx +
		vis_face(x, y).west * dyx +
		vis_face(x, y).south * dxy +  
		2 * vis_face(x, y).north * dxy;

	Sp_x(x, y) = 0.5 * dy * (p(x - 1, y) - p(x, y)) + 2 * u_lid * vis_face(x, y).north * dxy;
	Sp_y(x, y) = 0.5 * dx * (p(x, y - 1) - p(x, y));

	// BOTTOM FACE
	y = 0;
	for (x = 1; x < nx - 1; x++)
	{
		mflux_face(x, y).east = rho_face(x, y).east * vel_face(x, y).east;
		Ae(x, y) = -0.5 * dy * (abs(mflux_face(x, y).east) - mflux_face(x, y).east) -
			vis_face(x, y).east * dyx;

		mflux_face(x, y).west = rho_face(x, y).west * vel_face(x, y).west;
		Aw(x, y) = -0.5 * dy * (abs(mflux_face(x, y).west) + mflux_face(x, y).west) -
			vis_face(x, y).west * dyx;

		mflux_face(x, y).north = rho_face(x, y).north * vel_face(x, y).north;
		An(x, y) = -0.5 * dx * (abs(mflux_face(x, y).north) - mflux_face(x, y).north) -
			vis_face(x, y).north * dxy;

		mflux_face(x, y).south = 0.0;
		As(x, y) = 0.0;

		Ao(x, y) = 0.5 * dy * (abs(mflux_face(x, y).east) + mflux_face(x, y).east) +
			0.5 * dy * (abs(mflux_face(x, y).west) - mflux_face(x, y).west) +
			0.5 * dx * (abs(mflux_face(x, y).north) + mflux_face(x, y).north) +
			vis_face(x, y).east * dyx +
			vis_face(x, y).west * dyx +
			vis_face(x, y).north * dxy +
			2 * vis_face(x, y).south * dxy;

		Sp_x(x, y) = 0.5 * dy * (p(x - 1, y) - p(x + 1, y));
		Sp_y(x, y) = 0.5 * dx * (p(x, y) - p(x, y + 1));
	}

	// TOP FACE
	y = ny - 1;
	for (x = 1; x < nx - 1; x++)
	{
		mflux_face(x, y).east = rho_face(x, y).east * vel_face(x, y).east;
		Ae(x, y) = -0.5 * dy * (abs(mflux_face(x, y).east) - mflux_face(x, y).east) -
			vis_face(x, y).east * dyx;

		mflux_face(x, y).west = rho_face(x, y).west * vel_face(x, y).west;
		Aw(x, y) = -0.5 * dy * (abs(mflux_face(x, y).west) + mflux_face(x, y).west) -
			vis_face(x, y).west * dyx;

		mflux_face(x, y).north = 0.0;
		An(x, y) = 0.0;

		mflux_face(x, y).south = rho_face(x, y).south * vel_face(x, y).south;
		As(x, y) = -0.5 * dx * (abs(mflux_face(x, y).south) + mflux_face(x, y).south) -
			vis_face(x, y).south * dxy;

		Ao(x, y) = 0.5 * dy * (abs(mflux_face(x, y).east) + mflux_face(x, y).east) +
			0.5 * dy * (abs(mflux_face(x, y).west) - mflux_face(x, y).west) +
			0.5 * dx * (abs(mflux_face(x, y).south) - mflux_face(x, y).south) +
			vis_face(x, y).east * dyx +
			vis_face(x, y).west * dyx +
			2 * vis_face(x, y).north * dxy +
			vis_face(x, y).south * dxy;

		Sp_x(x, y) = 0.5 * dy * (p(x - 1, y) - p(x + 1, y)) + 2 * u_lid * vis_face(x, y).north * dxy;
		Sp_y(x, y) = 0.5 * dx * (p(x, y - 1) - p(x, y));
	}

	// LEFT FACE
	x = 0;
	for (y = 1; y < ny - 1; y++)
	{
		mflux_face(x, y).east = rho_face(x, y).east * vel_face(x, y).east;
		Ae(x, y) = -0.5 * dy * (abs(mflux_face(x, y).east) - mflux_face(x, y).east) -
			vis_face(x, y).east * dyx;

		mflux_face(x, y).west = 0.0;
		Aw(x, y) = 0.0;

		mflux_face(x, y).north = rho_face(x, y).north * vel_face(x, y).north;
		An(x, y) = -0.5 * dx * (abs(mflux_face(x, y).north) - mflux_face(x, y).north) -
			vis_face(x, y).north * dxy;

		mflux_face(x, y).south = rho_face(x, y).south * vel_face(x, y).south;
		As(x, y) = -0.5 * dx * (abs(mflux_face(x, y).south) + mflux_face(x, y).south) -
			vis_face(x, y).south * dxy;

		Ao(x, y) = 0.5 * dy * (abs(mflux_face(x, y).east) + mflux_face(x, y).east) +
			0.5 * dx * (abs(mflux_face(x, y).north) + mflux_face(x, y).north) +
			0.5 * dx * (abs(mflux_face(x, y).south) - mflux_face(x, y).south) +
			vis_face(x, y).east * dyx +
			2 * vis_face(x, y).west * dyx +
			vis_face(x, y).north * dxy +
			vis_face(x, y).south * dxy;

		Sp_x(x, y) = 0.5 * dy * (p(x, y) - p(x + 1, y));
		Sp_y(x, y) = 0.5 * dx * (p(x, y - 1) - p(x, y + 1));
	}

	// RIGHT FACE
	x = nx - 1;
	for (y = 1; y < ny - 1; y++)
	{
		mflux_face(x, y).east = 0.0;
		Ae(x, y) = 0.0;

		mflux_face(x, y).west = rho_face(x, y).west * vel_face(x, y).west;
		Aw(x, y) = -0.5 * dy * (abs(mflux_face(x, y).west) + mflux_face(x, y).west) -
			vis_face(x, y).west * dyx;

		mflux_face(x, y).north = rho_face(x, y).north * vel_face(x, y).north;
		An(x, y) = -0.5 * dx * (abs(mflux_face(x, y).north) - mflux_face(x, y).north) -
			vis_face(x, y).north * dxy;

		mflux_face(x, y).south = rho_face(x, y).south * vel_face(x, y).south;
		As(x, y) = -0.5 * dx * (abs(mflux_face(x, y).south) + mflux_face(x, y).south) -
			vis_face(x, y).south * dxy;

		Ao(x, y) = 0.5 * dy * (abs(mflux_face(x, y).west) - mflux_face(x, y).west) +
			0.5 * dx * (abs(mflux_face(x, y).north) + mflux_face(x, y).north) +
			0.5 * dx * (abs(mflux_face(x, y).south) - mflux_face(x, y).south) +
			2 * vis_face(x, y).east * dyx + 
			vis_face(x, y).west * dyx +
			vis_face(x, y).north * dxy +
			vis_face(x, y).south * dxy;

		Sp_x(x, y) = 0.5 * dy * (p(x - 1, y) - p(x, y));
		Sp_y(x, y) = 0.5 * dx * (p(x, y - 1) - p(x, y + 1));
	}

	// INTERIOR NODES
	for (x = 1; x < nx - 1; ++x)
	{
		for (y = 1; y < ny - 1; ++y)
		{
			mflux_face(x, y).east = rho_face(x, y).east * vel_face(x, y).east;
			Ae(x, y) = -0.5 * dy * (abs(mflux_face(x, y).east) - mflux_face(x, y).east) -
				vis_face(x, y).east * dyx;

			mflux_face(x, y).west = rho_face(x, y).west * vel_face(x, y).west;
			Aw(x, y) = -0.5 * dy * (abs(mflux_face(x, y).west) + mflux_face(x, y).west) - 
				vis_face(x, y).west * dyx;
			
			mflux_face(x, y).north = rho_face(x, y).north * vel_face(x, y).north;
			An(x, y) = -0.5 * dx * (abs(mflux_face(x, y).north) - mflux_face(x, y).north) -
				vis_face(x, y).north * dxy;

			mflux_face(x, y).south = rho_face(x, y).south * vel_face(x, y).south;
			As(x, y) = -0.5 * dx * (abs(mflux_face(x, y).south) + mflux_face(x, y).south) -
				vis_face(x, y).south * dxy;

			Ao(x, y) = 0.5 * dy * (abs(mflux_face(x, y).east) + mflux_face(x, y).east) +
				0.5 * dy * (abs(mflux_face(x, y).west) - mflux_face(x, y).west) +
				0.5 * dx * (abs(mflux_face(x, y).north) + mflux_face(x, y).north) +
				0.5 * dx * (abs(mflux_face(x, y).south) - mflux_face(x, y).south) +
				vis_face(x, y).east * dyx + 
				vis_face(x, y).west * dyx +
				vis_face(x, y).north * dxy + 
				vis_face(x, y).south * dxy;

			Sp_x(x, y) = 0.5 * dy * (p(x - 1, y) - p(x + 1, y));
			Sp_y(x, y) = 0.5 * dx * (p(x, y - 1) - p(x, y + 1));
		}
	}

	// END OF MESH
}
