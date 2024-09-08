#include "2DFluidSimulator/Equation_Diffusion.h"
#include "2DFluidSimulator/Mesh.h"
#include <fstream>
#include <sstream>
#include <filesystem>
#include <iostream>


void Equation_Diffusion::initialiseEquation(int numberOfXCells, int numberOfYCells)
{
	m_numberOfXCells = numberOfXCells;
	m_numberOfYCells = numberOfYCells;

	v.setConstant(numberOfXCells, numberOfYCells, 0);
	Ao.setConstant(numberOfXCells, numberOfYCells, 0);
	Ae.setConstant(numberOfXCells, numberOfYCells, 0);
	Aw.setConstant(numberOfXCells, numberOfYCells, 0);
	An.setConstant(numberOfXCells, numberOfYCells, 0);
	As.setConstant(numberOfXCells, numberOfYCells, 0);
	S.setConstant(numberOfXCells, numberOfYCells, 0);
}


void Equation_Diffusion::update()
{
	// Get constants
	int nx = m_numberOfXCells;
	int ny = m_numberOfYCells;
	double dx = p_mesh->getCells()[0][0].dx;
	double dy = p_mesh->getCells()[0][0].dy;
	double dxx = 1 / (dx * dx);
	double dyy = 1 / (dy * dy);

	// BOTTOM LEFT CORNER
	int x = 0;
	int y = 0;

	Ao(x, y) = 3 * dxx + 3 * dyy;
	Ae(x, y) = -dxx;
	Aw(x, y) = 0.0;
	An(x, y) = -dyy;
	As(x, y) = 0.0;
	S(x, y) = 2 * dxx * m_solver.leftBoundary + 2 * dyy * m_solver.bottomBoundary;

	// BOTTOM RIGHT CORNER
	x = nx - 1;
	y = 0;

	Ao(x, y) = 3 * dxx + 3 * dyy;
	Ae(x, y) = 0.0;
	Aw(x, y) = -dxx;
	An(x, y) = -dyy;
	As(x, y) = 0.0;
	S(x, y) = 2 * dxx * m_solver.rightBoundary + 2 * dyy * m_solver.bottomBoundary;

	// TOP LEFT CORNER
	x = 0;
	y = ny - 1;

	Ao(x, y) = 3 * dxx + 3 * dyy;
	Ae(x, y) = -dxx;
	Aw(x, y) = 0.0;
	An(x, y) = 0.0;
	As(x, y) = -dyy;
	S(x, y) = 2 * dxx * m_solver.leftBoundary + 2 * dyy * m_solver.topBoundary;

	// TOP RIGHT CORNER
	x = nx - 1;
	y = ny - 1;

	Ao(x, y) = 3 * dxx + 3 * dyy;
	Ae(x, y) = 0.0;
	Aw(x, y) = -dxx;
	An(x, y) = 0.0;
	As(x, y) = -dyy;
	S(x, y) = 2 * dxx * m_solver.rightBoundary + 2 * dyy * m_solver.topBoundary;

	// BOTTOM WALL
	y = 0;
	for (x = 1; x < nx - 1; ++x)
	{
		Ao(x, y) = 2 * dxx + 3 * dyy;
		Ae(x, y) = -dxx;
		Aw(x, y) = -dxx;
		An(x, y) = -dyy;
		As(x, y) = 0.0;
		S(x, y) = 2 * dyy * m_solver.bottomBoundary;
	}

	// TOP WALL 
	y = ny - 1;
	for (x = 1; x < nx - 1; ++x)
	{
		Ao(x, y) = 2 * dxx + 3 * dyy;
		Ae(x, y) = -dxx;
		Aw(x, y) = -dxx;
		An(x, y) = 0.0;
		As(x, y) = -dyy;
		S(x, y) = 2 * dyy * m_solver.topBoundary;
	}

	// LEFT WALL
	x = 0;
	for (y = 1; y < ny - 1; ++y)
	{
		Ao(x, y) = 3 * dxx + 2 * dyy;
		Ae(x, y) = -dxx;
		Aw(x, y) = 0.0;
		An(x, y) = -dyy;
		As(x, y) = -dyy;
		S(x, y) = 2 * dxx * m_solver.leftBoundary;
	}

	// RIGHT WALL
	x = nx - 1;
	for (y = 1; y < ny - 1; ++y)
	{
		Ao(x, y) = 3 * dxx + 2 * dyy;
		Ae(x, y) = 0.0;
		Aw(x, y) = -dxx;
		An(x, y) = -dyy;
		As(x, y) = -dyy;
		S(x, y) = 2 * dxx * m_solver.rightBoundary;
	}

	// INTERIOR CELLS	
	for (x = 1; x < nx - 1; ++x)
	{
		for (y = 1; y < ny - 1; ++y)
		{
			Ao(x, y) = 2 * dxx + 2 * dyy;
			Ae(x, y) = -dxx; 
			Aw(x, y) = -dxx; 
			An(x, y) = -dyy; 
			As(x, y) = -dyy; 
			S(x, y) = 0.0; 
		}
	}

	// Solve the equation
	m_solver.solve(v, Ao, Ae, Aw, An, As, S, normalisedResidual, 1e-6);
}


void Equation_Diffusion::logData(std::string resultsDirectory)
{
	std::filesystem::path filePath = std::filesystem::path(resultsDirectory) / "ScalarValue.csv";

	std::ofstream file(filePath);

	if (!file.is_open()) {
		std::cout << "Failed to open file: " << filePath << std::endl;
		return;
	}

	for (int y = m_numberOfYCells - 1; y >= 0; y--)
	{
		for (int x = 0; x < m_numberOfXCells; ++x)
		{
			file << v(x, y) << ",";
		}
		file << std::endl;
	}

	file.close();
}
