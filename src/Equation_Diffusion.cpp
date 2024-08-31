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
	int nx = m_numberOfXCells;
	int ny = m_numberOfYCells;

	double dx = p_mesh->getCells()[0][0].dx;
	double dy = p_mesh->getCells()[0][0].dy;
	double dx_squared_inv = 1 / (dx * dx);
	double dy_squared_inv = 1 / (dy * dy);

	// Update link coefficients
	for (int x = 0; x < nx; ++x)
	{
		for (int y = 0; y < ny; ++y)
		{
			Ao(x, y) = 2 * dx_squared_inv + 2 * dy_squared_inv;
			Ae(x, y) = -dx_squared_inv; //
			Aw(x, y) = -dx_squared_inv; //
			An(x, y) = -dy_squared_inv; //
			As(x, y) = -dy_squared_inv; //
			S(x, y) = 0.0; //
		}
	}
	
	// Solve the equation
	m_solver.solve(v, Ao, Ae, Aw, An, As, S);
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
