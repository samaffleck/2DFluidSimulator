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

	v.setConstant(numberOfXCells * numberOfYCells, 0);
}

void Equation_Diffusion::update()
{
	int N = m_numberOfXCells * m_numberOfYCells;

	double dx = p_mesh->getCells()[0][0].dx;
	double dy = p_mesh->getCells()[0][0].dy;
	double dx_squared_inv = 1 / (dx * dx);
	double dy_squared_inv = 1 / (dy * dy);

	// Update link coefficients
	for (int n = 0; n < N; ++n)
	{
		m_solver.Ao[n] = 2 * dx_squared_inv + 2 * dy_squared_inv;
		m_solver.Ae[n] = - dx_squared_inv; //
		m_solver.Aw[n] = - dx_squared_inv; //
		m_solver.An[n] = - dy_squared_inv; //
		m_solver.As[n] = - dy_squared_inv; //
		m_solver.S[n] = 0.0; //
	}

	// Apply boundary conditions?
	
	// Solve the equation
	m_solver.solve(v);
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
			file << v[x + y * m_numberOfXCells] << ",";
		}
		file << std::endl;
	}

	file.close();
}
