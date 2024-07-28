#include "2DFluidSimulator/Domain.h"
#include <cassert>
#include <filesystem>
#include <iostream>


void Domain::initialise()
{
	assert(m_equations.size() > 0);

	for (auto& equation : m_equations)
	{
		equation->initialise(m_mesh.getNumberOfXCells(), m_mesh.getNumberOfYCells(), &m_mesh);
	}

	m_mesh.generateMesh();
}

void Domain::addEquation(std::unique_ptr<IEquation> equation)
{
	m_equations.emplace_back(std::move(equation));
}

void Domain::run()
{
	initialise();

	for (auto& equation : m_equations)
	{
		equation->update();
	}

	// Log data
	std::string resultsPath = "Results/";
	std::filesystem::create_directory(resultsPath);
	
	for (auto& equation : m_equations)
	{
		equation->logData(resultsPath);
	}
}

