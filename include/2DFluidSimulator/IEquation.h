#pragma once

#include "2DFluidSimulator/ODESolver.h"
#include "2DFluidSimulator/Mesh.h"
#include "Eigen/Dense"


class IEquation
{
public:
	IEquation() = default;
	virtual ~IEquation() = default;

	// Virtual functions
	virtual void logData(std::string resultsDirectory) = 0;
	virtual void update() = 0;

	void initialise(int numberOfXCells, int numberOfYCells, Mesh* mesh);

protected:
	virtual void initialiseEquation(int numberOfXCells, int numberOfYCells) = 0;
	Mesh* p_mesh = nullptr; // Non owning pointer to the mesh
	ODESolver m_solver;
	
	// Boundary Conditions?
	// Initial Conditions?

};
