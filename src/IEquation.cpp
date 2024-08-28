#include "2DFluidSimulator/IEquation.h"
#include "2DFluidSimulator/Mesh.h"

void IEquation::initialise(int numberOfXCells, int numberOfYCells, Mesh* mesh)
{
	m_solver.initialise(numberOfXCells, numberOfYCells);
	initialiseEquation(numberOfXCells, numberOfYCells);
	p_mesh = mesh;
}

void IEquation::setLeftBoundaryCondition(double boundaryValue)
{
	m_solver.leftBoundary = boundaryValue;
}

void IEquation::setRightBoundaryCondition(double boundaryValue)
{
	m_solver.rightBoundary = boundaryValue;
}

void IEquation::setTopBoundaryCondition(double boundaryValue)
{
	m_solver.topBoundary = boundaryValue;
}

void IEquation::setBottomBoundaryCondition(double boundaryValue)
{
	m_solver.bottomBoundary = boundaryValue;
}
