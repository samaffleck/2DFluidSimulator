#include "2DFluidSimulator/IEquation.h"
#include "2DFluidSimulator/Mesh.h"

void IEquation::initialise(int numberOfXCells, int numberOfYCells, Mesh* mesh)
{
	int numberOfCells = numberOfXCells * numberOfYCells;
	m_solver.Ao.setConstant(numberOfCells, 0.0);
	m_solver.Ae.setConstant(numberOfCells, 0.0);
	m_solver.Aw.setConstant(numberOfCells, 0.0);
	m_solver.An.setConstant(numberOfCells, 0.0);
	m_solver.As.setConstant(numberOfCells, 0.0);
	m_solver.S.setConstant(numberOfCells, 0.0);

	initialiseEquation(numberOfXCells, numberOfYCells);
	p_mesh = mesh;
	m_solver.initialise(numberOfXCells);
}
