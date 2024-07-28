#include "2DFluidSimulator/Equation_NavierStokes.h"

void Equation_NavierStokes::initialiseEquation(int numberOfXCells, int numberOfYCells)
{
	m_numberOfXCells = numberOfXCells;
	m_numberOfYCells = numberOfYCells;

	u_x.setConstant(numberOfXCells * numberOfYCells, 0.0);
	u_face_x.setConstant(numberOfXCells * numberOfYCells, 0.0);
	u_y.setConstant(numberOfXCells * numberOfYCells, 0.0);
	u_face_y.setConstant(numberOfXCells * numberOfYCells, 0.0);
	p.setConstant(numberOfXCells * numberOfYCells, 0.0);
	p_corr.setConstant(numberOfXCells * numberOfYCells, 0.0);
}


void Equation_NavierStokes::update()
{
	updateLinkCoefficient();

}

void Equation_NavierStokes::updateLinkCoefficient()
{
	int N = m_numberOfXCells * m_numberOfYCells;

	for (int n = 0; n < N; ++n)
	{

	}
}
