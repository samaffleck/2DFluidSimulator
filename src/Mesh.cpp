#include "2DFluidSimulator/Mesh.h"
#include <cassert>


void Mesh::generateMesh()
{
	m_cells.resize(m_numberOfCellsX);
	
	float dx = m_lengthX / m_numberOfCellsX;
	float dy = m_lengthY / m_numberOfCellsY;

	for (int x = 0; x < m_numberOfCellsX; ++x)
	{
		m_cells[x].resize(m_numberOfCellsY);
		for (int y = 0; y < m_numberOfCellsY; ++y)
		{
			m_cells[x][y].dx = dx;
			m_cells[x][y].dy = dy;
			m_cells[x][y].x = x * dx + dx / 2;
			m_cells[x][y].y = y * dy + dy / 2;
		}
	}
}

void Mesh::setNumberOfXCells(int numberOfCells)
{
	assert(numberOfCells > 1);
	m_numberOfCellsX = numberOfCells;
}

void Mesh::setNumberOfYCells(int numberOfCells)
{
	assert(numberOfCells > 1);
	m_numberOfCellsY = numberOfCells;
}

void Mesh::setDomainLengthX(float length)
{
	assert(length > 0);
	m_lengthX = length;
}

void Mesh::setDomainLengthY(float length)
{
	assert(length > 0);
	m_lengthX = length;
}

