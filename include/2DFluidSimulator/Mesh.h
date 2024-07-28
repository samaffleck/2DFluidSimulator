#pragma once

#include "2DFluidSimulator/Cell.h"
#include <vector>

class Mesh
{
public:
	Mesh() = default;
	~Mesh() = default;

	// API
	void generateMesh();
	std::vector<std::vector<Cell>>& getCells() { return m_cells; }

	// Setters
	void setNumberOfXCells(int numberOfCells);
	void setNumberOfYCells(int numberOfCells);
	void setDomainLengthX(float length);
	void setDomainLengthY(float length);

	// Getters
	int getNumberOfXCells() { return m_numberOfCellsX; }
	int getNumberOfYCells() { return m_numberOfCellsY; }
	float getDomainLengthX() { return m_lengthX; }
	float getDomainLengthY() { return m_lengthY; }

private:
	int m_numberOfCellsX = 2;
	int m_numberOfCellsY = 2;
	float m_lengthX = 1;
	float m_lengthY = 1;
	std::vector<std::vector<Cell>> m_cells{};

};
