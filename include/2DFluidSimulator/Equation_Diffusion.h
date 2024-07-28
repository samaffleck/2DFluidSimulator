#pragma once

#include "IEquation.h"


class Equation_Diffusion : public IEquation
{
public: 
	void initialiseEquation(int numberOfXCells, int numberOfYCells) override;
	void update() override;
	void logData(std::string resultsDirectory) override;

public:
	double D = 1e-5;				// Diffusion coefficient

public:
	Eigen::VectorXd& getV() { return v; }

private:
	Eigen::VectorXd v{};	// Scalar field, e.g. Temperature
	int m_numberOfXCells;
	int m_numberOfYCells;

};
