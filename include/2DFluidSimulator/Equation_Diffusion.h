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

private:
	Eigen::MatrixXd v{};	// Scalar field, e.g. Temperature
	
	// Link coefficients
	Eigen::MatrixXd Ao{};
	Eigen::MatrixXd Ae{};
	Eigen::MatrixXd Aw{};
	Eigen::MatrixXd An{};
	Eigen::MatrixXd As{};
	Eigen::MatrixXd S{};		

	int m_numberOfXCells{};
	int m_numberOfYCells{};
	double normalisedResidual = 1;

};
