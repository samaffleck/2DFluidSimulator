#include <iostream>

#include "2DFluidSimulator/Domain.h"
#include "2DFluidSimulator/Equation_Diffusion.h"
#include "2DFluidSimulator/Equation_NavierStokes.h"

int main() 
{
	Domain domain;
	domain.getMesh().setNumberOfXCells(10);
	domain.getMesh().setNumberOfYCells(10);
	domain.getMesh().setDomainLengthX(1);
	domain.getMesh().setDomainLengthY(1);

	//auto diffusion = std::make_unique<Equation_Diffusion>();
	//domain.addEquation(std::move(diffusion));

	auto ns = std::make_unique<Equation_NavierStokes>();
	domain.addEquation(std::move(ns));

	domain.run();
	 
	return 0;
}
