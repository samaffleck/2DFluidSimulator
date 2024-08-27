#include <iostream>

#include "2DFluidSimulator/Domain.h"
#include "2DFluidSimulator/Equation_Diffusion.h"

int main() 
{
	Domain domain;
	domain.getMesh().setNumberOfXCells(10);
	domain.getMesh().setNumberOfYCells(10);
	domain.getMesh().setDomainLengthX(1);
	domain.getMesh().setDomainLengthY(1);

	auto diffusion = std::make_unique<Equation_Diffusion>();
	diffusion->setLeftBoundaryCondition(100);
	diffusion->setRightBoundaryCondition(0);
	diffusion->setTopBoundaryCondition(0);
	diffusion->setBottomBoundaryCondition(100);
	domain.addEquation(std::move(diffusion));

	domain.run();

	return 0;
}
