#include <iostream>

#include "2DFluidSimulator/Domain.h"
#include "2DFluidSimulator/Equation_Diffusion.h"
#include "2DFluidSimulator/Equation_NavierStokes.h"


bool testDiffusion()
{
	Domain domain;
	domain.getMesh().setNumberOfXCells(10);
	domain.getMesh().setNumberOfYCells(10);
	domain.getMesh().setDomainLengthX(1);
	domain.getMesh().setDomainLengthY(1);

	auto diffusion = std::make_unique<Equation_Diffusion>();
	diffusion->setLeftBoundaryCondition(100);
	diffusion->setBottomBoundaryCondition(100);
	diffusion->setRightBoundaryCondition(0);
	diffusion->setTopBoundaryCondition(0);
	domain.addEquation(std::move(diffusion));
	domain.run();

	return true;
}


bool testNavierStokes()
{
	Domain domain;
	domain.getMesh().setNumberOfXCells(10);
	domain.getMesh().setNumberOfYCells(10);
	domain.getMesh().setDomainLengthX(1);
	domain.getMesh().setDomainLengthY(1);

	auto navierStokes = std::make_unique<Equation_NavierStokes>();
	domain.addEquation(std::move(navierStokes));
	domain.run();

	return true;
}


int main() 
{
	//testDiffusion();
	testNavierStokes();

	return 0;
}
 