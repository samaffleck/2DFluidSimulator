#include "2DFluidSimulator/Equation_NavierStokes.h"


void Equation_NavierStokes::initialiseEquation(int numberOfXCells, int numberOfYCells)
{
	nx = numberOfXCells;
	ny = numberOfYCells;
	
	double viscosity = 1e-5;
	double density = 1.0;

	u.setConstant(nx, ny, 0.0);
	v.setConstant(nx, ny, 0.0);
	p.setConstant(nx, ny, 0.0);
	p_c.setConstant(nx, ny, 0.0);
	vis.setConstant(nx, ny, viscosity);
	rho.setConstant(nx, ny, density);
	Sp_x.setConstant(nx, ny, 0.0);
	Sp_y.setConstant(nx, ny, 0.0);

	vel_face.resize(nx, ny);
	rho_face.resize(nx, ny);
	vis_face.resize(nx, ny);
	mflux_face.resize(nx, ny);

	for (int x = 0; x < nx; ++x)
	{
		for (int y = 0; y < ny; ++y)
		{
			vis_face(x, y).east = viscosity;
			vis_face(x, y).west = viscosity;
			vis_face(x, y).north = viscosity;
			vis_face(x, y).south = viscosity;
		}
	}
}


void Equation_NavierStokes::update()
{
	int itt = 0;
	int maxItterations = 2e5;
	bool hasConverged = false;

	while (!hasConverged && itt < maxItterations)
	{
		// Update momentum links and sources for x and y-direction
		updateLinkCoefficient();

		// solve x-momentum
		m_solver.S = Sp_x;
		m_solver.solve(u);

		// solve y-momentum
		m_solver.S = Sp_y;
		m_solver.solve(v);

		// Calculate the face velocities using PWIM
		// updateFaceVelocities();

		// Calculate pressure links and mass imbalance
		// updatePressureLinks();

		// Solve pressure correction equation
		//m_solver.solve(p_corr);

		// Correct velocities and pressure
		// correctVelocitiesAndPressure();

		// Check for overall convergance
		// hasConverged = isConverged();

		itt++;
	}
}


void Equation_NavierStokes::updateLinkCoefficient()
{
	// Get constants
	double dx = p_mesh->getCells()[0][0].dx;
	double dy = p_mesh->getCells()[0][0].dy;

	double mflux = 0.0;

	// BOTTOM LEFT CORNER
	int x = 0;
	int y = 0;

	mflux = rho_face(x, y).east * vel_face(x, y).east;
	mflux_face(x, y).east = 0.5 * (abs(mflux) - mflux);
	m_solver.Ae(x, y) = -dy * mflux_face(x, y).east - vis_face(x, y).east * dy / dx;

	///
	///
	/// 


	// INTERIOR NODES
	for (x = 1; x < nx - 1; ++x)
	{
		for (y = 1; y < ny - 1; ++y)
		{
			mflux = rho_face(x, y).east * vel_face(x, y).east;
			mflux_face(x, y).east = 0.5 * (abs(mflux) - mflux);
			m_solver.Ae(x, y) = -dy * mflux_face(x, y).east - vis_face(x, y).east * dy / dx;

			mflux = rho_face(x, y).west * vel_face(x, y).west;
			mflux_face(x, y).west = 0.5 * (abs(mflux) - mflux);
			m_solver.Aw(x, y) = -dy * mflux_face(x, y).west - vis_face(x, y).west * dy / dx;
			
			mflux = rho_face(x, y).north * vel_face(x, y).north;
			mflux_face(x, y).north = 0.5 * (abs(mflux) - mflux);
			m_solver.An(x, y) = -dx * mflux_face(x, y).north - vis_face(x, y).north * dx / dy;

			mflux = rho_face(x, y).south * vel_face(x, y).south;
			mflux_face(x, y).south = 0.5 * (abs(mflux) - mflux);
			m_solver.As(x, y) = -dx * mflux_face(x, y).south - vis_face(x, y).south * dx / dy;

			m_solver.S(x, y) = 0.0;
		}
	}


}
