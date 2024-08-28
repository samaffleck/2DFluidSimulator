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
	double dyx = dy / dx;
	double dxy = dx / dy;

	double mflux = 0.0;

	// BOTTOM LEFT CORNER
	int x = 0;
	int y = 0;
	
	mflux = rho_face(x, y).east * vel_face(x, y).east;
	mflux_face(x, y).east = 0.5 * (abs(mflux) - mflux);
	m_solver.Ae(x, y) = -dy * mflux_face(x, y).east - vis_face(x, y).east * dyx - vis(x, y) * dyx / 3;

	m_solver.Aw(x, y) = 0.0;

	mflux = rho_face(x, y).north * vel_face(x, y).north;
	mflux_face(x, y).north = 0.5 * (abs(mflux) - mflux);
	m_solver.An(x, y) = -dx * mflux_face(x, y).north - vis_face(x, y).north * dxy - vis(x, y) * dxy / 3;
	
	m_solver.As(x, y) = 0.0;

	m_solver.Ao(x, y) = dy * mflux_face(x, y).east + dx * mflux_face(x, y).north + 
						vis_face(x, y).north * dxy + 3 * vis(x, y) * dxy + 
						vis_face(x, y).east * dyx + 3 * vis(x, y) * dyx;

	Sp_x(x, y) = 0.5 * dy * (p(x, y) - p(x + 1, y));
	Sp_y(x, y) = 0.5 * dx * (p(x, y) - p(x, y + 1));

	// BOTTOM RIGHT CORNER
	x = nx - 1;
	y = 0;

	m_solver.Ae(x, y) = 0.0;

	mflux = rho_face(x, y).west * vel_face(x, y).west;
	mflux_face(x, y).west = 0.5 * (abs(mflux) - mflux);
	m_solver.Aw(x, y) = -dy * mflux_face(x, y).west - vis_face(x, y).west * dyx + vis(x, y) * dyx / 3;

	mflux = rho_face(x, y).north * vel_face(x, y).north;
	mflux_face(x, y).north = 0.5 * (abs(mflux) - mflux);
	m_solver.An(x, y) = -dx * mflux_face(x, y).north - vis_face(x, y).north * dxy - vis(x, y) * dxy / 3;
	
	m_solver.As(x, y) = 0.0;

	m_solver.Ao(x, y) = dy * mflux_face(x, y).west + dx * mflux_face(x, y).north + 
						vis_face(x, y).north * dxy + 3 * vis(x, y) * dxy + 
						vis_face(x, y).west * dyx + 3 * vis(x, y) * dyx;

	Sp_x(x, y) = 0.5 * dy * (p(x - 1, y) - p(x, y));
	Sp_y(x, y) = 0.5 * dx * (p(x, y) - p(x, y + 1));

	// TOP LEFT CORNER
	x = 0;
	y = ny - 1;

	mflux = rho_face(x, y).east * vel_face(x, y).east;
	mflux_face(x, y).east = 0.5 * (abs(mflux) - mflux);
	m_solver.Ae(x, y) = -dy * mflux_face(x, y).east - vis_face(x, y).east * dyx - vis(x, y) * dyx / 3;

	m_solver.Aw(x, y) = 0.0;

	m_solver.An(x, y) = 0.0;

	mflux = rho_face(x, y).south * vel_face(x, y).south;
	mflux_face(x, y).south = 0.5 * (abs(mflux) - mflux);	
	m_solver.As(x, y) = -dx * mflux_face(x, y).south - vis_face(x, y).south * dxy;

	m_solver.Ao(x, y) = dy * mflux_face(x, y).east + dx * mflux_face(x, y).south + 
						vis_face(x, y).east * dyx + vis_face(x, y).south * dxy +
						3 * vis(x, y) * dyx;

	Sp_x(x, y) = 0.5 * dy * (p(x, y) - p(x + 1, y));
	Sp_y(x, y) = 0.5 * dx * (p(x, y - 1) - p(x, y)) + u_lid * vis(x, y) * dxy;

	// TOP RIGHT CORNER
	x = nx - 1;
	y = ny - 1;

	m_solver.Ae(x, y) = 0.0;

	mflux = rho_face(x, y).west * vel_face(x, y).west;
	mflux_face(x, y).west = 0.5 * (abs(mflux) - mflux);
	m_solver.Aw(x, y) = -dy * mflux_face(x, y).west - vis_face(x, y).west * dyx + vis(x, y) * dyx / 3;

	m_solver.An(x, y) = 0.0;

	mflux = rho_face(x, y).south * vel_face(x, y).south;
	mflux_face(x, y).south = 0.5 * (abs(mflux) - mflux);	
	m_solver.As(x, y) = -dx * mflux_face(x, y).south - vis_face(x, y).south * dxy;

	m_solver.Ao(x, y) = dy * mflux_face(x, y).west + dx * mflux_face(x, y).south + 
						vis_face(x, y).west * dyx + vis_face(x, y).south * dxy +
						3 * vis(x, y) * dyx;

	Sp_x(x, y) = 0.5 * dy * (p(x - 1, y) - p(x, y));
	Sp_y(x, y) = 0.5 * dx * (p(x, y - 1) - p(x, y)) + u_lid * vis(x, y) * dxy;

	// BOTTOM FACE
	y = 0;
	for (x = 1; x < nx - 1; x++)
	{
		mflux = rho_face(x, y).east * vel_face(x, y).east;
		mflux_face(x, y).east = 0.5 * (abs(mflux) - mflux);
		m_solver.Ae(x, y) = -dy * mflux_face(x, y).east - vis_face(x, y).east * dyx;

		mflux = rho_face(x, y).west * vel_face(x, y).west;
		mflux_face(x, y).west = 0.5 * (abs(mflux) - mflux);
		m_solver.Aw(x, y) = -dy * mflux_face(x, y).west - vis_face(x, y).west * dyx;

		mflux = rho_face(x, y).north * vel_face(x, y).north;
		mflux_face(x, y).north = 0.5 * (abs(mflux) - mflux);
		m_solver.An(x, y) = -dx * mflux_face(x, y).north - vis_face(x, y).north * dxy - vis(x, y) * dxy / 3;
		
		m_solver.As(x, y) = 0.0;

		m_solver.Ao(x, y) = dy * mflux_face(x, y).east + 
							dy * mflux_face(x, y).west + 
							dx * mflux_face(x, y).north + 
							vis_face(x, y).north * dxy + 
							3 * vis(x, y) * dxy + 
							vis_face(x, y).east * dyx + 
							vis_face(x, y).west * dyx;

		Sp_x(x, y) = 0.5 * dy * (p(x - 1, y) - p(x + 1, y));
		Sp_y(x, y) = 0.5 * dx * (p(x, y) - p(x, y + 1));
	}

	// TOP FACE
	y = ny - 1;
	for (x = 1; x < nx - 1; x++)
	{
		mflux = rho_face(x, y).east * vel_face(x, y).east;
		mflux_face(x, y).east = 0.5 * (abs(mflux) - mflux);
		m_solver.Ae(x, y) = -dy * mflux_face(x, y).east - vis_face(x, y).east * dyx;

		mflux = rho_face(x, y).west * vel_face(x, y).west;
		mflux_face(x, y).west = 0.5 * (abs(mflux) - mflux);
		m_solver.Aw(x, y) = -dy * mflux_face(x, y).west - vis_face(x, y).west * dyx;

		m_solver.An(x, y) = 0.0;
		
		mflux = rho_face(x, y).south * vel_face(x, y).south;
		mflux_face(x, y).south = 0.5 * (abs(mflux) - mflux);	
		m_solver.As(x, y) = -dx * mflux_face(x, y).south - vis_face(x, y).south * dxy;

		m_solver.Ao(x, y) = dy * mflux_face(x, y).east + 
							dy * mflux_face(x, y).west + 
							dx * mflux_face(x, y).south + 
							vis_face(x, y).south * dxy +
							vis_face(x, y).east * dyx + 
							vis_face(x, y).west * dyx;

		Sp_x(x, y) = 0.5 * dy * (p(x - 1, y) - p(x + 1, y));
		Sp_y(x, y) = 0.5 * dx * (p(x, y - 1) - p(x, y)) + u_lid * vis(x, y) * dxy;
	}

	// LEFT FACE
	x = 0;
	for (y = 1; y < ny - 1; y++)
	{
		mflux = rho_face(x, y).east * vel_face(x, y).east;
		mflux_face(x, y).east = 0.5 * (abs(mflux) - mflux);
		m_solver.Ae(x, y) = -dy * mflux_face(x, y).east - vis_face(x, y).east * dyx - vis(x, y) * dyx / 3;

		m_solver.Aw(x, y) = 0.0;
	
		mflux = rho_face(x, y).north * vel_face(x, y).north;
		mflux_face(x, y).north = 0.5 * (abs(mflux) - mflux);
		m_solver.An(x, y) = -dx * mflux_face(x, y).north - vis_face(x, y).north * dxy;

		mflux = rho_face(x, y).south * vel_face(x, y).south;
		mflux_face(x, y).south = 0.5 * (abs(mflux) - mflux);
		m_solver.As(x, y) = -dx * mflux_face(x, y).south - vis_face(x, y).south * dxy;

		m_solver.Ao(x, y) = dy * mflux_face(x, y).east + 
							dx * mflux_face(x, y).north + 
							dx * mflux_face(x, y).south + 
							vis_face(x, y).south * dxy +
							vis_face(x, y).north * dxy + 
							vis_face(x, y).east * dyx + 
							3 * vis(x, y) * dyx;

		Sp_x(x, y) = 0.5 * dy * (p(x, y) - p(x + 1, y));
		Sp_y(x, y) = 0.5 * dx * (p(x, y - 1) - p(x, y + 1));
	}

	// RIGHT FACE
	x = nx - 1;
	for (y = 1; y < ny - 1; y++)
	{
		m_solver.Ae(x, y) = 0.0;

		mflux = rho_face(x, y).west * vel_face(x, y).west;
		mflux_face(x, y).west = 0.5 * (abs(mflux) - mflux);
		m_solver.Aw(x, y) = -dy * mflux_face(x, y).west - vis_face(x, y).west * dyx + vis(x, y) * dyx / 3;
		
		mflux = rho_face(x, y).north * vel_face(x, y).north;
		mflux_face(x, y).north = 0.5 * (abs(mflux) - mflux);
		m_solver.An(x, y) = -dx * mflux_face(x, y).north - vis_face(x, y).north * dxy;

		mflux = rho_face(x, y).south * vel_face(x, y).south;
		mflux_face(x, y).south = 0.5 * (abs(mflux) - mflux);
		m_solver.As(x, y) = -dx * mflux_face(x, y).south - vis_face(x, y).south * dxy;

		m_solver.Ao(x, y) = dy * mflux_face(x, y).west + 
							dx * mflux_face(x, y).north + 
							dx * mflux_face(x, y).south + 
							vis_face(x, y).south * dxy +
							vis_face(x, y).north * dxy + 
							vis_face(x, y).west * dyx - 
							3 * vis(x, y) * dyx;

		Sp_x(x, y) = 0.5 * dy * (p(x - 1, y) - p(x, y));
		Sp_y(x, y) = 0.5 * dx * (p(x, y - 1) - p(x, y + 1));
	}


	// INTERIOR NODES
	for (x = 1; x < nx - 1; ++x)
	{
		for (y = 1; y < ny - 1; ++y)
		{
			mflux = rho_face(x, y).east * vel_face(x, y).east;
			mflux_face(x, y).east = 0.5 * (abs(mflux) - mflux);
			m_solver.Ae(x, y) = -dy * mflux_face(x, y).east - vis_face(x, y).east * dyx;

			mflux = rho_face(x, y).west * vel_face(x, y).west;
			mflux_face(x, y).west = 0.5 * (abs(mflux) - mflux);
			m_solver.Aw(x, y) = -dy * mflux_face(x, y).west - vis_face(x, y).west * dyx;
			
			mflux = rho_face(x, y).north * vel_face(x, y).north;
			mflux_face(x, y).north = 0.5 * (abs(mflux) - mflux);
			m_solver.An(x, y) = -dx * mflux_face(x, y).north - vis_face(x, y).north * dxy;

			mflux = rho_face(x, y).south * vel_face(x, y).south;
			mflux_face(x, y).south = 0.5 * (abs(mflux) - mflux);
			m_solver.As(x, y) = -dx * mflux_face(x, y).south - vis_face(x, y).south * dxy;

			m_solver.Ao(x, y) = mflux_face(x, y).east + mflux_face(x, y).west +
								mflux_face(x, y).north +  mflux_face(x, y).south + 
								vis_face(x, y).east * dyx + vis_face(x, y).west * dyx +
								vis_face(x, y).north * dxy + vis_face(x, y).south * dxy;

			Sp_x(x, y) = 0.5 * dy * (p(x - 1, y) - p(x + 1, y));
			Sp_y(x, y) = 0.5 * dx * (p(x, y - 1) - p(x, y + 1));
		}
	}

	// END OF MESH
}
