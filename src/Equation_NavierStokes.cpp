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
		updateFaceVelocities();

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


void Equation_NavierStokes::updateFaceVelocities()
{
	// Using PWIM Equations

	// Get constants
	double dx = p_mesh->getCells()[0][0].dx;
	double dy = p_mesh->getCells()[0][0].dy;

	// 1st LAYER
	// BOTTOM LEFT CORNER (1st LAYER)
	int x = 0;
	int y = 0;

	double ao_o = 1 / m_solver.Ao(x, y);
	double ao_e = 1 / m_solver.Ao(x + 1, y);
	double ao_w = 1 / m_solver.Ao(x - 1, y);
	double ao_n = 1 / m_solver.Ao(x, y + 1);
	double ao_s = 1 / m_solver.Ao(x, y - 1);

	vel_face(x, y).east = 0.5 * (u(x, y) + u(x + 1, y)) +
		0.25 * dy * ao_o * (p(x + 1, y) - p(x, y)) +
		0.25 * dy * ao_e * (p(x + 2, y) - p(x, y)) -
		0.5 * dy * (ao_o + ao_e) * (p(x + 1, y) - p(x, y));

	vel_face(x, y).west = 0.0;

	vel_face(x, y).north = 0.5 * (v(x, y) + v(x, y + 1)) +
		0.25 * dx * ao_o * (p(x, y + 1) - p(x, y)) +
		0.25 * dx * ao_n * (p(x, y + 2) - p(x, y)) -
		0.5 * dx * (ao_o + ao_n) * (p(x, y + 1) - p(x, y));

	vel_face(x, y).south = 0.0;

	// x = 1, y = 0
	x = 1;
	y = 0;

	ao_o = 1 / m_solver.Ao(x, y);
	ao_e = 1 / m_solver.Ao(x + 1, y);
	ao_w = 1 / m_solver.Ao(x - 1, y);
	ao_n = 1 / m_solver.Ao(x, y + 1);
	ao_s = 1 / m_solver.Ao(x, y - 1);




	// x = 0, y = 1
	x = 0;
	y = 1;

	ao_o = 1 / m_solver.Ao(x, y);
	ao_e = 1 / m_solver.Ao(x + 1, y);
	ao_w = 1 / m_solver.Ao(x - 1, y);
	ao_n = 1 / m_solver.Ao(x, y + 1);
	ao_s = 1 / m_solver.Ao(x, y - 1);





	// BOTTOM RIGHT CORNER (1st LAYER)
	x = nx - 1;
	y = 0;

	ao_o = 1 / m_solver.Ao(x, y);
	ao_e = 1 / m_solver.Ao(x + 1, y);
	ao_w = 1 / m_solver.Ao(x - 1, y);
	ao_n = 1 / m_solver.Ao(x, y + 1);
	ao_s = 1 / m_solver.Ao(x, y - 1);

	vel_face(x, y).east = 0.0;

	vel_face(x, y).west = 0.5 * (u(x, y) + u(x - 1, y)) +
		0.25 * dy * ao_o * (p(x, y) - p(x - 1, y)) +
		0.25 * dy * ao_w * (p(x, y) - p(x - 2, y)) -
		0.5 * dy * (ao_o + ao_w) * (p(x, y) - p(x - 1, y));

	vel_face(x, y).north = 0.5 * (v(x, y) + v(x, y + 1)) +
		0.25 * dx * ao_o * (p(x, y + 1) - p(x, y)) +
		0.25 * dx * ao_n * (p(x, y + 2) - p(x, y)) -
		0.5 * dx * (ao_o + ao_n) * (p(x, y + 1) - p(x, y));

	vel_face(x, y).south = 0.0;

	// x = nx - 2, y = 0
	x = nx - 2;
	y = 0;

	ao_o = 1 / m_solver.Ao(x, y);
	ao_e = 1 / m_solver.Ao(x + 1, y);
	ao_w = 1 / m_solver.Ao(x - 1, y);
	ao_n = 1 / m_solver.Ao(x, y + 1);
	ao_s = 1 / m_solver.Ao(x, y - 1);



	// x = nx - 1, y = 1
	x = nx - 1;
	y = 1;

	ao_o = 1 / m_solver.Ao(x, y);
	ao_e = 1 / m_solver.Ao(x + 1, y);
	ao_w = 1 / m_solver.Ao(x - 1, y);
	ao_n = 1 / m_solver.Ao(x, y + 1);
	ao_s = 1 / m_solver.Ao(x, y - 1);





	// TOP LEFT CORNER (1st LAYER)
	x = 0;
	y = ny - 1;

	ao_o = 1 / m_solver.Ao(x, y);
	ao_e = 1 / m_solver.Ao(x + 1, y);
	ao_w = 1 / m_solver.Ao(x - 1, y);
	ao_n = 1 / m_solver.Ao(x, y + 1);
	ao_s = 1 / m_solver.Ao(x, y - 1);

	vel_face(x, y).east = 0.5 * (u(x, y) + u(x + 1, y)) +
		0.25 * dy * ao_o * (p(x + 1, y) - p(x, y)) +
		0.25 * dy * ao_e * (p(x + 2, y) - p(x, y)) -
		0.5 * dy * (ao_o + ao_e) * (p(x + 1, y) - p(x, y));

	vel_face(x, y).west = 0.0;

	vel_face(x, y).north = 0.0;

	vel_face(x, y).south = 0.5 * (v(x, y) + v(x, y - 1)) +
		0.25 * dx * ao_o * (p(x, y) - p(x, y - 1)) +
		0.25 * dx * ao_s * (p(x, y) - p(x, y - 2)) -
		0.5 * dx * (ao_o + ao_s) * (p(x, y) - p(x, y - 1));

	// x = 0, y = ny - 2
	x = 0;
	y = ny - 2;

	ao_o = 1 / m_solver.Ao(x, y);
	ao_e = 1 / m_solver.Ao(x + 1, y);
	ao_w = 1 / m_solver.Ao(x - 1, y);
	ao_n = 1 / m_solver.Ao(x, y + 1);
	ao_s = 1 / m_solver.Ao(x, y - 1);



	// x = 1, y = ny - 1
	x = 1;
	y = ny - 1;

	ao_o = 1 / m_solver.Ao(x, y);
	ao_e = 1 / m_solver.Ao(x + 1, y);
	ao_w = 1 / m_solver.Ao(x - 1, y);
	ao_n = 1 / m_solver.Ao(x, y + 1);
	ao_s = 1 / m_solver.Ao(x, y - 1);





	// TOP RIGHT CORNER (1st LAYER)
	x = nx - 1;
	y = ny - 1;

	ao_o = 1 / m_solver.Ao(x, y);
	ao_e = 1 / m_solver.Ao(x + 1, y);
	ao_w = 1 / m_solver.Ao(x - 1, y);
	ao_n = 1 / m_solver.Ao(x, y + 1);
	ao_s = 1 / m_solver.Ao(x, y - 1);

	vel_face(x, y).east = 0.0;

	vel_face(x, y).west = 0.5 * (u(x, y) + u(x - 1, y)) +
		0.25 * dy * ao_o * (p(x, y) - p(x - 1, y)) +
		0.25 * dy * ao_w * (p(x, y) - p(x - 2, y)) -
		0.5 * dy * (ao_o + ao_w) * (p(x, y) - p(x - 1, y));

	vel_face(x, y).north = 0.0;

	vel_face(x, y).south = 0.5 * (v(x, y) + v(x, y - 1)) +
		0.25 * dx * ao_o * (p(x, y) - p(x, y - 1)) +
		0.25 * dx * ao_s * (p(x, y) - p(x, y - 2)) -
		0.5 * dx * (ao_o + ao_s) * (p(x, y) - p(x, y - 1));

	// x = nx - 2, y = ny - 1
	x = nx - 2;
	y = ny - 1;

	ao_o = 1 / m_solver.Ao(x, y);
	ao_e = 1 / m_solver.Ao(x + 1, y);
	ao_w = 1 / m_solver.Ao(x - 1, y);
	ao_n = 1 / m_solver.Ao(x, y + 1);
	ao_s = 1 / m_solver.Ao(x, y - 1);



	// x = nx - 1, y = ny - 2
	x = nx - 1;
	y = ny - 2;

	ao_o = 1 / m_solver.Ao(x, y);
	ao_e = 1 / m_solver.Ao(x + 1, y);
	ao_w = 1 / m_solver.Ao(x - 1, y);
	ao_n = 1 / m_solver.Ao(x, y + 1);
	ao_s = 1 / m_solver.Ao(x, y - 1);





	// BOTTOM FACE (1st LAYER)
	y = 0;
	for (x = 2; x < nx - 2; x++)
	{
		ao_o = 1 / m_solver.Ao(x, y);
		ao_e = 1 / m_solver.Ao(x + 1, y);
		ao_w = 1 / m_solver.Ao(x - 1, y);
		ao_n = 1 / m_solver.Ao(x, y + 1);
		ao_s = 1 / m_solver.Ao(x, y - 1);


	}

	// TOP FACE (1st LAYER)
	y = ny - 1;
	for (x = 2; x < nx - 2; x++)
	{
		ao_o = 1 / m_solver.Ao(x, y);
		ao_e = 1 / m_solver.Ao(x + 1, y);
		ao_w = 1 / m_solver.Ao(x - 1, y);
		ao_n = 1 / m_solver.Ao(x, y + 1);
		ao_s = 1 / m_solver.Ao(x, y - 1);


	}

	// LEFT FACE (1st LAYER)
	x = 0;
	for (y = 2; y < ny - 2; y++)
	{
		ao_o = 1 / m_solver.Ao(x, y);
		ao_e = 1 / m_solver.Ao(x + 1, y);
		ao_w = 1 / m_solver.Ao(x - 1, y);
		ao_n = 1 / m_solver.Ao(x, y + 1);
		ao_s = 1 / m_solver.Ao(x, y - 1);


	}

	// RIGHT FACE (1st LAYER)
	x = nx - 1;
	for (y = 2; y < ny - 2; y++)
	{
		ao_o = 1 / m_solver.Ao(x, y);
		ao_e = 1 / m_solver.Ao(x + 1, y);
		ao_w = 1 / m_solver.Ao(x - 1, y);
		ao_n = 1 / m_solver.Ao(x, y + 1);
		ao_s = 1 / m_solver.Ao(x, y - 1);


	}

	// 2nd LAYER
	// BOTTOM LEFT CORNER (2nd LAYER)
	int x = 1;
	int y = 1;

	ao_o = 1 / m_solver.Ao(x, y);
	ao_e = 1 / m_solver.Ao(x + 1, y);
	ao_w = 1 / m_solver.Ao(x - 1, y);
	ao_n = 1 / m_solver.Ao(x, y + 1);
	ao_s = 1 / m_solver.Ao(x, y - 1);


	// BOTTOM RIGHT CORNER (2nd LAYER)
	x = nx - 2;
	y = 1;

	ao_o = 1 / m_solver.Ao(x, y);
	ao_e = 1 / m_solver.Ao(x + 1, y);
	ao_w = 1 / m_solver.Ao(x - 1, y);
	ao_n = 1 / m_solver.Ao(x, y + 1);
	ao_s = 1 / m_solver.Ao(x, y - 1);


	// TOP LEFT CORNER (2nd LAYER)
	x = 1;
	y = ny - 2;

	ao_o = 1 / m_solver.Ao(x, y);
	ao_e = 1 / m_solver.Ao(x + 1, y);
	ao_w = 1 / m_solver.Ao(x - 1, y);
	ao_n = 1 / m_solver.Ao(x, y + 1);
	ao_s = 1 / m_solver.Ao(x, y - 1);


	// TOP RIGHT CORNER (2nd LAYER)
	x = nx - 2;
	y = ny - 2;

	ao_o = 1 / m_solver.Ao(x, y);
	ao_e = 1 / m_solver.Ao(x + 1, y);
	ao_w = 1 / m_solver.Ao(x - 1, y);
	ao_n = 1 / m_solver.Ao(x, y + 1);
	ao_s = 1 / m_solver.Ao(x, y - 1);


	// BOTTOM FACE (2nd LAYER)
	y = 1;
	for (x = 2; x < nx - 2; x++)
	{
		ao_o = 1 / m_solver.Ao(x, y);
		ao_e = 1 / m_solver.Ao(x + 1, y);
		ao_w = 1 / m_solver.Ao(x - 1, y);
		ao_n = 1 / m_solver.Ao(x, y + 1);
		ao_s = 1 / m_solver.Ao(x, y - 1);


	}

	// TOP FACE (2nd LAYER)
	y = ny - 2;
	for (x = 2; x < nx - 2; x++)
	{
		ao_o = 1 / m_solver.Ao(x, y);
		ao_e = 1 / m_solver.Ao(x + 1, y);
		ao_w = 1 / m_solver.Ao(x - 1, y);
		ao_n = 1 / m_solver.Ao(x, y + 1);
		ao_s = 1 / m_solver.Ao(x, y - 1);


	}

	// LEFT FACE (2nd LAYER)
	x = 1;
	for (y = 2; y < ny - 2; y++)
	{
		ao_o = 1 / m_solver.Ao(x, y);
		ao_e = 1 / m_solver.Ao(x + 1, y);
		ao_w = 1 / m_solver.Ao(x - 1, y);
		ao_n = 1 / m_solver.Ao(x, y + 1);
		ao_s = 1 / m_solver.Ao(x, y - 1);


	}

	// RIGHT FACE (2nd LAYER)
	x = nx - 2;
	for (y = 2; y < ny - 2; y++)
	{
		ao_o = 1 / m_solver.Ao(x, y);
		ao_e = 1 / m_solver.Ao(x + 1, y);
		ao_w = 1 / m_solver.Ao(x - 1, y);
		ao_n = 1 / m_solver.Ao(x, y + 1);
		ao_s = 1 / m_solver.Ao(x, y - 1);


	}

	// INTERIOR NODES
	for (x = 2; x < nx - 2; ++x)
	{
		for (y = 2; y < ny - 2; ++y)
		{
			ao_o = 1 / m_solver.Ao(x, y);
			ao_e = 1 / m_solver.Ao(x + 1, y);
			ao_w = 1 / m_solver.Ao(x - 1, y);
			ao_n = 1 / m_solver.Ao(x, y + 1);
			ao_s = 1 / m_solver.Ao(x, y - 1);

			vel_face(x, y).east = 0.5 * (u(x, y) + u(x + 1, y)) +
				0.25 * dy * ao_o * (p(x + 1, y) - p(x - 1, y)) +
				0.25 * dy * ao_e * (p(x + 2, y) - p(x, y)) -
				0.5 * dy * (ao_o + ao_e) * (p(x + 1, y) - p(x, y));

			vel_face(x, y).west = 0.5 * (u(x, y) + u(x - 1, y)) +
				0.25 * dy * ao_o * (p(x + 1, y) - p(x - 1, y)) +
				0.25 * dy * ao_w * (p(x, y) - p(x - 2, y)) -
				0.5 * dy * (ao_o + ao_w) * (p(x, y) - p(x - 1, y));

			vel_face(x, y).north = 0.5 * (v(x, y) + v(x, y + 1)) +
				0.25 * dx * ao_o * (p(x, y + 1) - p(x, y - 1)) +
				0.25 * dx * ao_n * (p(x, y + 2) - p(x, y)) -
				0.5 * dx * (ao_o + ao_n) * (p(x, y + 1) - p(x, y));

			vel_face(x, y).south = 0.5 * (v(x, y) + v(x, y - 1)) +
				0.25 * dx * ao_o * (p(x, y + 1) - p(x, y - 1)) +
				0.25 * dx * ao_s * (p(x, y) - p(x, y - 2)) -
				0.5 * dx * (ao_o + ao_s) * (p(x, y) - p(x, y - 1));

		}
	}

	// END OF MESH
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
