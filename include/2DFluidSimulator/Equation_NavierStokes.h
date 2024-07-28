#include "IEquation.h"

class Equation_NavierStokes : public IEquation
{
public:

	void initialiseEquation(int numberOfXCells, int numberOfYCells) override;
	void update() override;

	// Getters
	Eigen::VectorXd& getU_x() { return u_x; }
	Eigen::VectorXd& getU_y() { return u_y; }
	Eigen::VectorXd& getP() { return p; }
	double& U_x(int x, int y) { return u_x(x + m_numberOfXCells * y); }
	double& U_y(int x, int y) { return u_y(x + m_numberOfXCells * y); }
	double& P(int x, int y) { return p(x + m_numberOfXCells * y); }

private:
	Eigen::VectorXd u_x{};		// Cell center vecloity in x-direction
	Eigen::VectorXd u_face_x{};	// Cell face velocity in x-direction
	Eigen::VectorXd u_y{};		// Cell center vecloity in y-direction
	Eigen::VectorXd u_face_y{};	// Cell face velocity in y-direction
	Eigen::VectorXd p{};		// Absolute pressure
	Eigen::VectorXd p_corr{};	// Pressure correction

	int m_numberOfXCells{};
	int m_numberOfYCells{};

private:
	void updateLinkCoefficient();

};

