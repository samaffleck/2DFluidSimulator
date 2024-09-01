#pragma once

#include "Eigen/Dense"

class ODESolver
{
public:
	ODESolver() = default;
	~ODESolver() = default;

	void initialise(int numberOfXNodes, int numberOfYNodes);
	void solve(Eigen::MatrixXd& var,
		const Eigen::MatrixXd& Ao,
		const Eigen::MatrixXd& Ae,
		const Eigen::MatrixXd& Aw,
		const Eigen::MatrixXd& An,
		const Eigen::MatrixXd& As,
		const Eigen::MatrixXd& S,
		double damping_factor = 0.0);
	double getResidual(const Eigen::MatrixXd& var,
		const Eigen::MatrixXd& Ao,
		const Eigen::MatrixXd& Ae,
		const Eigen::MatrixXd& Aw,
		const Eigen::MatrixXd& An,
		const Eigen::MatrixXd& As,
		const Eigen::MatrixXd& S);

	double leftBoundary{};
	double rightBoundary{};
	double topBoundary{};
	double bottomBoundary{};

private:
	double m_tolerance = 1e-3;
	int nx{}, ny{};
	double ve{}, vw{}, vn{}, vs{};
	int m_damping_factor{};
	
private:
	void update(Eigen::MatrixXd& var, 
		const Eigen::MatrixXd& Ao,
		const Eigen::MatrixXd& Ae,
		const Eigen::MatrixXd& Aw,
		const Eigen::MatrixXd& An,
		const Eigen::MatrixXd& As,
		const Eigen::MatrixXd& S);
	void setNeighbourCells(const Eigen::MatrixXd& var, 
		int x, 
		int y);

};
