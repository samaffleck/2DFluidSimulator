#pragma once

#include "Eigen/Dense"

class ODESolver
{
public:
	ODESolver() = default;
	~ODESolver() = default;

	void initialise(int numberOfXNodes, int numberOfYNodes);
	void solve(Eigen::MatrixXd& var);

	Eigen::MatrixXd Ao{};
	Eigen::MatrixXd Ae{};
	Eigen::MatrixXd Aw{};
	Eigen::MatrixXd An{};
	Eigen::MatrixXd As{};
	Eigen::MatrixXd S{};

	double leftBoundary{};
	double rightBoundary{};
	double topBoundary{};
	double bottomBoundary{};

private:
	double m_tolerance = 1e-6;
	int nx, ny;
	double ve, vw, vn, vs;
	
private:
	void update(Eigen::MatrixXd& var);
	double getResidual(const Eigen::MatrixXd& var);
	void setNeighbourCells(const Eigen::MatrixXd& var, int x, int y);

};
