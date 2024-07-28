#pragma once

#include "Eigen/Dense"

class ODESolver
{
public:
	ODESolver() = default;
	~ODESolver() = default;

	void initialise(int numberOfNodes);
	void solve(Eigen::VectorXd& V);

	Eigen::VectorXd Ao{};
	Eigen::VectorXd Ae{};
	Eigen::VectorXd Aw{};
	Eigen::VectorXd An{};
	Eigen::VectorXd As{};
	Eigen::VectorXd S{};

private:
	double m_tolerance = 1e-6;
	int k{}; // number of x nodes

private:
	void update(Eigen::VectorXd& V);
	double getResudule(const Eigen::VectorXd& V);

};
