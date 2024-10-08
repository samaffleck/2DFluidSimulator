#include "IEquation.h"
#include "Eigen/Dense"

class Equation_NavierStokes : public IEquation
{
public:
	void initialiseEquation(int numberOfXCells, int numberOfYCells) override;
	void update() override;
	void logData(std::string resultsDirectory) override;

private:
	// Cell center values
	Eigen::MatrixXd u{};		// Vecloity in x-direction
	Eigen::MatrixXd u_c{};		// Vecloity correction in x-direction
	Eigen::MatrixXd v{};		// Vecloity in y-direction
	Eigen::MatrixXd v_c{};		// Vecloity correction in y-direction
	Eigen::MatrixXd p{};		// Absolute pressure
	Eigen::MatrixXd p_c{};		// Pressure correction
	Eigen::MatrixXd vis{};		// Viscosity [Pa-s]
	Eigen::MatrixXd rho{};		// Density [kg/m3]
	
	// Face values accessed by: FACE_VALUE(x, y).east/west/north/south
	Eigen::MatrixX<FaceValue> vel_face{};	// Vecloity [m/s]
	Eigen::MatrixX<FaceValue> vel_c_face{};	// Vecloity correction [m/s]
	Eigen::MatrixX<FaceValue> vis_face{};	// Viscosity [Pa-s]
	Eigen::MatrixX<FaceValue> rho_face{};	// Density [kg/m3]
	Eigen::MatrixX<FaceValue> mflux_face{};	// Mass Flux [kg/(m2-s)]

	// Momentum link coefficients
	Eigen::MatrixXd Ao{};
	Eigen::MatrixXd Ae{};
	Eigen::MatrixXd Aw{};
	Eigen::MatrixXd An{};
	Eigen::MatrixXd As{};
	Eigen::MatrixXd Sp_x{};		// Pressure source in the x-direction
	Eigen::MatrixXd Sp_y{};		// Pressure source in the y-direction

	// Pressure link coefficients
	Eigen::MatrixXd Ap_o{};
	Eigen::MatrixXd Ap_e{};
	Eigen::MatrixXd Ap_w{};
	Eigen::MatrixXd Ap_n{};
	Eigen::MatrixXd Ap_s{};
	Eigen::MatrixXd Sp{};

	// Tuning parameters
	double p_relax = 0.2;
	double vel_relax = 0.8; // 0.8
	double tolerance = 1e-10;

	// Resuduals
	double res_u;
	double res_v;
	double res_p;

	// Max residuals
	double res_u_max = tolerance;
	double res_v_max = tolerance;
	double res_p_max = tolerance;

	// Other
	int nx{};
	int ny{};
	double u_lid = 1.0;

private:
	void updateLinkCoefficient();
	void updateFaceVelocities();
	void updatePressureLinks();
	void updateVelocityCorrection();
	void correctVelocityAndPressure();
	bool isConverged();
	void logVelocity();

};

