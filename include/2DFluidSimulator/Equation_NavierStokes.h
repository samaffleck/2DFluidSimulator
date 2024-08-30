#include "IEquation.h"

class Equation_NavierStokes : public IEquation
{
public:

	void initialiseEquation(int numberOfXCells, int numberOfYCells) override;
	void update() override;

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

	// Resuduals
	double res_u;
	double res_v;
	double res_p;

	// Tuning parameters
	double p_scalar = 0.2;
	double vel_scalar = 0.8;

	// Other
	int nx{};
	int ny{};
	double u_lid{};

private:
	void updateLinkCoefficient();
	void updateFaceVelocities();
	void updatePressureLinks();
	void updateVelocityCorrection();
	void correctVelocityAndPressure();

};

