/*
 * BedConfig.h
 *
 *  Created on: 24 nov. 2021
 *      Author: impor
 */

#ifndef SRC_BEDCONFIG_H_
#define SRC_BEDCONFIG_H_


class Config;

class BedConfig {
public:
	const int DebugOutput;
	const int Npx;			// number of grid points in horizontal x direction
	const int AllowAvalanching;
	const int transport_eq;	// flag for which sediment transport equation to use
	// 1: MPM (original)
	// 2: N&T
	// 3: MPM with linear relaxation equation (LRE)
	const int alpha_varies;	//for N&T and LRE
	//0: alpha = constant
	//1: Sekine&Kikkawa (with N&T)
	//2: Shimizu et al. (with N&T)
	//3: Shimizu et al., adjusted by van Duin et al. 2021 (with N&T)
	const double alpha_lag;	// For alpha_varies=0 with LRE or N&T
	const int moeilijkdoen;	// different guess flux(0) method for lag with linear relaxation
	const int Npsl_min;		//for N&T
	const int stle_factor;	//for N&T
	const double D50;		// (uniform) grain size
	const double thetacr;	// N&T: critical Shields parameter
	//should be 0.040 for D=0.28mm, 0.0391 for D=0.29mm
	const int nd;			// number of dunes in domain
	const double sepcritangle;// bed angle at which flow separation sets in
	const double g;			// acceleration of gravity
	const double ii;		// global bed slope
	const double tt;		// splits up the timestep in tt parts, to make the transport calculations more stable
	const double delta;		// delta parameter
	const double epsilonp;	// porosity parameter
	const double repose;	// angle of repose
	const double alpha;		// alpha in the sediment transport equation
	const double be;		// power of theta in the sediment transport equation
	const double l2;		// parameter in sediment transport equation
	// = 1/tan(-repose/grad_2_deg)
	const double F0_dim;	// parameter and corrected parameter in N&T sediment transport equation
	const double meanstle;	// mean step length in lag sediment transport equation

	//Sekine & Kikkawa parameters
	const double alpha_2;
	const double D_star;
	const double w_s;
	double u_star_cr;
	const double alpha_min_SK;
	const double alpha_max_SK;

	//Shimizu et al. step length parameters
	const double alpha_min_S;
	const double alpha_max_S;
	const double theta_min_S;
	const double theta_max_S;
	const double H_ref;
//	const int keepsgrowing; //

	BedConfig(const Config& cfg);
	BedConfig(const BedConfig&) = default;
private:
	BedConfig() = delete;

};

#endif /* SRC_BEDCONFIG_H_ */
