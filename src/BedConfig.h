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
	int DebugOutput;
	int Npx;			// number of grid points in horizontal x direction
	int AllowAvalanching;
	int transport_eq;	// flag for which sediment transport equation to use
	// 1: MPM (original)
	// 2: N&T
	// 3: MPM with linear relaxation equation (LRE)
	int alpha_varies;	//for N&T and LRE
	//0: alpha = constant
	//1: Sekine&Kikkawa (with N&T)
	//2: Shimizu et al. (with N&T)
	double alpha_lag;	// For alpha_varies=0 with LRE or N&T
	int moeilijkdoen;	// different guess flux(0) method for lag with linear relaxation
	int Npsl_min;		//for N&T
	int stle_factor;	//for N&T
	double D50;			// (uniform) grain size
	double thetacr;		// N&T: critical Shields parameter
	//should be 0.040 for D=0.28mm, 0.0391 for D=0.29mm
	int nd;				// number of dunes in domain
	double sepcritangle;// bed angle at which flow separation sets in
	double g;			// acceleration of gravity
	double tt;			// splits up the timestep in tt parts, to make the transport calculations more stable
	double delta;		// delta parameter
	double epsilonp;	// porosity parameter
	double repose;		// angle of repose
	double alpha;		// alpha in the sediment transport equation
	double be;			// power of theta in the sediment transport equation
	double l2;			// parameter in sediment transport equation
	// = 1/tan(-repose/grad_2_deg)
	double F0_dim;		// parameter and corrected parameter in N&T sediment transport equation
	double meanstle;		// mean step length in lag sediment transport equation

	//Sekine & Kikkawa parameters
	double alpha_2;
	double D_star;
	double w_s;
	double u_star_cr;
	double alpha_min_SK;
	double alpha_max_SK;

	//Shimizu et al. step length parameters
	double alpha_min_S;
	double alpha_max_S;
	double theta_min_S;
	double theta_max_S;
	double H_ref;
	int keepsgrowing;

	BedConfig(const Config& cfg);
	BedConfig(const BedConfig&) = default;
private:
	BedConfig() = delete;

};

#endif /* SRC_BEDCONFIG_H_ */
