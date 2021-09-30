// admin.h

#ifndef _ADMIN_H
#define _ADMIN_H

#include <cmath>
#include <iostream>

extern double tijd; // JW wordt eigenlijk alleen gebruikt voor logging (en q_in interpolatie)
extern double dz; // JW per definitie H/Npz
extern double Av;
extern double S;
extern double H;
extern double L;
extern double dx;
extern double dt;
//JW extern double q_in;
extern std::ofstream outlog;

namespace admin{
	extern int DebugOutput; //2012 10 03 Olav

	/*
	 * General parameters
     */
	extern int Npx;
		// number of grid points in horizontal x direction (120 was original)
	extern double dtr; // OLAV: 1.0 was original
		// time step during a simulation (s)
	extern double dt_write; // OLAV: 15.0 was original
		// time step at which data should be written (s)
    extern double tend;
		// end time of a simulation (s)
    
	extern double ampbeds_factor; //Default: 0.1 OLAV: 2012 12 18
		// factor to determine initial amplitude of a  sinusoidal bed disturbance
	extern int AllowFlowSep; //OLAV: 2012 09 14
	//const int AllowOneFSZ  =1; //OLAV: 2013 04 16 
	extern int AllowAvalanching; //OLAV: 2013 04 16
	extern int SimpleLength; //OLAV: 2012 09 06
                          //0: do stability anlysis for dune length (original)
                          //1: don't do stab, assume length = factor*H
    extern double SimpleLengthFactor; //For when SimpleLength =1
    extern int numStab; //For when SimpleLength = 0
	    //OLAV: 10=default
		// number of intervals for stability analysis
	extern int Hifactor; //For when SimpleLength = 0
	    //OLAV: added 2011 04 01 (10=default)
        // factor with which H is multiplied to find 
        // maximum wave length for stability analysis
	extern double Hcrit_global; //OLAV 2013 02 01 (was 0.25)
    
    extern int transport_eq; //OLAV 2011 02 28
        // flag for which sediment transport equation to use 
        // 1: MPM (original)
        // 2: N&T
        // 3: MPM with linear relaxation equation (LRE)
	extern int alpha_varies; //for N&T and LRE
								   //0: alpha = constant
								   //1: Sekine&Kikkawa (with N&T)
								   //2: Shimizu et al. (with N&T)
    extern double alpha_lag; // For alpha_varies=0 with LRE or N&T
	extern int moeilijkdoen; //OLAV: different guess flux(0) method for lag with linear relaxation
    extern double correction_NT; // correction factor sediment flux N&T (TEST OLAV 2011 03 17)
	extern int Npsl_min; //for N&T
	extern int stle_factor; //for N&T
	
    /* 
	 * Simulation dependent parameters
	 */
	extern double q_in1;// OLAV: 0.07716 was original (0.32)
		// constant discharge in simulation (m^2/s)
	extern double H0; // OLAV: 0.152 was original
		// initial flow depth (m)
	extern double ii;
		// bed slope (m/m, -)~  	 	 12.e-4;
	extern double D50; // OLAV: 0.50e-3 was original
		// (uniform) grain size
	extern double thetacr; //N&T: 0.035 //was 0.05 //should be 0.040 for D=0.28mm, 0.0391 for D=0.29mm
		// critical Shields parameter [-]
	extern double dts;
		// time step for stability analysis (s)
	extern int Npz;
		// number of grid points in vertical z direction
	extern int nd;
		// number of dunes in domain (-)			
	extern int readbed;
		// 1: read bed from file (inp_bottom.inp)
		// 0: start from scratch 
	extern int readfw;
		// 1: read floodwave from file (floodwave.inp)
		// 0: use q_in1

	/*
	 * Numerical parameters
	 */
	extern double sepcritangle;
		// bed angle at which flow separation sets in (degrees) -10
	const double grad_2_deg = 360./(2.*M_PI);
		// help variable for grad to degress
	extern double g;
		// acceleration of gravity
	extern double F;
		// Forcing term in the momentum equation	
	extern double kappa;
		// Von Karman constant
	extern double tt; // default is 100
		// splits up the timestep in tt parts, to make the transport calculations more stable
	extern double tresh;
		// threshold accuracy for flow solver
	extern int max_it;
		// maximum number of iterations for flow solver

	/*
	* Water (behaviour) parameters
	*/
	extern double denswater;
        // density of water [kg/m3] OLAV 2011 02 24	
	extern double nu;
	    // kinematic viscosity of water [m2/s]
	extern double BETA1;
	extern double BETA2;
		// 17-10-2016 RD, toegevoegd om te kunnen kalibreren met BETA1 en BETA2
			
	/*
	 * Sediment (behaviour) parameters
	 */
    extern double denssand;
        // density of sand [kg/m3] OLAV 2011 02 24
    extern double sgsand;
        // specific gravity of sand [-] OLAV 2011 02 24
    extern double delta;
        // delta parameter [-] OLAV 2011 02 24 
    extern double ampbeds; //Default: 0.1*D50
		// initial amplitude of a sinusoidal bed disturbance
	extern double epsilonp;
		// porosity parameter
	extern double repose;
		// angle of repose in degrees
	extern double m;
		// coefficient of alpha
	extern double alpha; // changed OLAV 2011 02 24
		// alpha in the sediment transport equation
	extern double be;
		// power of theta in the sediment transport equation
	extern double l1;
		// parameter in sediment transport equation (without critical shear)
	extern double l2; // = 1/tan(-repose/grad_2_deg)
		// parameter in sediment transport equation
    extern double F0 ; //0.144
	extern double F0_dim;
        // parameter and corrected parameter in N&T sediment transport equation [-] OLAV 2011 02 24
   	extern double meanstle;
        // mean step length in lag sediment transport equation [-] OLAV 2011 02 24
	extern double A2_geom; //
	    // sediment geometry factor (N&T)
    extern double A3_geom; //
	    // sediment geometry factor (N&T)
    extern double k2;
        // parameter in N&T sediment transport equation [-] OLAV 2011 04 11
		
	//Sekine & Kikkawa parameters
	extern double alpha_2;
	extern double D_star;
	extern double w_s;
	extern double u_star_cr;
	extern double alpha_min_SK;
	extern double alpha_max_SK;
	
	//Shimizu et al. step length parameters
	extern double alpha_min_S;
	extern double alpha_max_S;
	extern double theta_min_S;
	extern double theta_max_S;
	extern double H_ref;
	extern int keepsgrowing;
	
	/*
	 * Array defintion
	 */
//	extern int nt;
//	extern int nf;
//	extern int nf2;
//	int o(int j_ex,int i_ex,int v);
	int o2(int i_ex);
//	int o3(int i_in);
}
#endif
