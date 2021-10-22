// admin.cpp

#include "admin.h"
#include <cmath>
#include <iostream>

using namespace std;

namespace admin {

int DebugOutput = 0;

int Npx = 120;
double dtr = 240.0;
double dt_write = 240.;
// JW const double tend=100000.; // OLAV: 3.0*60.*60. was original
double tend = 3.0*60.*60.;
double ampbeds_factor = 100.1;
int AllowFlowSep = 0;
int AllowAvalanching = 1;
int SimpleLength = 0;
double SimpleLengthFactor = 7;
int numStab = 30;
int Hifactor = 15;
double Hcrit_global = 5;
int transport_eq = 2;
int alpha_varies = 2;
double alpha_lag =100.0;
int moeilijkdoen = 0;
double correction_NT = 1.0;
int Npsl_min = 40;
int stle_factor =5;

double q_in1 = 6.25;
double H0 = 9.1;
double ii = 1.1e-4;
double D50 = 0.18e-3;
double thetacr = 0.053;
double dts = 0.01;
int Npz = 25;
int nd=1;
int readbed=0;
int readfw=1;

double sepcritangle = -10.0;
double g = 9.81;
double F = g*ii;
double kappa = 0.407;
double tt = 100.0;
double tresh = 1e-8;
int max_it=8;

double denswater = 1000.0;
double nu = 1e-6;
double BETA1 = 0.245;
double BETA2 = 0.245;

double denssand = 2650.0;
double sgsand = denssand / denswater;
double delta = sgsand-1;
double ampbeds = ampbeds_factor * D50;
double epsilonp = 0.4;
double repose = -30.0;
double m = 4.0;
double alpha = m / (delta * g);
double be = 1.5;
double l1 = 1.9 * D50;
double l2 = 1.73;
double F0 = 0.025;
double F0_dim = correction_NT * F0 * sqrt(g * delta / D50);
double meanstle = alpha_lag * D50;
double A2_geom = M_PI / 4;
double A3_geom = M_PI / 6;
double k2 = 0.7;

double alpha_2 = 3000.0;
double D_star  = D50 * cbrt(g * delta / (nu * nu));
double w_s = nu / D50 * (pow(pow(10.36, 2) + 1.049 * pow(D_star, 3),(1./2.)) - 10.36);
double u_star_cr = sqrt(thetacr * g * delta * D50);
double alpha_min_SK = 50;
double alpha_max_SK = 10000;

double alpha_min_S = 50.0;
double alpha_max_S = 400.0;
double theta_min_S = 0.75;
double theta_max_S = 1.20;
double H_ref = 0.1166;
int keepsgrowing= 1;

//int nt=(Npx*(Npz+1));
//int nf = (int(ceil(Npx/10+1)))*7+2;
//int nf2 = int(ceil(Npx/10+1));
}

int admin::o2(int i_ex){
	/*adres vertaling voor periodieke rvw ten behoeve van dz*/
	int uit;
	if(i_ex>=0&&i_ex<Npx)uit=i_ex;
	else{
		if(i_ex==-1)uit=Npx-1;
		else{
			if(i_ex==Npx)uit=0;
			else {
				cerr<<"  ERROR: o2() buffer overflow (i_in="<<i_ex<<")"<<endl;
				cout<<"T="<<tijd<<" - ERROR: o2() buffer overflow (i_in="<<i_ex<<")"<<endl;
                uit=0;
			}
		}
	}
	return uit;
}
