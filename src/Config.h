/*
 * config.h
 *
 *  Created on: Aug 10, 2021
 *      Author: jw
 */

#ifndef CONFIG_H_
#define CONFIG_H_

#include <string>
#include <iosfwd>
#include <cmath>

class Config {
	std::ofstream* cfglog;
	const double D2R = M_PI / 180.0;
	std::string line;
	unsigned int lineno = 0u;
	template <typename T>
	T scale(T val, const std::string& what) const;
	void warn(const std::string& message) const;
	template <typename T>
	T getValue(const std::string& name, const std::string& val);
public:
	bool DebugOutput = false;
	std::string FileName = "dude.log";
	std::string FileLevel = "info";
	std::string ConsoleLevel = "debug";

	int Npx = 120;
	int Npz = 25;
	double dtr = 240;
	double dt_write = 240;
	double tend = 3 * 60 * 60;
	double ampbeds_factor = 100.1;
	double bedResetFac = 0.5;
	bool AllowFlowSep = false;
	bool AllowAvalanching = true;
	int SimpleLength = 0;
	int SimpleLengthFactor = 7;
	int numStab = 30;
	int Hifactor = 15;
	int Minfactor = 3;
	double Hcrit_global = 5;
	int transport_eq = 2;
	int alpha_varies = 3;
	double alpha_lag = 100;
	bool moeilijkdoen = false;
	double correction_NT = 1;
	int Npsl_min = 40;
	int stle_factor = 5;
	bool write_velocities = false;

	bool use_H_only = false;
	double q_in1 = 6.25;
	double H0 = 9.1;
	double Lini = 1;
	double ii = 1.1e-4;
	double D50 = 0.18e-3;
	double thetacr = 0.053;
	double dts = 0.01;
	int nd = 1;
	std::string readbed;
	std::string readfw;

	double sepcritangle = 10 * D2R;
	double g = 9.81;
	double kappa = -0.407;
	double tt = 100;
	double thresh = 1e-8;
	int max_it = 8;

	double denswater = 1000;
	double nu = 1e-6;
	double BETA1 = 0.245;
	double BETA2 = 0.245;
	bool S_Av_const = 0;

	double denssand = 2650;
	double epsilonp = 0.4;
	double repose = -30 * D2R;
	double m = 4;
	double be = 1.5;
	double F0 = 0.025;
	double A2_geom = 45 * D2R;
	double A3_geom = 30 * D2R;
	double k2 = 0.7;

	double alpha_2 = 3000;
	double alpha_min_SK = 50;
	double alpha_max_SK = 10000;

	double alpha_min_S = 50;
	double alpha_max_S = 400;
	double theta_min_S = 0.9;
	double theta_max_S = 1.35;
	double H_ref = 0.1166;
//	bool keepsgrowing;
	bool Lrangefix = false;

	Config(const std::string& path = "config.cfg");
	virtual ~Config() {}
};

#endif /* CONFIG_H_ */
