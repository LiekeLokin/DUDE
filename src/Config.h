/*
 * config.h
 *
 *  Created on: Aug 10, 2021
 *      Author: jw
 */

#ifndef CONFIG_H_
#define CONFIG_H_

#include <string>

class Config {
	std::string line;
	unsigned int lineno = 0u;
	template <typename T>
	T scale(T val, const std::string& what) const;
	void warn(const std::string& message) const;
public:
	const bool DebugOutput;

	const int Npx;
	const int Npz;
	const double dtr;
	const double dt_write;
	const double tend;
	const double ampbeds_factor;
	const bool AllowFlowSep;
	const bool AllowAvalanching;
	const int SimpleLength;
	const int SimpleLengthFactor;
	const int numStab;
	const int Hifactor;
	const double Hcrit_global;
	const int transport_eq;
	const int alpha_varies;
	const double alpha_lag;
	const bool moeilijkdoen;
	const double correction_NT;
	const int Npsl_min;
	const int stle_factor;
	const bool write_velocities;

	const double q_in1;
	const double H0;
	const double ii;
	const double D50;
	const double thetacr;
	const double dts;
	const int nd;
	const std::string readbed;
	const std::string readfw;

	const double sepcritangle;
	const double g;
	const double kappa;
	const double tt;
	const double thresh;
	const int max_it;

	const double denswater;
	const double nu;
	const double BETA1;
	const double BETA2;

	const double denssand;
	const double epsilonp;
	const double repose;
	const double m;
	const double be;
	const double F0;
	const double A2_geom;
	const double A3_geom;
	const double k2;

	const double alpha_2;
	const double alpha_min_SK;
	const double alpha_max_SK;

	const double alpha_min_S;
	const double alpha_max_S;
	const double theta_min_S;
	const double theta_max_S;
	const double H_ref;
	const bool keepsgrowing;

	Config(const std::string& path = "config.cfg");
	virtual ~Config() {}
};

#endif /* CONFIG_H_ */
