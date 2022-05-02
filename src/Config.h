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

class Config {
	std::ofstream* cfglog;
	std::string line;
	unsigned int lineno = 0u;
	template <typename T>
	T scale(T val, const std::string& what) const;
	void warn(const std::string& message) const;
	template <typename T>
	T getValue(const std::string& name, const std::string& val);
public:
	bool DebugOutput;
	std::string FileName;
	std::string FileLevel;
	std::string ConsoleLevel;

	int Npx;
	int Npz;
	double dtr;
	double dt_write;
	double tend;
	double ampbeds_factor;
	bool AllowFlowSep;
	bool AllowAvalanching;
	int SimpleLength;
	int SimpleLengthFactor;
	int numStab;
	int Hifactor;
	int Minfactor;
	double Hcrit_global;
	int transport_eq;
	int alpha_varies;
	double alpha_lag;
	bool moeilijkdoen;
	double correction_NT;
	int Npsl_min;
	int stle_factor;
	bool write_velocities;

	double q_in1;
	double H0;
	double ii;
	double D50;
	double thetacr;
	double dts;
	int nd;
	std::string readbed;
	std::string readfw;

	double sepcritangle;
	double g;
	double kappa;
	double tt;
	double thresh;
	int max_it;

	double denswater;
	double nu;
	double BETA1;
	double BETA2;

	double denssand;
	double epsilonp;
	double repose;
	double m;
	double be;
	double F0;
	double A2_geom;
	double A3_geom;
	double k2;

	double alpha_2;
	double alpha_min_SK;
	double alpha_max_SK;

	double alpha_min_S;
	double alpha_max_S;
	double theta_min_S;
	double theta_max_S;
	double H_ref;
//	bool keepsgrowing;
	bool Lrangefix;

	Config(const std::string& path = "config.cfg");
	virtual ~Config() {}
};

#endif /* CONFIG_H_ */
