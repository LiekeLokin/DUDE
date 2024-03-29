/*
 * config.cpp
 *
 *  Created on: Aug 10, 2021
 *      Author: jw
 */

#include "Config.h"

#include <fstream>
#include <iostream>
#include <cassert>
#include <cmath>
#include <boost/algorithm/string.hpp>
#include <regex>

namespace {
std::regex group_expr{"^\\[(.+)]"};
std::regex number_param_expr{"^(\\S+)\\s*=\\s*([-+]?[0-9]*\\.?[0-9]+([eE][-+]?[0-9]+)?)\\s*(\\[([^\\]]*)\\])?"};
std::regex string_param_expr{"^(\\S+)\\s*=\\s*([^#]+)#?"};

bool equal(std::string a, std::string b) {
	std::transform(a.begin(), a.end(), a.begin(), ::tolower);
	std::transform(b.begin(), b.end(), b.begin(), ::tolower);
	return a == b;
}
} // anonymous

void Config::warn(const std::string& message) const {
	std::cerr << " ** " << message << " at line " << lineno << ": " << line << std::endl;

}

template <typename T>
T Config::scale(T val, const std::string& what) const {
	if (what.empty())
		return val;
	if (what == "s" || what == "sec" || what == "m" || what == "-" ||
			what == "m/m" || what == "m2/s" || what == "m/s2" ||
			what == "kg/m3")
		return val;
	if (what == "day")
		return val * 24 * 3600.0;
	if (what == "hr")
		return val * 3600.0;
	if (what == "min")
		return val * 60.0;
	if (what == "deg")
		return val * M_PI / 180.0;
	if (what == "mm")
		return val * 1e-3;
	warn("Unknown unit");
	return val;
}

template <typename T>
T getValue(const std::string& name, const std::string& val) {
	const char* env = getenv(name.c_str());
	std::stringstream ss(env ? env : val.c_str());
	std::cout << (env ? "ENV " : "    ");
	T value;
	ss >> value;
	return value;
}

#define ASSIGNDOUBLE(param) \
	if (equal(what[1], #param)) { \
		auto value = getValue<double>(what[1], what[2]); \
		param = scale(value, what[5]); \
		std::cout << what[1] << " = " << value << " [" << what[5] << "]" << std::endl; \
		continue; \
	}

#define ASSIGNINT(param) \
	if (equal(what[1], #param)) { \
		param = getValue<int>(what[1], what[2]); \
		std::cout << what[1] << " = " << param << std::endl; \
		continue; \
	}

#define ASSIGNSTRING(param) \
	if (equal(what[1], #param)) { \
		param = getValue<std::string>(what[1], what[2]); \
		std::cout << what[1] << " = '" << param << "'" << std::endl; \
		continue; \
	}

#define ASSIGNBOOL(param) \
	if (equal(what[1], #param)) { \
		auto value = getValue<char>(what[1], what[2]); \
		param = std::string("1TtYy").find(value) != std::string::npos; \
		std::cout << what[1] << " = " << std::boolalpha << param << std::endl; \
		continue; \
	}


Config::Config(const std::string& path) {
	std::ifstream in(path);
	assert(in && "config file not found");
	while (in) {
		char s[255];
		in.getline(s, sizeof s);
		lineno++;
		line = s;
		boost::trim(line);
		if (line.empty() || line.find("#") == 0)
			continue;

		std::cout << line << std::endl;

		std::smatch what;
		if (std::regex_search(line, what, group_expr)) {
			std::cout << " ** Group " << what[1] << " found" << std::endl;
			continue;
		}
		if (std::regex_search(line, what, number_param_expr)) {
			//std::cout << " ++ Param " << what[1] << " = " << what[2] << " [" << what[5] << "]" << std::endl;
			//for (auto w : what) std::cout << w << std::endl;
			ASSIGNBOOL(DebugOutput);
			ASSIGNINT(Npx);
			ASSIGNINT(Npz);
			ASSIGNDOUBLE(dtr);
			ASSIGNDOUBLE(dt_write);
			ASSIGNDOUBLE(tend);
			ASSIGNDOUBLE(ampbeds_factor);
			ASSIGNBOOL(AllowFlowSep);
			ASSIGNBOOL(AllowAvalanching);
			ASSIGNINT(SimpleLength);
			ASSIGNDOUBLE(SimpleLengthFactor);
			ASSIGNINT(numStab);
			ASSIGNINT(Hifactor);
			ASSIGNINT(Minfactor);
			ASSIGNDOUBLE(Hcrit_global);
			ASSIGNINT(transport_eq);
			ASSIGNINT(alpha_varies);
			ASSIGNDOUBLE(alpha_lag);
			ASSIGNBOOL(moeilijkdoen);
			ASSIGNDOUBLE(correction_NT);
			ASSIGNINT(Npsl_min);
			ASSIGNINT(stle_factor);
			ASSIGNBOOL(write_velocities);

			ASSIGNDOUBLE(q_in1);
			ASSIGNDOUBLE(H0);
			ASSIGNDOUBLE(ii);
			ASSIGNDOUBLE(D50);
			ASSIGNDOUBLE(thetacr);
			ASSIGNDOUBLE(dts);
			ASSIGNINT(nd);

			ASSIGNDOUBLE(sepcritangle);
			ASSIGNDOUBLE(g);
			ASSIGNDOUBLE(kappa);
			ASSIGNDOUBLE(tt);
			ASSIGNDOUBLE(thresh);
			ASSIGNINT(max_it);

			ASSIGNDOUBLE(denswater);
			ASSIGNDOUBLE(nu);
			ASSIGNDOUBLE(BETA1);
			ASSIGNDOUBLE(BETA2);

			ASSIGNDOUBLE(denssand);
			ASSIGNDOUBLE(epsilonp);
			ASSIGNDOUBLE(repose);
			ASSIGNDOUBLE(m);
			ASSIGNDOUBLE(be);
			ASSIGNDOUBLE(F0);
			ASSIGNDOUBLE(A2_geom);
			ASSIGNDOUBLE(A3_geom);
			ASSIGNDOUBLE(k2);

			ASSIGNDOUBLE(alpha_2);
			ASSIGNDOUBLE(alpha_min_SK);
			ASSIGNDOUBLE(alpha_max_SK);

			ASSIGNDOUBLE(alpha_min_S);
			ASSIGNDOUBLE(alpha_max_S);
			ASSIGNDOUBLE(theta_min_S);
			ASSIGNDOUBLE(theta_max_S);
			ASSIGNDOUBLE(H_ref);
			//ASSIGNBOOL(keepsgrowing);

			// Fall through
			warn("Unknown number parameter");
			continue;
		} else if (std::regex_search(line, what, string_param_expr)) {
			ASSIGNSTRING(FileName);
			ASSIGNSTRING(FileLevel);
			ASSIGNSTRING(ConsoleLevel);
			ASSIGNSTRING(readbed);
			ASSIGNSTRING(readfw);

			// Fall through
			warn("Unknown string parameter");
			continue;
		}

		// Fall through
		warn("Can't parse");
	}
}

