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
std::regex param_expr{"^(\\S+)\\s*=\\s*([-+]?[0-9]*\\.?[0-9]+([eE][-+]?[0-9]+)?)\\s*(\\[([^\\]]*)\\])?"};

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
	if (what == "hr")
		return val * 3600.0;
	if (what == "min")
		return val * 60.0;
	if (what == "deg")
		return val * M_PI / 180.0;
	warn("Unknown unit");
	return val;
}

#define ASSIGN(param) \
	if (equal(what[1], #param)) { \
		param = scale(value, what[5]); \
		continue; \
	}

#define ASSIGNBOOL(param) \
	if (equal(what[1], #param)) { \
		param = value != 0; \
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
		if (std::regex_search(line, what, param_expr)) {
			//std::cout << " ++ Param " << what[1] << " = " << what[2] << " [" << what[5] << "]" << std::endl;
			//for (auto w : what) std::cout << w << std::endl;
			std::stringstream ss(what[2]);
			double value;
			ss >> value;
			//std::cout << "value=" << value << std::endl;
			ASSIGNBOOL(DebugOutput);
			ASSIGN(Npx);
			ASSIGN(Npz);
			ASSIGN(dtr);
			ASSIGN(dt_write);
			ASSIGN(tend);
			ASSIGN(ampbeds_factor);
			ASSIGNBOOL(AllowFlowSep);
			ASSIGNBOOL(AllowAvalanching);
			ASSIGN(SimpleLength);
			ASSIGN(SimpleLengthFactor);
			ASSIGN(numStab);
			ASSIGN(Hifactor);
			ASSIGN(Hcrit_global);
			ASSIGN(transport_eq);
			ASSIGN(alpha_varies);
			ASSIGN(alpha_lag);
			ASSIGNBOOL(moeilijkdoen);
			ASSIGN(correction_NT);
			ASSIGN(Npsl_min);
			ASSIGN(stle_factor);

			ASSIGN(q_in1);
			ASSIGN(H0);
			ASSIGN(ii);
			ASSIGN(D50);
			ASSIGN(thetacr);
			ASSIGN(dts);
			ASSIGN(nd);
			ASSIGNBOOL(readbed);
			ASSIGNBOOL(readfw);

			ASSIGN(sepcritangle);
			ASSIGN(g);
			ASSIGN(kappa);
			ASSIGN(tt);
			ASSIGN(thresh);
			ASSIGN(max_it);

			ASSIGN(denswater);
			ASSIGN(nu);
			ASSIGN(BETA1);
			ASSIGN(BETA2);

			ASSIGN(denssand);
			ASSIGN(epsilonp);
			ASSIGN(repose);
			ASSIGN(m);
			ASSIGN(be);
			ASSIGN(F0);
			ASSIGN(A2_geom);
			ASSIGN(A3_geom);
			ASSIGN(k2);

			ASSIGN(alpha_2);
			ASSIGN(alpha_min_SK);
			ASSIGN(alpha_max_SK);

			ASSIGN(alpha_min_S);
			ASSIGN(alpha_max_S);
			ASSIGN(theta_min_S);
			ASSIGN(theta_max_S);
			ASSIGN(H_ref);
			ASSIGNBOOL(keepsgrowing);

			// Fall through
			warn("Unknown parameter");
			continue;
		}

		// Fall through
		warn("Can't parse");
	}
}

