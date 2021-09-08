/*
 * test_config.cpp
 *
 *  Created on: Aug 10, 2021
 *      Author: jw
 */

#include "iostream"
#include "Config.h"

#define REPORT(what) \
	std::cout << #what" = " << std::boolalpha << cfg.what << std::endl


int main(int argc, char* argv[]) {
	std::string filename = (argc == 1) ? "config.cfg" : argv[1];
	Config cfg(filename);
	REPORT(DebugOutput);
	REPORT(Npx);
	REPORT(Npz);
	REPORT(dtr);
	REPORT(dt_write);
	REPORT(tend);
	REPORT(ampbeds_factor);
	REPORT(AllowFlowSep);
	REPORT(AllowAvalanching);
	REPORT(SimpleLength);
	REPORT(SimpleLengthFactor);
	REPORT(numStab);
	REPORT(Hifactor);
	REPORT(Hcrit_global);
	REPORT(transport_eq);
	REPORT(alpha_varies);
	REPORT(alpha_lag);
	REPORT(moeilijkdoen);
	REPORT(correction_NT);
	REPORT(Npsl_min);
	REPORT(stle_factor);

	REPORT(q_in1);
	REPORT(H0);
	REPORT(D50);
	REPORT(thetacr);
	REPORT(dts);
	REPORT(nd);
	REPORT(readbed);
	REPORT(readfw);

	REPORT(sepcritangle);
	REPORT(g);
	REPORT(kappa);
	REPORT(tt);
	REPORT(thresh);
	REPORT(max_it);

	REPORT(denswater);
	REPORT(nu);
	REPORT(BETA1);
	REPORT(BETA2);

	REPORT(denssand);
	REPORT(epsilonp);
	REPORT(repose);
	REPORT(m);
	REPORT(be);
	REPORT(F0);
	REPORT(A2_geom);
	REPORT(A3_geom);
	REPORT(k2);

	REPORT(alpha_2);
	REPORT(alpha_min_SK);
	REPORT(alpha_max_SK);

	REPORT(alpha_min_S);
	REPORT(alpha_max_S);
	REPORT(theta_min_S);
	REPORT(theta_max_S);
	REPORT(H_ref);
	REPORT(keepsgrowing);

}

