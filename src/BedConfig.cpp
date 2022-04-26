/*
 * BedConfig.cpp
 *
 *  Created on: 24 nov. 2021
 *      Author: impor
 */
#include "BedConfig.h"
#include "Config.h"

#include <cmath>

BedConfig::BedConfig(const Config& cfg) :
	DebugOutput(cfg.DebugOutput),
	Npx(cfg.Npx),
	AllowAvalanching(cfg.AllowAvalanching),
	transport_eq(cfg.transport_eq),
	alpha_varies(cfg.alpha_varies),
	alpha_lag(cfg.alpha_lag),
	moeilijkdoen(cfg.moeilijkdoen),
	Npsl_min(cfg.Npsl_min),
	stle_factor(cfg.stle_factor),
	D50(cfg.D50),
	thetacr(cfg.thetacr),
	nd(cfg.nd),
	sepcritangle(cfg.sepcritangle),
	g(cfg.g),
	tt(cfg.tt),
	delta(cfg.denssand / cfg.denswater - 1),
	epsilonp(cfg.epsilonp),
	repose(cfg.repose),
	alpha(cfg.m / (delta * g)),
	be(cfg.be),
	//l2(1 / tan(-repose)),
	l2(1.73),//TODO moet die oude zijn met de tan van de repose ering
	F0_dim(cfg.correction_NT * cfg.F0 * sqrt(g * delta / D50)),
	meanstle(alpha_lag * cfg.D50),
	alpha_2(cfg.alpha_2),
	D_star(D50 * cbrt(g * delta / (cfg.nu * cfg.nu))),
	w_s(cfg.nu / D50 * (pow(pow(10.36, 2) + 1.049 * pow(D_star, 3),(1./2.)) - 10.36)),
	u_star_cr(sqrt(thetacr * g * delta * D50)),
	alpha_min_SK(cfg.alpha_min_SK),
	alpha_max_SK(cfg.alpha_max_SK),
	alpha_min_S(cfg.alpha_min_S),
	alpha_max_S(cfg.alpha_max_S),
	theta_min_S(cfg.theta_min_S),
	theta_max_S(cfg.theta_max_S),
	H_ref(cfg.H_ref)
//	keepsgrowing(cfg.keepsgrowing)
{}
