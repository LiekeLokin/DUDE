/*
 * FlowConfig.cpp
 *
 *  Created on: 24 nov. 2021
 *      Author: impor
 */
#include "FlowConfig.h"
#include "Config.h"

FlowConfig::FlowConfig(const Config& cfg) :
	Npx(cfg.Npx),
	Npz(cfg.Npz),
	g(cfg.g),
	F(g * cfg.ii),
	tresh(cfg.thresh),
	max_it(cfg.max_it)
{}
