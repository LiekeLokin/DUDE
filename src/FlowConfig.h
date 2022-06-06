/*
 * FlowConfig.h
 *
 *  Created on: 24 nov. 2021
 *      Author: impor
 */

#ifndef SRC_FLOWCONFIG_H_
#define SRC_FLOWCONFIG_H_

class Config;

class FlowConfig {
public:
	const int Npx;		// number of grid points in horizontal x direction
	const int Npz;		// number of grid points in vertical z direction
	const double g;		// acceleration of gravity
	const double F;		// Forcing term in the momentum equation
	const double tresh;	// threshold accuracy for flow solver
	const int max_it;	// maximum number of iterations for flow solver

	FlowConfig(const Config& cfg);
	FlowConfig(const FlowConfig&) = default;
private:
	FlowConfig() = delete;
};

#endif /* SRC_FLOWCONFIG_H_ */
