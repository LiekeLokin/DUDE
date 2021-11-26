// bottom.h

#ifndef _BOTTOM_H
#define _BOTTOM_H

#include "BedConfig.h"
#include "linalg.h"
#include <vector>
#include <complex>

class bottom{
private:
	const BedConfig cfg;
	const int Npx;
	const int nf;
	const int nf2;
	vec *b;
	vec *bp;
	vec *flux;
	vec *x;
	vec *Sr;
	std::vector<int> *fsz;
	int o3(int i_in) const;
	void detQ(vec ub, vec &dhdx);
	void detQcr(vec ub, vec &dhdx);
	void sep_migr_lee(vec flux, vec oldb);
	void sep_sort_fsz(int num);
	vec sep_tau_distr(vec ub);
	std::vector<int> setFSZ(int xsi, int nfsz, int wavelet);
	int paramSepline(int xsi, int xti, int xci, int nfsz);
	vec paramFindNeighbors(double x_p, int xi);
	vec crossPoint_migrlee(double xl, double xr, int max_it,int dir, double b1, double b2, double tol, int xi, int j);
	double area2D_Polygon(int n, vec xarr, vec yarr );
	int findTrough(int xsi, vec bed);
	int findCrest(int xsi, vec bed);
	vec filter(int np, vec inp_arr);
	void smooth_param(int np, int num);
	vec fftBed(vec bed, int fftnum);
	std::vector<std::complex<double> > fftbot(vec bed);
	int maxloc_complex(std::vector<std::complex<double> > du);
	void avalanche(); //OLAV 2014 01 30
	double detAlphaLag(vec ub, int method,int suppressoutput);
	vec detDistributeFunc(double alpha_lag1,double deltax);

public:
	bottom() = delete;
	bottom(const BedConfig& cfg);
	~bottom();
	vec update(vec ub, vec &bss1, vec &fluxtot, vec &dhdx);
	vec update_flowsep(vec ub, vec &bss1, vec &bss2, vec &fluxtot, vec &dhdx);
	vec readBottomInp(const std::string& readbed);
	void writeBottom();
	void checkFlowsep();
	void write_flowsep();
	void setShape(vec);
	void setSin(double amp);
	void setSin(double amp,int n);
#if 0
	void setRand(double amp);
	void setRand(double amp,int seed);
#endif
	void setWave(int xwi, int xcin);
#if 0
	void setDist(double amp_dist);
	void setDistSin (double ampdist,int n);
#endif
	void setMidSin(double amp, double length);
	void setCustom(double amp, int n);
	vec getShape(int sepflag);
	std::vector<int> getFsz();
	std::vector<double> getSr();
	double detint1(vec bed);
	double detint2(vec bed);
	vec detNd(vec bot);
	vec detNd_fft(vec bot, int fftnum);
	double detMigr(vec current,vec next);
	vec smooth(vec bed_in);
	double meanval(vec vr, int Np);
	double minval(vec vr, int Np);
};


#endif
