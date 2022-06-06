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
	vec b;
	vec bp;
	vec flux;
	vec x;
	vec Sr;
	std::vector<int> fsz;
	int o3(int i_in) const;
	//void detQ(vec ub, vec &dhdx);
	void detQcr(const vec& ub, const vec& Umean, vec &dhdx);
	void sep_migr_lee(const vec& fluxtot, const vec& oldb);
	void sep_sort_fsz(int num);
	vec sep_tau_distr(vec ub);
	std::vector<int> setFSZ(int xsi, int nfsz, int wavelet);
	int paramSepline(int xsi, int xti, int xci, int nfsz);
	vec paramFindNeighbors(double x_p, int xi) const;
	vec crossPoint_migrlee(double xl, double xr, int max_it,int dir, double b1, double b2, double tol, int xi, int j) const;
	double area2D_Polygon(int n, const vec& xarr, const vec& yarr) const;
	int findTrough(int xsi, const vec& bed) const ;
	int findCrest(int xsi, const vec& bed) const;
	vec filter(int np, const vec& inp_arr) const;
	void smooth_param(int np, int num);
	vec fftBed(const vec& bed, int fftnum) const;
	std::vector<std::complex<double> > fftbot(const vec& bed) const;
	int maxloc_complex(const std::vector<std::complex<double> >& du) const;
	void avalanche(); //OLAV 2014 01 30
	double detAlphaLag(const vec& ub, int method,int suppressoutput) const;
	vec detDistributeFunc(double alpha_lag1,double deltax) const;

public:
	bottom() = delete;
	bottom(const BedConfig& cfg);
	~bottom();
	vec update(const vec& ub,const vec& Umean, vec &bss1, vec &fluxtot, vec &dhdx);
	vec update_flowsep(const vec& ub,const vec& Umean, vec &bss1, vec &bss2, vec &fluxtot, vec &dhdx);
	vec readBottomInp(const std::string& readbed);
	void writeBottom() const;
	void checkFlowsep();
	void write_flowsep() const;
	void setShape(const vec& b_in);
	void setSin(double amp,int n);
	void setRand(double amp);
	void setRand(double amp,int seed);
	void setWave(int xwi, int xcin);
	void setDist(double amp_dist);
	void setDistSin (double ampdist,int n);
	void setMidSin(double amp, double length);
	void setCustom(double amp, int n);
	vec get_dhdx() const;
	vec getShape(int sepflag) const;
	std::vector<int> getFsz() const;
	std::vector<double> getSr()const;
	double detint1(const vec& bed) const;
	double detint2(const vec& bed) const;
	vec detNd(const vec& bot) const;
	vec detNd_fft(const vec& bot, int fftnum) const;
	double detMigr(const vec& current, const vec& next) const;
	vec smooth(const vec& bed_in) const;
	double meanval(const vec& vr, int Np) const;
	double minval(const vec& vr, int Np) const;
};


#endif
