// fft.h

#ifndef _FFT_H
#define _FFT_H

#include <complex>
#include <vector>
#include <fftw3.h>
#include "linalg.h"
using namespace std;

class fft{
private:
	fftw_plan heen_plan;
	fftw_plan terug_plan;
	fftw_complex *f;
	double *u;
	int Npxf;
	int Npx;
public:
	fft(int Npx_in);
	~fft();
	void heen(vector<double> u_in,vector<complex<double> > &f_uit);
	void terug(vector< complex<double> >f_in,vector<double> &u_uit);
};

#endif
