// fft.h

#ifndef _FFT_H
#define _FFT_H

#include <complex>
#include <vector>
#include <fftw3.h>
#include "linalg.h"

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
	void heen(std::vector<double> u_in, std::vector<std::complex<double> > &f_uit);
	void terug(std::vector<std::complex<double> >f_in, std::vector<double> &u_uit);
};

#endif
