// fft.cpp

#include "fft.h"

using namespace std;

fft::fft(int Npx_in){
	Npx=Npx_in;
	Npxf=Npx/2+1;
	f=(fftw_complex*)fftw_malloc(sizeof(fftw_complex)*Npxf);
	u=(double*)fftw_malloc(sizeof(double)*Npx);
	heen_plan=fftw_plan_dft_r2c_1d(Npx,u,f,FFTW_MEASURE);
	terug_plan=fftw_plan_dft_c2r_1d(Npx,f,u,FFTW_MEASURE);
}
fft::~fft(){
	fftw_destroy_plan(heen_plan);
	fftw_destroy_plan(terug_plan);
	fftw_free(u);fftw_free(f);
}

void fft::heen(vector<double> u_in,vector<complex<double> > &f_uit){
	for(int i=0;i<Npx;i++)u[i]=u_in[i];
	fftw_execute(heen_plan);
	for(int i=0;i<Npxf;i++){
		f_uit[i]=complex<double>(f[i][0],f[i][1])/(0.5*Npx);
	}
}

void fft::terug(vector<complex<double> >f_in,vector<double> &u_uit){
	for(int i=0;i<Npxf;i++){
		f[i][0]=f_in[i].real();
		f[i][1]=f_in[i].imag();
	}
	fftw_execute(terug_plan);
	for(int i=0;i<Npx;i++)u_uit[i]=u[i]*(0.5);
}
