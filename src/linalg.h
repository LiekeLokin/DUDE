//linalg.h

#ifndef _LINALG_H
#define _LINALG_H
#include <vector>
#include <map>
#include <iostream>
#include <cmath>

#define mat std::vector<std::vector<double> >
#define vec std::vector<double>

std::ostream& operator<<(std::ostream&, const vec&);
class crLU{
	public:
		crLU(const mat &fm,int s);
		crLU(const crLU &orig);
		crLU(int Llvc_ex,int Ulvc_ex,double Lval_ex[],double Uval_ex[],int Lci_ex[],int Uci_ex[],int Lri_ex[],int Uri_ex[],double diag_ex[]);
		~crLU();
		void bf(vec &v);
	private:
		int Llvc;
		int Ulvc;
		int lri;
		double *Lval;
		double *Uval;
		double *diag;
		int *Lci;
		int *Uci;
		int *Lri;
		int *Uri;
};

class spMat{
	public:
		spMat(int s,int blockWidth);
		double& e(int j,int i);
		double& e(int j,int b,int i);
		void softEmpty();
		~spMat();
		void rmrow(int rj);
		crLU LU();//destroys original matrix contents.
		std::map<int,double* > *r;
		int bw;
		int nrows;
};
class crMat{
	public:
		crMat(const mat &fm,int s);
		crMat(const spMat &fm,int s);
		~crMat();
		vec matvec(const vec &v);
		vec colmatvec(const mat &v,int s,int c);
		friend std::ostream& operator<<(std::ostream&, const crMat&);
	private:
		int lvc;
		int lri;
		double *val;
		int *ci;
		int *ri;
};

#if 0
vec matvec(const mat &m,const vec &v);
vec colmatvec(const mat &m,const mat &v,int s,int c);
double L2err(const vec &v);
#endif
double L2(const vec &v);
#if 0
mat LUc(const mat &m);
void LU(mat &m);
void bf(const mat &m,vec &v);
double dot(const vec &v1,const vec &v2);
double coldot(const vec &v1,const mat &v2,int s,int c);
double gmres(vec x0,vec &b,mat &A,int m);
double gmres(vec x0,vec &b,mat &A,mat &P,int m);
double gmres(vec x0,vec &b,crMat &A,int m);
#endif
double gmres(vec x0,vec &b,crMat &A,crLU &P,int m);
double interpolate(vec x,vec y, double x0);

#endif
