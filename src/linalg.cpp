// linalg.cpp

#include "linalg.h"
#include <iostream>
#include <cmath>

using namespace std;

double interpolate(vec x,vec y,double x0){
	int in_range = -999;
	double y0;
	
	for(int j=0;j<x.size();j++){
		if(x[j]>=x0) {
			in_range =j;
			break;
		}
	}

	if(in_range>0){ //x0 is between two values in the range
		y0=y[in_range-1]+(y[in_range]-y[in_range-1])*(x0-x[in_range-1])/(x[in_range]-x[in_range-1]);
	}
	else if(in_range==0){ //the first x is already larger than x0
		y0=y[0];
	}
	else if(in_range==-999){ //the last x is still smaller than x0
		y0=y[x.size()-1];
	}
	
	return y0;
}

#if 0
vec matvec(const mat &m,const vec &v){
	vec p(v.size());
	for(int j=0;j<v.size();j++){
		p[j]=0.0;
		for(int i=0;i<v.size();i++){
			p[j]+=m[j][i]*v[i];
		}
	}
	return p;
}
vec colmatvec(const mat &m,const mat &v,int s,int c){
	vec p(s);
	for(int j=0;j<s;j++){
		p[j]=0.0;
		for(int i=0;i<s;i++){
			p[j]+=m[j][i]*v[i][c];
		}
	}
	return p;
}
#endif

double L2(const vec &v){
	double sum=v[0]*v[0];
	for(int i=1;i<v.size();i++){
		sum+=v[i]*v[i];
	}
	sum/=(double)v.size();
	sum=sqrt(sum);
	return sum;
}

double L2err(const vec &v){
	double sum=v[0]*v[0];
	for(int i=1;i<v.size();i++){
		sum+=v[i]*v[i];
	}
	//sum/=(double)v.size();
	sum=sqrt(sum);
	return sum;
}

#if 0
void LU(mat &m){
	int n=m.size();
	int nm=m.size()-1;
	for(int k=0;k<nm;k++){
		if(m[k][k]==0.0)cerr<<"pivot zero on row "<<k<<endl;
		for(int j=k+1;j<n;j++){
			if(m[j][k]!=0.0){
				double f=m[j][k]/m[k][k];
				m[j][k]=f;
				for(int i=k+1;i<n;i++){
					m[j][i]-=f*m[k][i];
				}
			}
		}
	}
}
mat LUc(const mat &m){
	mat mc=m;
	LU(mc);
	return mc;
}

void bf(const mat &m, vec &v){
	int n=v.size();
	for(int i=1;i<n;i++)
		for(int j=0;j<i;j++)v[i]-=m[i][j]*v[j];
	for(int i=n-1;i>=0;i--){
		for(int j=i+1;j<n;j++)v[i]-=m[i][j]*v[j];
		v[i]/=m[i][i];
	}
}
double dot(const vec &v1,const vec &v2){
	int n=v1.size();
	double d=0.0;
	for(int i=0;i<n;i++)d+=v1[i]*v2[i];
	return d;
}
#endif

double coldot(const vec &v1,const mat &v2,int s,int c){
	int n=s;
	double d=0.0;
	for(int i=0;i<n;i++)d+=v1[i]*v2[i][c];
	return d;
}

#if 0
double gmres(vec x0,vec &b,mat &A,int m){
	crMat spA(A,x0.size());
	return gmres(x0,b,spA,m);
}
double gmres(vec x0,vec &b,mat &A,mat &P,int m){
	crMat spA(A,x0.size());
	crLU spP(P,x0.size());
	return gmres(x0,b,spA,spP,m);
}

double gmres(vec x0,vec &b,crMat &A,int m){
	int n=b.size();
	// JW mat H(m+1,m+1);
	// JW mat V(n,m+1);
    mat H(m+1,vec(m+1));
	mat V(n,vec(m+1));
	double beta;
	vec betae1(m+1,0.0);
	vec r0(n,0.0);
	vec wj(n,0.0);
	vec mv(n,0.0);
	mv=A.matvec(x0);
	for(int i=0;i<n;i++)r0[i]=b[i]-mv[i];
	beta=L2err(r0);
	for(int i=0;i<n;i++)V[i][0]=r0[i]/beta;

	for(int j=0;j<m;j++){
		wj=A.colmatvec(V,n,j);
		for(int i=0;i<=j;i++){
			H[i][j]=coldot(wj,V,n,i);
			for(int il=0;il<n;il++)wj[il]-=H[i][j]*V[il][i];
		}
		H[j+1][j]=L2err(wj);
		if(fabs(H[j+1][j])<1.e-20){
			cerr<<"Lucky breakdown!!!"<<endl;
			break;
		}
		for(int il=0;il<n;il++)V[il][j+1]=wj[il]/H[j+1][j];
	}
	//givens
	betae1[0]=beta;
	vec row0(m,0.0);
	vec row1(m,0.0);
	for(int j=0;j<m;j++){
		double c=H[j][j]/sqrt(H[j][j]*H[j][j]+H[j+1][j]*H[j+1][j]);
		double s=-H[j+1][j]/sqrt(H[j][j]*H[j][j]+H[j+1][j]*H[j+1][j]);
		for(int i=0;i<m;i++){row0[i]=0.0;row1[i]=0.0;}
		for(int i=j;i<m;i++){
			row0[i]=c*H[j][i]-s*H[j+1][i];
			row1[i]=s*H[j][i]+c*H[j+1][i];
		}
		double val[2];
		val[0]=c*betae1[j]-s*betae1[j+1];
		val[1]=s*betae1[j]+c*betae1[j+1];
		betae1[j]=val[0];
		betae1[j+1]=val[1];
		for(int i=0;i<m;i++){
			H[j][i]=row0[i];
			H[j+1][i]=row1[i];
		}
	}
	for(int i=m-1;i>=0;i--){
		for(int j=i+1;j<m;j++)betae1[i]-=H[i][j]*betae1[j];
		betae1[i]/=H[i][i];
	}
	for(int j=0;j<n;j++){
		for(int i=0;i<m;i++)x0[j]+=V[j][i]*betae1[i];
	}
	mv=A.matvec(x0);
	for(int i=0;i<n;i++)r0[i]=b[i]-mv[i];
	for(int i=0;i<n;i++)b[i]=x0[i];
	return L2err(r0);
}
#endif

double gmres(vec x0, vec &b, crMat &A, crLU &P, int m) {
	auto n = b.size();
	// JW mat H(m+1,m+1);
	// JW mat V(n,m+1);
	mat H(m + 1, vec(m + 1));
	mat V(n, vec(m + 1));
	vec betae1(m + 1, 0.0);
	vec r0(n, 0.0);
	//vec wj(n, 0.0);
	auto mv = A.matvec(x0);
	for (auto i = 0; i < n; i++)
		r0[i] = b[i] - mv[i];
	P.bf(r0);
	auto beta = L2err(r0);
	for (auto i = 0; i < n; i++)
		V[i][0] = r0[i] / beta;

	for (auto j = 0; j < m; j++) {
		auto wj = A.colmatvec(V, n, j);
		P.bf(wj);
		for (auto i = 0; i <= j; i++) {
			H[i][j] = coldot(wj, V, n, i);
			for (auto il = 0; il < n; il++)
				wj[il] -= H[i][j] * V[il][i];
		}
		H[j + 1][j] = L2err(wj);
		if (fabs(H[j + 1][j]) < 1.e-20) {
			cerr << "Lucky breakdown!!!" << endl;
			break;
		}
		for (auto il = 0; il < n; il++)
			V[il][j + 1] = wj[il] / H[j + 1][j];
	}
	//givens
	betae1[0] = beta;
	vec row0(m, 0.0);
	vec row1(m, 0.0);
	for (auto j = 0; j < m; j++) {
		auto c = H[j][j]
				/ sqrt(H[j][j] * H[j][j] + H[j + 1][j] * H[j + 1][j]);
		auto s = -H[j + 1][j]
				/ sqrt(H[j][j] * H[j][j] + H[j + 1][j] * H[j + 1][j]);
		for (auto i = 0; i < m; i++) {
			row0[i] = 0.0;
			row1[i] = 0.0;
		}
		for (auto i = j; i < m; i++) {
			row0[i] = c * H[j][i] - s * H[j + 1][i];
			row1[i] = s * H[j][i] + c * H[j + 1][i];
		}
		double val[2];
		val[0] = c * betae1[j] - s * betae1[j + 1];
		val[1] = s * betae1[j] + c * betae1[j + 1];
		betae1[j] = val[0];
		betae1[j + 1] = val[1];
		for (auto i = 0; i < m; i++) {
			H[j][i] = row0[i];
			H[j + 1][i] = row1[i];
		}
	}
	for (auto i = m - 1; i >= 0; i--) {
		for (auto j = i + 1; j < m; j++)
			betae1[i] -= H[i][j] * betae1[j];
		betae1[i] /= H[i][i];
	}
	for (auto j = 0; j < n; j++) {
		for (auto i = 0; i < m; i++)
			x0[j] += V[j][i] * betae1[i];
	}
	mv = A.matvec(x0);
	for (auto i = 0; i < n; i++)
		r0[i] = b[i] - mv[i];
	for (auto i = 0; i < n; i++)
		b[i] = x0[i];
	return L2err(r0);
}

ostream & operator<<(ostream & s,const vec &v){
	for(int j=0;j<v.size();j++){
		s<<v[j];
	}
	s<<endl;
	return s;
}
ostream & operator<<(ostream & s,const crMat & M){
	for(int j=0;j<M.lri-1;j++){
		for(int i=M.ri[j];i<M.ri[j+1];i++){
			s<<j+1<<" "<<M.ci[i]+1<<" "<<M.val[i]<<endl;
		}
	}
	return s;
}

crMat::crMat(const spMat &fm,int s){
	lri=s+1;
	lvc=0;
	for(int j=0;j<s;j++){
		for(map<int,double*>::const_iterator p=fm.r[j].begin();
				p!=fm.r[j].end();p++){
			for(int i=0;i<fm.bw;i++){
				if(fabs((p->second)[i])>1.e-16)lvc++;
			}
		}
	}
	val=new double[lvc];
	ci=new int[lvc];
	ri=new int[lri];
	int rowi=0;
	for(int j=0;j<s;j++){
		ri[j]=rowi;
		for(map<int,double*>::const_iterator p=fm.r[j].begin();
				p!=fm.r[j].end();p++){
			for(int i=0;i<fm.bw;i++){
				if(fabs((p->second)[i])>1.e-16){
					val[rowi]=p->second[i];
					ci[rowi]=p->first*fm.bw+i;
					rowi++;
				}
			}
		}
	}
	ri[s]=lvc;
	/*cout<<"val"<<endl;
	for(int i=0;i<lvc;i++)cout<<val[i]<<" ";
	cout<<endl<<"ci"<<endl;
	for(int i=0;i<lvc;i++)cout<<ci[i]<<" ";
	cout<<endl<<"ri"<<endl;
	for(int i=0;i<s;i++)cout<<ri[i]<<" ";
	cout<<endl;*/
}

crMat::crMat(const mat &fm,int s){
	lri=s+1;
	lvc=0;
	for(int j=0;j<s;j++){
		for(int i=0;i<s;i++){
			if(fabs(fm[j][i])>1.e-16){
				lvc++;
			}
		}
	}
	val=new double[lvc];
	ci=new int[lvc];
	ri=new int[lri];
	int rowi=0;
	for(int j=0;j<s;j++){
		ri[j]=rowi;
		for(int i=0;i<s;i++){
			if(fabs(fm[j][i])>1.e-16){
				val[rowi]=fm[j][i];
				ci[rowi]=i;
				rowi++;
			}
		}
	}
	ri[s]=lvc;
	/*cout<<"val"<<endl;
	for(int i=0;i<lvc;i++)cout<<val[i]<<" ";
	cout<<endl<<"ci"<<endl;
	for(int i=0;i<lvc;i++)cout<<ci[i]<<" ";
	cout<<endl<<"ri"<<endl;
	for(int i=0;i<s;i++)cout<<ri[i]<<" ";
	cout<<endl;*/
}

crMat::~crMat(){
	delete[] val;
	delete[] ci;
	delete[] ri;
}
vec	crMat::matvec(const vec &v){
	vec p(v.size());
	for(int j=0;j<v.size();j++){
		p[j]=0.0;
		for(int i=ri[j];i<ri[j+1];i++){
			p[j]+=val[i]*v[ci[i]];
		}
	}
	return p;
}
vec	crMat::colmatvec(const mat &v,int s,int c){
	vec p(s);
	for(int j=0;j<s;j++){
		p[j]=0.0;
		for(int i=ri[j];i<ri[j+1];i++){
			p[j]+=val[i]*v[ci[i]][c];
		}
	}
	return p;
}

void crLU::bf(vec &v){
	int n=v.size();
	for(int j=1;j<n;j++){
		for(int i=Lri[j];i<Lri[j+1];i++) {
			if(i>=Llvc)cerr<<"zoveel elementen heb ik niet"<<endl;
			v[j]-=Lval[i]*v[Lci[i]];
		}
	}
	for(int j=n-1;j>=0;j--){
		for(int i=Uri[j];i<Uri[j+1];i++){
			if(i>=Ulvc)cerr<<"zoveel elementen heb ik niet"<<endl;
			v[j]-=Uval[i]*v[Uci[i]];
		}
		v[j]/=diag[j];
	}
}

crLU::crLU(int Llvc_ex,int Ulvc_ex,double Lval_ex[],double Uval_ex[],int Lci_ex[],int Uci_ex[],int Lri_ex[],int Uri_ex[],double diag_ex[]){
	Llvc=Llvc_ex;
	Ulvc=Ulvc_ex;
	Lci=Lci_ex;
	Uci=Uci_ex;
	Lval=Lval_ex;
	Uval=Uval_ex;
	Lri=Lri_ex;
	Uri=Uri_ex;
	diag=diag_ex;
	/*
	cout<<"Llvc: "<<Llvc<<" Ulvc: "<<Ulvc<<endl;
	cout<<"crLU up"<<endl;
	cout<<"Uval[i]: "; for(int i=0;i<Uri[9];i++)cout<<Uval[i]<<" "; cout<<endl;
	cout<<"Uci[i]:  "; for(int i=0;i<Ulvc;i++)cout<<Uci[i]<<" ";	cout<<endl;
	cout<<"diag"<<endl;
	cout<<"diag[i]: "; for(int i=0;i<9;i++)cout<< diag[i]<<" "; cout<<endl;
	cout<<"lo"<<endl;
	cout<<"Lval[i]: "; for(int i=0;i<Lri[9];i++) cout<<Lval[i]<<" "; cout<<endl;
	cout<<"Lci[i]:  "; for(int i=0;i<Llvc;i++) cout<<Lci[i]<<" "; cout<<endl;
	cout<<"ri"<<endl;
	cout<<"Uri[i]:  "; for(int i=0;i<=1100;i++)cout<<Uri[i]<<" "; cout<<endl;
	cout<<"Lri[i]:  "; for(int i=0;i<=1100;i++)cout<<Lri[i]<<" "; cout<<endl<<endl;
	*/
}

crLU::crLU(const crLU &orig){
	cerr<<"kom ik hier?"<<endl;
	lri=orig.lri;
	Llvc=orig.Llvc;
	Ulvc=orig.Ulvc;
	Lval=new double[Llvc];
	Uval=new double[Ulvc];
	diag=new double[lri];
	Lci=new int[Llvc];
	Uci=new int[Ulvc];
	Lri=new int[lri];
	Uri=new int[lri];
	for(int i=0;i<Llvc;i++){
		Lval[i]=orig.Lval[i];
		Lci[i]=orig.Lci[i];
	}
	for(int i=0;i<lri;i++)Lri[i]=orig.Lri[i];
	for(int i=0;i<Ulvc;i++){
		Uval[i]=orig.Uval[i];
		Uci[i]=orig.Uci[i];
	}
	for(int i=0;i<lri;i++)Uri[i]=orig.Uri[i];
}

crLU::crLU(const mat &fm,int s){
	cerr<<"kom ik hier?"<<endl;
	lri=s+1;
	Llvc=0;
	Ulvc=0;
	for(int j=0;j<s;j++){
		for(int i=0;i<j;i++){
			if(fabs(fm[j][i])>1.e-16)Llvc++;
		}
		for(int i=j+1;i<s;i++){
			if(fabs(fm[j][i])>1.e-16)Ulvc++;
		}
	}
	Lval=new double[Llvc];
	Uval=new double[Ulvc];
	diag=new double[lri];
	Lci=new int[Llvc];
	Uci=new int[Ulvc];
	Lri=new int[lri];
	Uri=new int[lri];
	int Lrowi=0;
	int Urowi=0;
	for(int j=0;j<s;j++){
		Lri[j]=Lrowi;
		Uri[j]=Urowi;
		for(int i=0;i<j;i++){
			if(fabs(fm[j][i])>1.e-16){
				Lval[Lrowi]=fm[j][i];
				Lci[Lrowi]=i;
				Lrowi++;
			}
		}
		diag[j]=fm[j][j];
		for(int i=j+1;i<s;i++){
			if(fabs(fm[j][i])>1.e-16){
				Uval[Urowi]=fm[j][i];
				Uci[Urowi]=i;
				Urowi++;
			}
		}
	}
	Lri[s]=Llvc;
	Uri[s]=Ulvc;
}

crLU::~crLU(){
	delete[] Lval;
	delete[] Uval;
	delete[] diag;
	delete[] Lci;
	delete[] Uci;
	delete[] Lri;
	delete[] Uri;
}

spMat::spMat(int s,int blockWidth){
	nrows=s;
	bw=blockWidth;
	r=new map<int,double*>[nrows];
}
spMat::~spMat(){
	for(int j=0;j<nrows;j++){
		for(map<int,double*>::iterator p=(r[j].begin());
				p!=r[j].end();p++){
			delete[] p->second;
		}
	}
	delete[] r;
}
double& spMat::e(int j,int i){
	int ni=i%bw;
	int b=i/bw;
	return e(j,b,ni);
}
double& spMat::e(int j,int b,int i){
	if(NULL==(r[j])[b]){
		(r[j])[b]=new double[bw];
		for(int ii=0;ii<bw;ii++)(r[j])[b][ii]=0.0;
	}
	return ((r[j])[b])[i];
}

void spMat::softEmpty(){
	for(int j=0;j<nrows;j++){
		for(map<int,double*>::iterator p=(r[j].begin());
				p!=r[j].end();p++){
			for(int i=0;i<bw;i++)(p->second)[i]=0.0;
		}
	}
}
void spMat::rmrow(int rj){
	for(map<int,double*>::iterator p=(r[rj].begin());
				p!=r[rj].end();p++){
		delete[] p->second;
	}
	r[rj].clear();
}

crLU spMat::LU(){
	int n=nrows;
	int nm=nrows-1;
	int Llvc=0;
	int Ulvc=0;
	spMat lo(nrows,bw);
	double *diag=new double[nrows];
	for(int i=0;i<nrows;i++)diag[i]=0.0;
	//cerr<<"construct LU"<<endl;
	for(int k=0;k<nm;k++){
		double pivot=e(k,k);
		diag[k]=pivot;
		if(pivot==0.0)cerr<<"ERROR: pivot zero on row "<<k<<endl;
		for(int j=k+1;j<n;j++){
			if(r[j].begin()->first*bw<=k && (r[j].begin()->second)[k%bw]!=0.0){
				double f=(r[j].begin()->second)[k%bw]/pivot;
				lo.e(j,k)=f;
				//if(fabs(f)<=1.e-16)cerr<<"Betrapt!"<<endl;
				if(fabs(f)>1.e-16)Llvc++;
				for(int i=k%bw+1;i<bw;i++)
					(r[j].begin()->second)[i]-=f*(r[k])[k/bw][i];
				for(map<int,double*>::const_iterator p=(++r[k].begin());
						p!=r[k].end();p++){
					for(int i=0;i<bw;i++){
						e(j,p->first,i)-=f*(p->second)[i];
					}
				}
			}
			if(k%bw==bw-1){
				delete[] r[j][k/bw];
				r[j].erase(k/bw);
			}
		}
	}
	//cerr<<"copy LU to friendly structure"<<endl;
	diag[nm]=e(nm,nm);
	double* Lval=new double[Llvc];
	int* Lci=new int[Llvc];
	int* Lri=new int[nrows+1];
	for(int j=0;j<n;j++){
		for(int i=j%bw+1;i<bw;i++){
			if(fabs((r[j].begin()->second)[i])>1.e-16)Ulvc++;
		}
		for(map<int,double*>::const_iterator p=(++r[j].begin());
				p!=r[j].end();p++){
			for(int i=0;i<bw;i++){
				if(fabs((p->second)[i])>1.e-16)Ulvc++;
			}
		}
	}
	double* Uval=new double[Ulvc];
	int* Uci=new int[Ulvc];
	int* Uri=new int[nrows+1];
	int Lrowi=0;
	int Urowi=0;
	for(int j=0;j<n;j++){
		Lri[j]=Lrowi;
		Uri[j]=Urowi;
		for(map<int,double*>::const_iterator p=lo.r[j].begin();
				p!=lo.r[j].end();p++){
			for(int i=0;i<bw;i++){
				if(fabs((p->second)[i])>1.e-16){
					Lval[Lrowi]=(p->second)[i];
					Lci[Lrowi]=p->first*bw+i;
					Lrowi++;
				}
			}
		}
		for(int i=j%bw+1;i<bw;i++){
			if(fabs((r[j].begin()->second)[i])>1.e-16){
				Uval[Urowi]=(r[j].begin()->second)[i];
				Uci[Urowi]=r[j].begin()->first*bw+i;
				Urowi++;
			}
		}
		for(map<int,double*>::const_iterator p=(++r[j].begin());
				p!=r[j].end();p++){
			for(int i=0;i<bw;i++){
				if(fabs((p->second)[i])>1.e-16){
					Uval[Urowi]=(p->second)[i];
					Uci[Urowi]=p->first*bw+i;
					Urowi++;
				}
			}
		}
	}
	Lri[nrows]=Llvc;
	Uri[nrows]=Ulvc;
	return crLU(Llvc,Ulvc,Lval,Uval,Lci,Uci,Lri,Uri,diag);
}

