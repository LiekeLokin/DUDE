// flow.cpp

#include "flow.h"
#include "admin.h"

using namespace std;
using namespace admin;


flow::flow(int Npx, int Npz) : Npx(Npx), Npz(Npz), nt(Npx*(Npz+1)) {
	A=new spMat(nt,Npz+1);
	prec=0;
	b=new vec(nt,0.0);
	iu=new vec(nt,0.0);
	alpha=new vec(Npx*Npz,0.0);
	beta=new vec(2*Npx,0.0);
	u0=new vec(Npz,0.0);
	Avx=new vec(Npx,Av);
	Sx=new vec(Npx,S);

	/*
	Am=new spMat(Npz,Npz);

	//0 oplossing numeriek uitrekenen
	//viscositeit linkerlid u vgl, op bodem, partial slip
	Am->e(0,1)-=Av/(dz*dz);
	double a2=Av;
	double b2=S*dz;
	Am->e(0,0)-=((2.0*a2-b2)/(2.0*a2+b2)-2.0)*Av/(dz*dz);
	for(int j=1;j<Npz-1;j++){
		Am->e(j,(j+1))-=Av/(dz*dz);
		Am->e(j,j)+=2.0*Av/(dz*dz);
		Am->e(j,(j-1))-=Av/(dz*dz);
	}
	//viscositeit met du/dz=0 aan oppervlakte
	Am->e((Npz-1),(Npz-1))+=Av/(dz*dz);
	Am->e((Npz-1),(Npz-2))-=Av/(dz*dz);
	for(int j=0;j<Npz;j++){
		(*u0)[j]=F;
	}
	//Am->GaussElim(*u0);
	crLU lu=A->LU();
	lu.bf(*u0);
	*/

}

flow::~flow(){
	delete A;
	delete b;
	delete iu;
	delete alpha;
	delete beta;
	delete u0;
	delete Avx;
	delete Sx;
	/*delete Am;*/
}
void flow::reprec(){
	prec=0;
}

int flow::o(int j_ex,int i_ex,int v) const {
	/*geeft het adres van component a, van var. v bij x is i en z is j*/
	int j,i = 0;
	if(v==1)cerr<<"w staat niet op de kaart !!!"<<endl;
	if(v==2){
		if(j_ex!=0)cout<<"zeta te hoog gegrepen"<<endl;
		j=0;
		v=1;
	}
	else j=j_ex;
	if(j<0||j>=Npz)cout<<"j mag natuurlijk geen "<<j<<"zijn"<<endl;
	if(i_ex>=0&&i_ex<Npx)i=i_ex;
	else{
		if(i_ex==-1)i=Npx-1;
		else{
			if(i_ex==Npx)i=0;
			else cout<<"Te groot "<<i_ex<<" , "<<j<<endl;
		}
	}
	return int(i*(Npz+1)+j+v*Npz);
}
void flow::det_AvS(vec bottom_state){
	//bsp = bottom_state with parameterization
	for(int i=0;i<Npx;i++){
		double loc_wd=-bottom_state[i]+H+(*iu)[o(0,i,2)]/(*beta)[2*i+1];
		(*Avx)[i]=Av;
		(*Sx)[i]=S;
		//cerr<<"loc_wd: "<<loc_wd<<"; ustar: "<<ustar<<"; Av: "<<(*Avx)[i]<<"; S: "<<(*Sx)[i]<<endl;
	}
}

vec flow::getiu(){
	/*
	vec u(nt,0.0);
	for(int i=0;i<Npx;i++){
		u[o(0,i,2)]=(*iu)[o(0,i,2)]/(*beta)[2*i+1];
		for(int j=0;j<Npz;j++){
			u[o(j,i,0)]=(*iu)[o(j,i,0)]*(*beta)[2*i];
		}
	}
	return u;
	*/
	return *iu;
}

double flow::check_qsp(){
	vec u_gem(Npx,0.0);
	for(int i=0;i<Npx;i++){
	  for(int j=0;j<Npz;j++){
		u_gem[i]+=(*iu)[o(j,i,0)]/Npz;
	  }
	  u_gem[i]*=H;
	  u_gem[i]+=(*iu)[o(Npz-1,i,0)]*((*iu)[o(0,i,2)]+(*iu)[o(0,o2(i-1),2)])/2.;
	}
	double q_sp = 0.0;
	for(int i=0;i<Npx;i++) q_sp+=u_gem[i];
	q_sp/=double(Npx);

	//cerr<<"qsp check: "<<q_sp<<" ("<<u_gem[0]<<" "<<u_gem[5]<<" "<<u_gem[10]<<" "<<u_gem[20]<<" "<<u_gem[24]<<")"<<endl;
	return q_sp;
}

void flow::resetIu(){
	for(int i=0;i<nt;i++)(*iu)[i]=0.0;
}

void flow::resetIu(vec u){
	for(int i=0;i<nt;i++)(*iu)[i]=u[i];
}

void flow::initIu(){
	for(int i=0;i<Npx;i++){
		for(int j=0;j<Npz;j++){
			(*iu)[o(j,i,0)]=(*u0)[j];
		}
		(*iu)[o(0,i,2)]=0.0;
	}
}

void flow::u_b(vec &u0){
	for(int i=0;i<Npx;i++){
		double a=2.*(*beta)[2*i]*(*Avx)[i];
		double fac=a/(a+(*Sx)[i]*dz);
		u0[i]=fac*(*iu)[o(0,i,0)];
	}
}

void flow::write_velocities(double tijd, vec bottom_state, vec u0_b){
	ofstream outvelub;
	ofstream outvelu;
	ofstream outvelw;
	outvelub.open ("out_velub.txt", ofstream::out | ofstream::app);
	outvelu.open  ("out_velu.txt", ofstream::out | ofstream::app);
	outvelw.open  ("out_velw.txt", ofstream::out | ofstream::app);
	outvelub.precision(10);
	outvelu.precision(10);
	outvelw.precision(10);

	outvelub<<tijd<<" ";
	outvelu<<tijd<<" ";
	outvelw<<tijd<<" ";

	for(int i=0;i<Npx;i++){
		outvelub<<u0_b[i]<<" ";
		if (i==Npx-1) outvelub<<endl;
	}

	double uu; double ww;
	for(int i=0;i<Npx;i++){
		for(int j=0;j<Npz;j++){
			//Gecolloceerd intern. Voor de onderste w heb je alpha op de bodem en w op de bodem nodig. w op de bodem is 0, dus dat is makkelijk.
			uu=0.5*((*iu)[o(j,i,0)]+(*iu)[o(j,o2(i+1),0)])*(*beta)[2*i];
			if (j>0) {
				ww=0.5*(w_from_u(j,i)+w_from_u(j-1,i)) -0.5*( (*alpha)[i*Npz+j]+(*alpha)[i*Npz+j-1]) *0.5*( (*iu)[o(j,i,0)]+(*iu)[o(j,i+1,0)] ); }
			else if (j==0) {
				double dh=(bottom_state[o2(i+1)]-bottom_state[i])/dx;
				double h=0.5*(bottom_state[o2(i+1)]+bottom_state[i]);
				double alpha_b=(-H)*dh/(H-h);
				ww=0.5*(w_from_u(j,i)+0) -0.5*( (*alpha)[i*Npz+j]+alpha_b) *0.5*( (*iu)[o(j,i,0)]+(*iu)[o(j,i+1,0)] ); }
			outvelu<<uu<<" ";
			outvelw<<ww<<" ";
			if (i==Npx-1 && j==Npz-1) {
				outvelu<<endl;
				outvelw<<endl;
			}
		}
	}


	outvelub.close();
	outvelu.close();
	outvelw.close();
}

/* replaced by writing via main
void flow::write_zeta(double tijd){
	ofstream outzeta;
	outzeta.open  ("out_zeta.txt", ofstream::out | ofstream::app);
	outzeta.precision(10);

	outzeta<<tijd<<" ";

	for(int i=0;i<Npx;i++){
		if ((*beta)[2*i+1]==0.) {outzeta<<0<<" ";}
		else {outzeta<<(*iu)[o(0,i,2)]/(*beta)[2*i+1]<<" ";}
	}
	outzeta<<endl;

	outzeta.close();
}
*/

vec flow::getZeta(){
	vec zta(Npx,0.0);
	for(int i=0;i<Npx;i++){
		if ((*beta)[2*i+1]==0.) {zta[i]=0.;}
		else {zta[i]=(*iu)[o(0,i,2)]/(*beta)[2*i+1];}
	}
	return zta;
}

double flow::zetaint1(){
	double zetaint = 0.0;
	for(int i=0;i<Npx;i++){
		zetaint+=dx*(*beta)[2*i+1]*(*iu)[o(0,i,2)];
	}
	return zetaint;
}

double flow::zetaint2(){
	double zetaint = 0.0;
	double minzeta = (*beta)[2*0+1]*(*iu)[o(0,0,2)];
	for(int i=1;i<Npx;i++){
		if ((*beta)[2*i+1]*(*iu)[o(0,i,2)]<minzeta) minzeta=(*beta)[2*i+1]*(*iu)[o(0,i,2)];}
	for(int i=0;i<Npx;i++){
		zetaint+=dx*((*beta)[2*i+1]*(*iu)[o(0,i,2)]-minzeta);
	}
	return zetaint;
}

/*
void flow::testNewton(vec bottom_state,double eps){
	//vergelijken van numeriek bepaalde en analytische Newton iterator
	vec b0(nt,0.0);
	vec iu0=(*iu);
	vec db(nt,0.0);
	spMtx An(nt);
	spMtx Av(nt);
	dzs_init(bottom_state);
	det_AvS(bottom_state);
	iu0=(*iu);
	vulb();
	vulA();
	b0=(*b);
	for(int i=0;i<nt;i++){
		(*iu)=iu0;
		(*iu)[i]+=eps;
		vulb();
		db=b0-(*b);
		for(int j=0;j<nt;j++){
			if(fabs(db[j])>1.e-14){
				An.e(j,i)=db[j]/eps;
				if(fabs((*A).e(j,i)-db[j]/eps)>1e-7)Av.e(j,i)=(*A).e(j,i)-db[j]/eps;
			}
		}
	}
	ofstream mtx("mtx.txt");
//cout<<"Av==========================="<<endl;
//cout<<Av<<endl;
//cout<<"An==========================="<<endl;
//cout<<An<<endl;
//cout<<"A============================"<<endl;
	mtx<<(*A)<<endl;
	mtx.close();
	ofstream of("rhs.txt");
//A->backForthLUi(*b);
	of.precision(16);
	for(int i=0;i<nt;i++)of<<(*b)[i]<<endl;
	of.close();
}
*/

int flow::solve_gm(vec bottom_state,int gmn){
	//flowsolver die gebruikt maakt van gmres, gmn is het aantal iteratiestappen
	double eta=0.1;
	double gmtresh=tresh;
	int teller=0;
	double resid=0.;
	dzs_init(bottom_state);
	det_AvS(bottom_state);
	//initIu(); //starten vanaf de vlakke bodem oplossing i.p.v. de vorige
	vulb();
	vulA();
	vec rhs=*b;
	vec x0=*iu;
	if(!prec){
		prec=1;
		prLU=new crLU(A->LU());
		vulA();
	}
	crMat crA(*A,nt);
	int gmrestel=0;
	//gmtresh=eta*b->L2(); //aanzetten voor incomplete newton
//	cerr<<"ik kom hier - 1"<<endl;
	double residu=gmres(x0,rhs,crA,*prLU,gmn);
//	cerr<<"ik kom hier - 2"<<endl;
	while(residu>gmtresh&&gmrestel<3){
		x0=rhs;
		rhs=*b;
		residu=gmres(x0,rhs,crA,*prLU,gmn);
		gmrestel++;
	}
	cerr<<"gmres restarted "<<gmrestel<<" times"<<endl;
	if(gmrestel>2){
		cerr<<"  WARNING: solve_gm, part 1, gmrestel="<<gmrestel<<endl;
		cerr<<"  now routine SOLVE called"<<endl;
		outlog<<"T="<<tijd<<" - WARNING: gmrestel="<<gmrestel<<" in solve_gm, part 1; routine solve is called!"<<endl;
		delete prLU;
		//prLU=new crLU(A->LU());
		//prLU->bf(*b);
		//cerr<<prec<<endl;
		teller=-1;
		prec=0;
		return teller;
	}
	else *b=rhs;
	vulu();
	teller++;
	while(L2(*b)>tresh&&teller<max_it){
		det_AvS(bottom_state);
		vulb();
		vulA();
		int gmrestel=0;
		x0=*iu;
		rhs=*b;
		//gmtresh=eta*b->L2(); //aanzetten voor incomplete newton
		double residu=gmres(x0,rhs,crA,*prLU,gmn);
		while(residu>gmtresh&&gmrestel<3){
			x0=rhs;
			rhs=*b;
			residu=gmres(x0,rhs,crA,*prLU,gmn);
			gmrestel++;
		}
		cerr<<"gmres restarted "<<gmrestel<<" times"<<endl;
		if(gmrestel>2){
			cerr<<"  WARNING.. solve_gm, part 2, gmrestel="<<gmrestel<<endl;
			cerr<<"  now routine SOLVE called"<<endl;
			outlog<<"T="<<tijd<<" - WARNING: gmrestel="<<gmrestel<<" in solve_gm, part 2; routine solve is called!"<<endl;
			delete prLU;
			//prLU=new crLU(A->LU());
			//prLU->bf(*b);
			//cerr<<prec<<endl;
			teller=-1;
			prec=0;
			return teller;
		}
		else *b=rhs;
		resid=L2(*b);
		cerr<<"Newton residu "<<resid<<endl;
		vulu();
		teller++;
	}
	if (resid>gmtresh){
		cerr<<"  WARNING: solve_gm residu > threshold, residu="<<resid<<endl;
		cerr<<"  now routine SOLVE called"<<endl;
		outlog<<"T="<<tijd<<" - WARNING: solve_gm residu > threshold, residu="<<resid<<endl;
		teller=-1;
	}
	return teller;
}

int flow::solve(vec bottom_state){
	//Flow solver met behulp van Gaus eliminatie
	int teller=0;
	//resetIu(); //*iu initialiseren op nul
	dzs_init(bottom_state);
	det_AvS(bottom_state);
	//initIu(); //starten vanaf de vlakke bodem oplossing i.p.v de vorige
	vulb();
	vulA();
	(A->LU()).bf(*b);
	//crLU lu=A->LU();
	//lu.bf(*b);
	vulu();
	teller++;
	cerr<<"Newton "<<teller<<", L2: "<<L2(*b)<<endl;
	while(L2(*b)>tresh&&teller<max_it){
		det_AvS(bottom_state);
		vulb();
		vulA();
		(A->LU()).bf(*b);
		//crLU lu=A->LU();
		//lu.bf(*b);
		vulu();
		teller++;
		cerr<<"Newton "<<teller<<", L2: "<<L2(*b)<<endl;
	}
	return teller;
}

void flow::vulu(){
	//Newton update
	for(int i=0;i<nt;i++)(*iu)[i]+=(*b)[i];
}

void flow::dzs_init(vec bottom_state){
	//vullen van de diverse bodemafhankelijke schaalparameters
	for(int i=0;i<Npx;i++){
		(*beta)[2*i]=1./(1.-(bottom_state[i])/H);
		(*beta)[2*i+1]=1./(1.-0.5*(bottom_state[o2(i+1)]+bottom_state[i])/H);
		double dh=(bottom_state[o2(i+1)]-bottom_state[i])/dx;
		double h=0.5*(bottom_state[o2(i+1)]+bottom_state[i]);
		for(int j=0;j<Npz;j++)(*alpha)[i*Npz+j]=((j+1.0)*dz-H)*dh/(H-h);
	}
}

void flow::zl(int r,int i){
	//vullen van het linkerlid van de zeta vgl
	double fac=1.0/(2.0*dx);
	w_from_uA(r,Npz-1,i,-1.0);
	A->e(r,o(Npz-1,i+1,0))+=fac*((*iu)[o(0,i+1,2)]+(*iu)[o(0,i,2)]);
	A->e(r,o(0,i+1,2))+=fac*(*iu)[o(Npz-1,i+1,0)];
	A->e(r,o(0,i,2))+=fac*(*iu)[o(Npz-1,i+1,0)];
	A->e(r,o(Npz-1,i,0))+=-fac*((*iu)[o(0,i-1,2)]+(*iu)[o(0,i,2)]);
	A->e(r,o(0,i-1,2))+=-fac*(*iu)[o(Npz-1,i,0)];
	A->e(r,o(0,i,2))+=-fac*(*iu)[o(Npz-1,i,0)];

}

void flow::zr(int r,int i){
	//vullen van het rechterlid van de zeta vgl
	double fac=1.0/(2.0*dx);
	(*b)[r]+=1.0*w_from_u(Npz-1,i);
	(*b)[r]-=fac*(*iu)[o(Npz-1,i+1,0)]*((*iu)[o(0,i+1,2)]+(*iu)[o(0,i,2)]);
	(*b)[r]-=fac*-(*iu)[o(Npz-1,i,0)]*((*iu)[o(0,i-1,2)]+(*iu)[o(0,i,2)]);
}

void flow::w_from_uA(int r, int j_ex,int i_ex,double co){
	for(int j=0;j<=j_ex;j++){
		A->e(r,o(j,i_ex,0))+=co*dz/dx;
		A->e(r,o(j,i_ex+1,0))-=co*dz/dx;
	}
}
double flow::w_from_u(int j_ex,int i_ex){
	double w=0.;
	for(int j=0;j<=j_ex;j++){
		w+=dz/dx*((*iu)[o(j,i_ex,0)]-(*iu)[o(j,i_ex+1,0)]);
	}
	return w;
}


void flow::a0l(int r,int j,int i){
	/*eerste advectieterm (du^2/dx) in linkerkant u vgl*/
	double fac=1.0/(2.0*dx);
	A->e(r,o(j,i,0))+=fac*(*beta)[2*i+1]*((*iu)[o(j,i,0)]+(*iu)[o(j,i+1,0)]);
	A->e(r,o(j,i+1,0))+=fac*(*beta)[2*i+1]*((*iu)[o(j,i,0)]+(*iu)[o(j,i+1,0)]);

	A->e(r,o(j,i-1,0))-=fac*(*beta)[2*o2(i-1)+1]*((*iu)[o(j,i-1,0)]+(*iu)[o(j,i,0)]);
	A->e(r,o(j,i,0))-=fac*(*beta)[2*o2(i-1)+1]*((*iu)[o(j,i-1,0)]+(*iu)[o(j,i,0)]);
}

void flow::a0r(int r,int j,int i){
	/*eerste advectieterm (du^2/dx) in rechterkant u vgl*/
	double fac=1.0/(4.0*dx);
	(*b)[r]-=fac*(*beta)[2*i+1]*((*iu)[o(j,i,0)]+(*iu)[o(j,i+1,0)])*((*iu)[o(j,i,0)]+(*iu)[o(j,i+1,0)]);
	(*b)[r]-=-fac*(*beta)[2*o2(i-1)+1]*((*iu)[o(j,i-1,0)]+(*iu)[o(j,i,0)])*((*iu)[o(j,i-1,0)]+(*iu)[o(j,i,0)]);
}

void flow::a1l(int r,int j,int i){
	/*tweede advectieterm (dwu/dz) in linkerkant van u vgl*/
	double fac=1.0/(4.0*dz);
	fac*=(*beta)[2*i];
	w_from_uA(r,j,i,fac*((*iu)[o(j+1,i,0)]+(*iu)[o(j,i,0)]));
	w_from_uA(r,j,i-1,fac*((*iu)[o(j+1,i,0)]+(*iu)[o(j,i,0)]));
	w_from_uA(r,j-1,i,-fac*((*iu)[o(j,i,0)]+(*iu)[o(j-1,i,0)]));
	w_from_uA(r,j-1,i-1,-fac*((*iu)[o(j,i,0)]+(*iu)[o(j-1,i,0)]));


	A->e(r,o(j+1,i,0))+=fac*(w_from_u(j,i)+w_from_u(j,i-1));
	A->e(r,o(j,i,0))+=fac*(w_from_u(j,i)+w_from_u(j,i-1)-w_from_u(j-1,i)-w_from_u(j-1,i-1));
	A->e(r,o(j-1,i,0))-=fac*(w_from_u(j-1,i)+w_from_u(j-1,i-1));

}

void flow::a1r(int r,int j,int i){
	/*tweede advectieterm (dwu/dz) in rechterkant van u vgl*/
	double fac=1.0/(4.0*dz);
	fac*=(*beta)[2*i];
	(*b)[r]-=fac*((*iu)[o(j+1,i,0)]+(*iu)[o(j,i,0)])*(w_from_u(j,i)+w_from_u(j,i-1));
	(*b)[r]+=fac*((*iu)[o(j,i,0)]+(*iu)[o(j-1,i,0)])*(w_from_u(j-1,i)+w_from_u(j-1,i-1));
}

void flow::a1l_opperv(int r,int i){
	/*tweede advectieterm met du/dz=0 aan opp. linkerkant u vgl*/
	double fac=1.0/(4.0*dz);
	fac*=(*beta)[2*i];
	w_from_uA(r,Npz-1,i,fac*(2.0*(*iu)[o(Npz-1,i,0)]));
	w_from_uA(r,Npz-1,i-1,fac*(2.0*(*iu)[o(Npz-1,i,0)]));
	w_from_uA(r,Npz-2,i,-fac*((*iu)[o(Npz-1,i,0)]+(*iu)[o(Npz-2,i,0)]));
	w_from_uA(r,Npz-2,i-1,-fac*((*iu)[o(Npz-1,i,0)]+(*iu)[o(Npz-2,i,0)]));

	A->e(r,o(Npz-1,i,0))+=fac*(w_from_u(Npz-1,i)+w_from_u(Npz-1,i-1));
	A->e(r,o(Npz-1,i,0))+=fac*(w_from_u(Npz-1,i)+w_from_u(Npz-1,i-1)-w_from_u(Npz-2,i)-w_from_u(Npz-2,i-1));
	A->e(r,o(Npz-2,i,0))-=fac*(w_from_u(Npz-2,i)+w_from_u(Npz-2,i-1));
}

void flow::a1r_opperv(int r,int i){
	/*tweede advectieterm met du/dz=0 aan opp. rechterkant u vgl*/
	double fac=1.0/(4.0*dz);
	fac*=(*beta)[2*i];
	(*b)[r]-=fac*(2.0*(*iu)[o(Npz-1,i,0)])*(w_from_u(Npz-1,i)+w_from_u(Npz-1,i-1));
	(*b)[r]+=fac*((*iu)[o(Npz-1,i,0)]+(*iu)[o(Npz-2,i,0)])*(w_from_u(Npz-2,i)+w_from_u(Npz-2,i-1));
}

void flow::a1l_bodem(int r,int i){
	/*dwu/dz met partial slip op bodem, linkerkant u vgl*/
	double fac=1.0/(4.0*dz);
	fac*=(*beta)[2*i];
	w_from_uA(r,0,i,fac*((*iu)[o(1,i,0)]+(*iu)[o(0,i,0)]));
	w_from_uA(r,0,i-1,fac*((*iu)[o(1,i,0)]+(*iu)[o(0,i,0)]));

	A->e(r,o(1,i,0))+=fac*(w_from_u(0,i)+w_from_u(0,i-1));
	A->e(r,o(0,i,0))+=fac*(w_from_u(0,i)+w_from_u(0,i-1));
}

void flow::a1r_bodem(int r,int i){
	/*dwu/dz met no slip op bodem, rechterkant u vgl*/
	double fac=1.0/(4.0*dz);
	fac*=(*beta)[2*i];
	(*b)[r]-=fac*((*iu)[o(1,i,0)]+(*iu)[o(0,i,0)])*(w_from_u(0,i)+w_from_u(0,i-1));
}

void flow::gzl(int r,int i){
	/*g dzeta/dx linkerkant u vgl*/
	//A->e(r,o(0,i,2))+=g/((*beta)[2*i]*dx)/(*beta)[2*i+1];
	//A->e(r,o(0,i-1,2))-=g/((*beta)[2*i]*dx)/(*beta)[2*o2(i-1)+1];
	A->e(r,o(0,i,2))+=g/((*beta)[2*i]*dx);
	A->e(r,o(0,i-1,2))-=g/((*beta)[2*i]*dx);
}

void flow::gzr(int r,int i){
	/*g dzeta/dx rechterkant u vgl*/
	//(*b)[r]-=g/((*beta)[2*i]*dx)*((*iu)[o(0,i,2)]/(*beta)[2*i+1]-(*iu)[o(0,i-1,2)]/(*beta)[2*o2(i-1)+1]);
	(*b)[r]-=g/((*beta)[2*i]*dx)*((*iu)[o(0,i,2)]-(*iu)[o(0,i-1,2)]);
}

void flow::vl(int r,int j,int i){
	/*viscositeit linkerlid u vgl*/
	A->e(r,o(j+1,i,0))-=(*beta)[2*i]*(*beta)[2*i]*(*Avx)[i]/(dz*dz);
	A->e(r,o(j,i,0))+=2.0*(*beta)[2*i]*(*beta)[2*i]*(*Avx)[i]/(dz*dz);
	A->e(r,o(j-1,i,0))-=(*beta)[2*i]*(*beta)[2*i]*(*Avx)[i]/(dz*dz);
}

void flow::vr(int r,int j,int i){
	/*viscositeit rechterlid u vgl*/
	(*b)[r]+=(*beta)[2*i]*(*beta)[2*i]*(*Avx)[i]*(*iu)[o(j+1,i,0)]/(dz*dz);
	(*b)[r]-=2.0*(*beta)[2*i]*(*beta)[2*i]*(*Avx)[i]*(*iu)[o(j,i,0)]/(dz*dz);
	(*b)[r]+=(*beta)[2*i]*(*beta)[2*i]*(*Avx)[i]*(*iu)[o(j-1,i,0)]/(dz*dz);
}

void flow::vl_bodem(int r,int i){
	/*viscositeit linkerlid u vgl, op bodem, partial slip*/
	A->e(r,o(1,i,0))-=(*beta)[2*i]*(*beta)[2*i]*(*Avx)[i]/(dz*dz);
	double a2=(*Avx)[i]*(*beta)[2*i];
	double b2=(*Sx)[i]*dz;
	A->e(r,o(0,i,0))-=(*beta)[2*i]*(*beta)[2*i]*((2.0*a2-b2)/(2.0*a2+b2)-2.0)*(*Avx)[i]/(dz*dz);
}

void flow::vr_bodem(int r,int i){
	/*viscositeit rechterlid u vgl, op bodem, partial slip*/
	(*b)[r]+=(*beta)[2*i]*(*beta)[2*i]*(*Avx)[i]*(*iu)[o(1,i,0)]/(dz*dz);
	double a2=(*Avx)[i]*(*beta)[2*i];
	double b2=(*Sx)[i]*dz;
	(*b)[r]+=(*beta)[2*i]*(*beta)[2*i]*(*Avx)[i]*(*iu)[o(0,i,0)]*((2.0*a2-b2)/(2.0*a2+b2)-2.0)/(dz*dz);
}

void flow::vl_opperv(int r,int i){
	/*viscositeit met du/dz=0 aan oppervlakte*/
	A->e(r,o(Npz-1,i,0))+=(*beta)[2*i]*(*beta)[2*i]*(*Avx)[i]/(dz*dz);
	A->e(r,o(Npz-2,i,0))-=(*beta)[2*i]*(*beta)[2*i]*(*Avx)[i]/(dz*dz);
}

void flow::vr_opperv(int r,int i){
	/*viscositeit met du/dz=0 aan oppervlakte*/
	(*b)[r]-=(*beta)[2*i]*(*beta)[2*i]*(*Avx)[i]*(*iu)[o(Npz-1,i,0)]/(dz*dz);
	(*b)[r]+=(*beta)[2*i]*(*beta)[2*i]*(*Avx)[i]*(*iu)[o(Npz-2,i,0)]/(dz*dz);
}

void flow::last_zetar(){
  /*de 0 component van zeta moet geintegreerd over domein 0 zijn*/
   (*b)[o(0,Npx-1,2)]=0.0;
   //for(int i=0;i<Npx;i++)(*b)[o(0,Npx-1,2)]+=(*iu)[o(0,i,2)];
	for(int i=0;i<Npx;i++)(*b)[o(0,Npx-1,2)]+=(*beta)[2*i+1]*(*iu)[o(0,i,2)];
}


void flow::last_zetal(){
	int r=o(0,Npx-1,2);
	A->rmrow(r);
	//for(int i=0;i<Npx;i++)A->e(r,o(0,i,2))=-1.0;
	for(int i=0;i<Npx;i++)A->e(r,o(0,i,2))=-(*beta)[2*i+1];
}

void flow::vulb(){
	/*vullen van de rechterkant van [A]*{x}={b}*/
	int r;
	for(int i=0;i<nt;i++)(*b)[i]=0.0;
	for(int i=0;i<Npx;i++){
		r=o(0,i,0);
		a0r(r,0,i);
		a1r_bodem(r,i);
		gzr(r,i);
		vr_bodem(r,i);
		(*b)[r]+=F/(*beta)[2*i];
		for(int j=1;j<(Npz-1);j++){
			r=o(j,i,0);
			a0r(r,j,i);
			a1r(r,j,i);
			gzr(r,i);
			vr(r,j,i);
			(*b)[r]+=F/(*beta)[2*i];
		}
		r=o(Npz-1,i,0);
		a0r(r,Npz-1,i);
		a1r_opperv(r,i);
		gzr(r,i);
		vr_opperv(r,i);
		zr(o(0,i,2),i);
		(*b)[r]+=F/(*beta)[2*i];
	}
	last_zetar();
}

void flow::vulA(){
	/*vullen van de linkerkant van [A]*{x}={b}*/
	int r;
	A->softEmpty();
	for(int i=0;i<Npx;i++){
		r=o(0,i,0);
		a0l(r,0,i);
		a1l_bodem(r,i);
		gzl(r,i);
		vl_bodem(r,i);
		for(int j=1;j<(Npz-1);j++){
			r=o(j,i,0);
			a0l(r,j,i);
			a1l(r,j,i);
			gzl(r,i);
			vl(r,j,i);
		}
		r=o(Npz-1,i,0);
		a0l(r,Npz-1,i);
		a1l_opperv(r,i);
		gzl(r,i);
		vl_opperv(r,i);
		zl(o(0,i,2),i);
	}
	last_zetal();
}

