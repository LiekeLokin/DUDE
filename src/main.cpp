// main.cpp

#include <iostream>
#include <fstream>
#include <vector>
#include "bottom.h"
#include "flow.h"
#include "linalg.h"
#include "admin.h"
#include "Config.h"
#include <ctime>
#include <cstring>


#include <sstream>
#include <string>

using namespace admin;
using namespace std;

ofstream outlog("out_log.txt");

double H;
double dt;
double L;
double dx;
double dz;
double tijd;
double Av;
double S;

void copyConfigToAdmin(const Config& cfg);
void doStabAnalysis(int stabWrite, flow& H2O, bottom& sand, double& q_in);
void doCheckQsp(vec bedflow, flow& H2O, double& q_in);
void setS_Av();
double maxval(vec vinp);

int main (int argc, char * const argv[]) {

std::string filename = (argc == 1) ? "config.cfg" : argv[1];
Config cfg(filename);

copyConfigToAdmin(cfg);

// Initialize global variables
H = H0;
dt = dtr;
L = 1.0;
dx = L / Npx;
dz = H / Npz;
tijd = 0.0;
Av = 0.0;
S = 0.0;

double q_in = q_in1;

flow H2O(Npx, Npz);
bottom sand(Npx);

if (dt_write==1.) {cerr<<endl<<endl<<endl<<endl<<endl<<"         ------ NOTE!! DT_WRITE==1!! -------"<<endl<<endl<<endl<<endl<<endl;}

for (int p=1;p<=1;p++){				//superloop!!!!!!!!!!!!

	cerr.precision(16);
	//vec current(Npx,0.0);
	//vec bedflow(Npx,0.0);
	//vec next(Npx,0.0);
	//vec u(nt,0.0); // hier zit "echte" snelheid in
	//vec u0_b(Npx,0.0);
	//vec bss1(Npx,0.0);
	//vec bss2(Npx,0.0);
	//vec fluxtot(Npx,0.0);
	//vec dhdx(Npx,0.0);
	int n_it_fl=0;
	//int Nd=0; // number of dunes in domain
	//double norm=0.0;
	//double bint1; double bint2; double zetaint1; double zetaint2;
	//int solve_method=0; // determine in CheckFlowsep (=1 if bottom changes strongly)
	int sepflag=0;
	//int nfsz=0;
	//cerr<<nt<<" is de dimensie van de matrix"<<endl;
	cerr<<F<<endl;
	sand.setSin(ampbeds,1);
//	sand.setRand(0.1*D50,(unsigned)time(0));
	//sand.setRand(1e-8,(unsigned)time(0));
	//sand.setShape(current);
	setS_Av();
	H2O.resetIu();
	if (readbed==1) {cerr<<endl<<endl<<endl<<endl<<endl<<"         ------>> NOTE !! Bed elevation read from file !! <<-------"<<endl<<endl<<endl<<endl<<endl;}
	int iinit1=0;
	//vec inp(3,0.0); //is input if a bottom profile is read
	
	// NEW STUFF, WRITTEN BY OLAV 2014 for FLOODWAVE
	vec fw_t(2,0.0);
	vec fw_q(2,0.0);
	if (readfw==1){
		double temp;
		ifstream in1("floodwave.inp");
		in1>>temp;
		int np = int(temp);
		fw_t.resize(np);
		fw_q.resize(np);
		cerr<< "Number of points in floodwave: " << np<<endl;
		for (int j=0;j<np;j++){	in1>>temp;
								fw_t[j]=double(temp);
								in1>>temp;
								fw_q[j]=double(temp);
								cerr << fw_t[j] << " " << fw_q[j] << endl; }
	}// NEW STUFF, WRITTEN BY OLAV 2014 for FLOODWAVE
	
	if (readbed==1) {
		auto inp=sand.readBottomInp();
		H=inp[1];
				
		dz=H/Npz;
		tijd=inp[0];
		
		if (readfw==1){q_in=interpolate(fw_t,fw_q,tijd);} // NEW STUFF, WRITTEN BY OLAV 2014 for FLOODWAVE
		
		iinit1=int(tijd/dt);
		L=inp[2];
						
		dx=L/Npx;
		setS_Av();
		cerr<<"read check:"<<endl;
		cerr<<"Timeprevious: " <<tijd<<endl;
		cerr<<"H: "<<H<<endl;
		cerr<<"L: "<<L<<endl;
		cerr<<"first bed point: "<<sand.getShape(0)[0]<<endl;
		cerr<<"last bed point: "<<sand.getShape(0)[Npx-1]<<endl;
		cerr<<"dx: "<<dx<<endl;
		//cerr<<"Av: "<<Av<<endl;
		//cerr<<"S : "<<S<<endl;
		
		//avalanching protocol
		// if(transport_eq == 2 && AllowAvalanching == 1){
		// avalanche(); 
		// }
	}
	else {
		if (readfw==1){
			q_in=interpolate(fw_t,fw_q,tijd);
		}
	}// NEW STUFF, WRITTEN BY OLAV 2014 for FLOODWAVE
	
	cerr << endl << "first discharge: " << q_in << endl;

	ostringstream tmpbot;
	tmpbot << "out_bottom"<< p << ".txt";
	string ofname = tmpbot.str();
	ofstream outbot(ofname.c_str(),ios_base::out);

	ostringstream tmpbot1;
	tmpbot1 << "out_general"<< p << ".txt";
	string ofname1 = tmpbot1.str();
	ofstream data(ofname1.c_str(),ios_base::out);

	ostringstream tmpbot2;
	tmpbot2 << "out_int"<< p << ".txt";
	string ofname2 = tmpbot2.str();
	ofstream outint(ofname2.c_str(),ios_base::out);

	ostringstream tmpbot3;
	tmpbot3 << "out_bss"<< p << ".txt";
	string ofname3 = tmpbot3.str();
	ofstream outbss(ofname3.c_str(),ios_base::out);

	ostringstream tmpbot4;
	tmpbot4 << "out_fsz"<< p << ".txt";
	string ofname4 = tmpbot4.str();
	ofstream outfsz(ofname4.c_str(),ios_base::out);

	ostringstream tmpbot5;
	tmpbot5 << "out_flux"<< p << ".txt";
	string ofname5 = tmpbot5.str();
	ofstream outflux(ofname5.c_str(),ios_base::out);

	ostringstream tmpbot6;
	tmpbot6 << "out_dhdx"<< p << ".txt";
	string ofname6 = tmpbot6.str();
	ofstream outdhdx(ofname6.c_str(),ios_base::out);

	ostringstream tmpbot7;
	tmpbot7 << "out_Sround"<< p << ".txt";
	string ofname7 = tmpbot7.str();
	ofstream outSround(ofname7.c_str(),ios_base::out);

	ostringstream tmpbot8;
	tmpbot8 << "out_zeta"<< p << ".txt";
	string ofname8 = tmpbot8.str();
	ofstream outzeta(ofname8.c_str(),ios_base::out);

	cout.precision(16);
	outbss.precision(10);
	outflux.precision(10);
	outdhdx.precision(10);
	outint.precision(16);
	outSround.precision(16);
	outzeta.precision(10);
	
	// ADDED 2011 2 25 (OLAV)
	//    ostringstream tmpbot10;
	//			tmpbot10 << "out_debug"<< p << ".txt";
	//			string ofname10 = tmpbot10.str();
	//	ofstream outdebug(ofname10.c_str(),ios_base::out);
	//	outdebug.precision(10);
	//	outdebug.close();
	// end ADDED 2011 2 25 (OLAV)

	data<<H<<endl<<L<<endl<<Npx<<endl<<dx<<endl<<Npz<<endl<<dz<<endl<<dt<<endl<<dt_write<<endl<<Av<<endl<<S<<endl<<q_in<<endl<<F<<endl<<1<<endl<<alpha<<endl<<be<<endl<<l1<<endl<<l2<<endl<<nd<<endl<<readfw;
	data.close();

	auto current=sand.getShape(sepflag);

	int write_teller = int(dt_write/dt);
	int cor=1; if (dt_write==dt) cor=0;
	int stabWrite=1; // this is set to zero if initial results is written to file
	/* OLD STUFF, WRITTEN BEFORE OLAV AND SULEYMAN
	// THIS PIECE FOR FLOODWAVE SIMULATION
	// read floodwave
	double temp;
	ifstream in1("floodwave.inp");
	in1>>temp;
	int np = int(temp);
	vec qa(np,0.0);
	cerr<<np<<endl;
	if (double(np-1)*dt != tend) {
		cerr<<"something is wrong with floodwave"<<endl;
		cerr<<double(np-1)*dt<<" "<<tend<<endl;
		return 1;}
	else {// read floodwave
		//for (int j=0;j<np;j++) in1>>qa[j];
		for (int j=0;j<np;j++) qa[j]=q_in;
		cerr<<qa[0]<<" "<<qa[1]<<" "<<qa.maxval()<<" "<<qa[np-2]<<" "<<qa[np-1]<<endl;}
	//end reading floodwave
	//OLD STUFF, WRITTEN BEFORE OLAV AND SULEYMAN
	*/

	auto bedflow=current;

  	struct tm  ts;
  	time_t     now;
  	char       buf[80];
  	std::time(&now);
	ts = *std::localtime(&now);
	std::strftime(buf, sizeof(buf), "%a %Y-%b-%d %H:%M:%S", &ts);
  	outlog<<"Simulation started at "<<buf<<endl;

	for(int i=iinit1;i<=tend/dt;i++){ // 40 minutes
		
		if (readfw==1){q_in=interpolate(fw_t,fw_q,tijd);} // NEW STUFF, WRITTEN BY OLAV 2014

		//OLAV: 2011 02 21 changed from 
        //Hdiff=(H-H0)/H0;
		const auto Hdiff=abs((H-H0)/H0); //equals 0 when exactly the same, 1 when the difference is 100%
		const auto Hcrit = Hcrit_global;
		cerr<<"Hdiff = " <<Hdiff <<" (Hcrit = "<<Hcrit<<")"<<endl;

		bool doStab=0;
		if (Hdiff>=Hcrit ) doStab=1; //OLAV: 2011 02 21 changed from 
        //if (Hdiff>=1.+Hcrit) doStab=1;
        
		/* start with initial stability analysis, or if H is sufficiently changed */
		//if ( (i==iinit1&&readbed1==0) | doStab==1 | i==1000) { // OLAV TEST 2011 2 23
		if (SimpleLength==0) { // OLAV 2012 09 06: added simple length implementation
		   if ( (i==iinit1&&readbed==0) || doStab==1) {
			  outlog<<"T="<<tijd<<" - WARNING: Stability Analysis. (Hdiff="<<Hdiff<<")"<<endl;
			  doStabAnalysis(stabWrite, H2O, sand, q_in);
			  stabWrite=0;
			  //OLAV: 2013 02 06 added doStab=0;
			  doStab=0; 
			  }
         }
        else if (SimpleLength==1) {
			//L=1.6;
            L = H*SimpleLengthFactor;
			//L=319.;
			//L=4.;
			
            dz=H/Npz;
	        dx=L/Npx;
			//cerr << "dx: " << dx << endl;
	        setS_Av();
	        dt=dtr;
			doStab=0;
        }

		//q_in=qa[i];
	  	
	  	//2012 09 13 Olav 
		if(AllowFlowSep == 1){ 
			sand.checkFlowsep(); 
		}

		//2012 09 13 Olav
		const auto& stateFsz=sand.getFsz();
		const auto nf = stateFsz.size();
		int solve_method=stateFsz[nf-3];
		sepflag=stateFsz[nf-2];
		int nfsz=stateFsz[nf-1];//}

		current=sand.getShape(0); // for determination of migration rate
		bedflow=sand.getShape(sepflag);

		if (i==iinit1) {
			cerr<<"Eerste keer stroming oplossing met SOLVE"<<endl;
			outlog<<"Eerste keer stroming oplossing met SOLVE"<<endl;
			n_it_fl=H2O.solve(bedflow);
	  		cerr<<"number of flow iteration required: "<<n_it_fl<<endl;}
		else if (doStab==1) {
			cerr<<"Stroming oplossing met SOLVE vanwege Stability Analysis"<<endl;
			outlog<<"T="<<tijd<<" - WARNING: SOLVE vanwege Stability Analysis. (Hdiff="<<Hdiff<<")"<<endl;
			H2O.resetIu();
			n_it_fl=H2O.solve(bedflow);
	  		cerr<<"number of flow iteration required: "<<n_it_fl<<endl;}
		else if (solve_method>=1) {
			cerr<<"Stroming oplossing met SOLVE vanwege sterke 'bodem' veranderingen"<<endl;
			cerr<<"solve_method= "<<solve_method<<endl;
			H2O.resetIu();
			n_it_fl=H2O.solve(bedflow);
	  		cerr<<"number of flow iteration required: "<<n_it_fl<<endl;}
		else {
		n_it_fl=H2O.solve_gm(bedflow,20); }
		if (n_it_fl==-1) {
			ofstream outtemp("out_temp.txt");
			vec bedtemp(Npx,0.);
			bedtemp=sand.getShape(0);
			for(int i=0;i<bedtemp.size();i++) { outtemp<<bedtemp[i]<<" "; } outtemp<<endl;
			bedtemp=sand.getShape(1);
			for(int i=0;i<bedtemp.size();i++) { outtemp<<bedtemp[i]<<" "; } outtemp<<endl;
			outtemp.close();
			H2O.resetIu();
			H2O.solve(bedflow);
		}

		doCheckQsp(bedflow, H2O, q_in);
		vec u0_b(Npx);
	  	H2O.u_b(u0_b);
		
		//cerr << "current: " << current[0] << " bedflow: " << bedflow[0] << endl; //OLAV 2014 03 31

	  	if(write_teller==dt_write/dt){
			//H2O.write_velocities(tijd,sand.getShape(sepflag),u0_b); //TEST: uncommented at 2011 03 21 (OLAV)
			//H2O.write_zeta(tijd); 
			const auto& stateZeta=H2O.getZeta();
			outzeta<<tijd<<" ";	for(int i=0;i<stateZeta.size();i++)outzeta<<stateZeta[i]<<" "; outzeta<<endl;
			//2014 01 27: changed q_sp to q_in
			outbot<<tijd<<" "<<sepflag<<" "<<nfsz<<" "<<q_in<<" "<<H<<" "<<L<<" "; for(int i=0;i<current.size();i++)outbot<<current[i]<<" "; outbot<<endl;
			vec bot=sand.getShape(sepflag);
			//2014 01 27: changed q_sp to q_in
			outbot<<tijd+0.1<<" "<<sepflag<<" "<<nfsz<<" "<<q_in<<" "<<H<<" "<<L<<" "; for(int i=0;i<bot.size();i++)outbot<<bot[i]<<" "; outbot<<endl;
		} 
		
		//cerr << "bint1: " << bint1 << " bint2: " << bint2 << endl; //OLAV 2014 03 31

		auto bint1=sand.detint1(current);
		auto bint2=sand.detint2(current);
		auto zetaint1=H2O.zetaint1();
		auto zetaint2=H2O.zetaint2();
		const auto& Dc=sand.detNd_fft(sand.getShape(0),2); // dune characteristics
		auto Nd=int(Dc[0]);
		double cr=Dc[1];
		double tr=Dc[2];
		double Hav=cr-tr;
		double Lav=L; 
		
		//cerr << "Hav: " << Hav << endl; //OLAV 2014 03 31

		vec next;
		vec bss1(Npx);
		vec bss2(Npx);
		vec fluxtot(Npx);
		vec dhdx(Npx);
		if (transport_eq==1 || transport_eq==3 ){
			if (sepflag==0){
				next=sand.update(u0_b,bss1,fluxtot,dhdx);
				bss2=bss1;}
			else if (sepflag==1){
				next=sand.update_flowsep(u0_b,bss1,bss2,fluxtot,dhdx);}
			}
		else if (transport_eq==2) {
			if (sepflag==0){
				next=sand.update(u0_b,bss1,fluxtot,dhdx);
				bss2=bss1;}
			else if (sepflag==1){
				next=sand.update_flowsep(u0_b,bss1,bss2,fluxtot,dhdx);}
		}
		
		//cerr << "next: " << next[1] << " bint2: " << bint2 << endl; //OLAV 2014 03 31
		
		double flux_av=0.;
		for(int i=0;i<fluxtot.size();i++) flux_av+=fluxtot[i];
		flux_av/=Npx;

		double mig=sand.detMigr(current,next);
		//outint<<tijd<<" "<<H<<" "<<bint1<<" "<<bint2<<" "<<zetaint1<<" "<<zetaint2<<endl;

	  	if(write_teller==dt_write/dt-cor){
			const auto& stateFsz=sand.getFsz();
			const auto& stateSr=sand.getSr();
			//wegschrijven fsz:
			outfsz<<tijd<<" ";		for(int i=0;i<stateFsz.size();i++)outfsz<<stateFsz[i]<<" "; outfsz<<endl;
			outSround<<tijd<<" "; for(int i=0;i<stateSr.size();i++)outSround<<stateSr[i]<<" "; outSround<<endl;
		}

	  	if(write_teller==dt_write/dt){
			outint<<tijd<<" "<<H<<" "<<bint1<<" "<<bint2<<" "<<zetaint1<<" "<<zetaint2<<" "<<cr<<" "<<tr<<" "<<Nd<<" "<<mig<<" "<<flux_av<<endl;
			outflux<<tijd<<" ";    for(int i=0;i<fluxtot.size();i++)outflux<<fluxtot[i]<<" "; outflux <<endl;
			outdhdx<<tijd<<" ";    for(int i=0;i<dhdx.size();   i++)outdhdx<<dhdx[i]   <<" "; outdhdx <<endl;
			outbss<<tijd<<" ";     for(int i=0;i<bss1.size();   i++)outbss <<bss1[i]   <<" "; outbss  <<endl;
			outbss<<tijd+0.1<<" "; for(int i=0;i<bss2.size();   i++)outbss <<bss2[i]   <<" "; outbss  <<endl;
			write_teller = 0;
		}

	  	const auto nt = next.size();
	  	vec verschil(nt);
	  	for(int k=0;k<nt;k++)verschil[k]=next[k]-current[k];
	  	auto norm=L2(verschil);
	  	sand.setShape(next);

	  	current=next;
	  	tijd+=dt;

	  	cerr<<"flowsolver "<<tijd<<" seconden ("<<tijd/60.<<" min)"<<"	onderweg"<<endl;
	  	cerr<<"number of flow iteration required: "<<n_it_fl<<endl;
	  	cerr<<"sepflag: "<<sepflag<<"; nfsz: "<<nfsz<<endl;
	  	cerr<<"Nd: "<<Nd;
	  	cerr<<"; wd: "<<H;
	  	cerr<<"; Lav: "<<Lav;
	  	cerr<<"; Hav: "<<Hav<<endl;
	  	cerr<<"integral of bed: "<<bint1<<endl;
	  	cerr<<"bodem ge update met (L2) : "<<norm<<" tot (L2) : "<<L2(next)<<endl<<endl;

	  	write_teller+=1;

	  	if(Hav<2.*ampbeds) {
	  		sand.setSin(ampbeds,1);
	  		cerr<<"Hav very low, bed set to initial disturbance."<<endl<<endl;
	  	}
	}
	
    sand.writeBottom();

	outfsz.close();
	outSround.close();
	outbss.close();
	outflux.close();
	outdhdx.close();
	outint.close();
	outzeta.close();
//	// ADDED 2011 2 22 (OLAV)
//	outdebug.close();
//	// end ADDED 2011 2 22 (OLAV)

	std::time(&now);
	ts = *std::localtime(&now);
	std::strftime(buf, sizeof(buf), "%a %Y-%b-%d %H:%M:%S", &ts);
	outlog<<"Simulation ended at "<<buf<<endl;

	outlog.close();

	char ofname9[15]="out_log";
	char pstr [10];
	sprintf(pstr,"%d",p);
	strcat ( ofname9, pstr );
	strcat ( ofname9, ".txt" );
	cerr<<ofname9<<endl;

  char oldname[] ="out_log.txt";
  rename( oldname , ofname9 );
	outlog.open("out_log.txt");
  outlog<<"run initialized; this file may be safely deleted"<<endl;

}
	return 0;
}

void copyConfigToAdmin(const Config& cfg) {
	DebugOutput = cfg.DebugOutput;
	Npx = cfg.Npx;
	Npz = cfg.Npz;
	dtr = cfg.dtr;
	dt_write = cfg.dt_write;
	tend = cfg.tend;
	ampbeds_factor = cfg.ampbeds_factor;
	AllowFlowSep = cfg.AllowFlowSep;
	AllowAvalanching = cfg.AllowAvalanching;
	SimpleLength = cfg.SimpleLength;
	SimpleLengthFactor = cfg.SimpleLengthFactor;
	numStab = cfg.numStab;
	Hifactor = cfg.Hifactor;
	Hcrit_global = cfg.Hcrit_global;
	transport_eq = cfg.transport_eq;
	alpha_varies = cfg.alpha_varies;
	alpha_lag = cfg.alpha_lag;
	moeilijkdoen = cfg.moeilijkdoen;
	correction_NT = cfg.correction_NT;
	Npsl_min = cfg.Npsl_min;
	stle_factor = cfg.stle_factor;

	q_in1 = cfg.q_in1;
	H0 = cfg.H0;
	ii = cfg.ii;
	D50 = cfg.D50;
	thetacr = cfg.thetacr;
	dts = cfg.dts;
	nd = cfg.nd;
	readbed = cfg.readbed;
	readfw = cfg.readfw;

	sepcritangle = cfg.sepcritangle * grad_2_deg;
	g = cfg.g;
	F = g * ii;
	kappa = cfg.kappa;
	tt = cfg.tt;
	tresh = cfg.thresh;
	max_it = cfg.max_it;

	denswater = cfg.denswater;
	sgsand = denssand / denswater;
	delta = sgsand - 1;
	ampbeds = ampbeds_factor * D50;
	epsilonp = cfg.epsilonp;
	repose = cfg.repose * grad_2_deg;
	m = cfg.m;
	alpha = m / (delta * g);
	be = cfg.be;
	l1 = 1.9 * D50;
	//l2 = 1 / tan(-repose);
	l2 = 1.73;
	F0 = cfg.F0;
	F0_dim = correction_NT * F0 * sqrt(g * delta / D50);
	meanstle = alpha_lag * D50;
	A2_geom = cfg.A2_geom;
	A3_geom = cfg.A3_geom;
	k2 = cfg.k2;

	alpha_2 = cfg.alpha_2;
	D_star  = D50 * cbrt(g * delta / (nu * nu));
	w_s = nu / D50 * (pow(pow(10.36, 2) + 1.049 * pow(D_star, 3),(1./2.)) - 10.36);
	u_star_cr = sqrt(thetacr * g * delta * D50);
	alpha_min_SK = cfg.alpha_min_SK;
	alpha_max_SK = cfg.alpha_max_SK;

	alpha_min_S = cfg.alpha_min_S;
	alpha_max_S = cfg.alpha_max_S;
	theta_min_S = cfg.theta_min_S;
	theta_max_S = cfg.theta_max_S;
	H_ref = cfg.H_ref;
	keepsgrowing = cfg.keepsgrowing;
}

double maxval(vec vinp) {
	int len = vinp.size();
	double nm = vinp[0];
	for (int i = 1; i < len; i++) {
		double vi = vinp[i];
		if (nm < vi) nm = vi;
	}
	return nm;
}

void doStabAnalysis(int stabWrite, flow& H2O, bottom& sand, double& q_in){
	int num=numStab; int cols=4;
	// JW vector<vector<double> > dta(num+1,cols);
	vector<vector<double> > dta(num+1,vector<double>(cols));
	
	//double Lmax=Hi*Hifactor; //OLAV: changed 2011 04 01 (was 10*Hi)
	double minfactor = 3.; 
	double Lmax=H*Hifactor; //OLAV: changed 2014 01 31
	double Lstep=(Lmax-H*minfactor)/num; //was double Lstep=Lmax/num-H*3; 
	//OLAV: 2012 09 17: shouldnt this be H? 
	
	dt=dts;
	vec bedstab(Npx,0.0);
	vec newbed(Npx,0.0);
	vec ubed(Npx,0.0);
	vec dump1(Npx,0.); vec dump2(Npx,0.); vec dump3(Npx,0.);
	for(int i=0;i<Npx;i++){
		bedstab[i]=ampbeds*sin(1*2.0*M_PI/Npx*(i));
	} 
	for (int p=0;p<=num;p++){
		//L=Hi/10+Lstep*(p+1);  //2013 1 31: OLAV (was L=Hi/10+Lstep*(p);)
		//L=Hi/numStab+Lstep*(p); //2012 09 17: OLAV (was with /10., now with numStab)
		//L=Hi*5+Lstep*(p); //2012 09 17: OLAV test
		L=H*minfactor+Lstep*(p);  //OLAV: changed 2014 01 31 was L=Hi*5+Lstep*(p);
		dx=L/Npx;

		cerr<<p<<" "<<L<<" "<<H<<" "<<dx<<endl;
		setS_Av();
		H2O.resetIu();
		H2O.solve(bedstab);
		if (p==0) {
			doCheckQsp(bedstab, H2O, q_in);
			cerr<<"H = "<<H<<"m"<<endl;
		}
		H2O.u_b(ubed);
		newbed=sand.update(ubed,dump1,dump2,dump3);
		double gri=(1/dt)*log(maxval(newbed)/ampbeds);;
		cerr<<"gri: "<<gri<<endl;
		double mig=sand.detMigr(bedstab,newbed);
		cerr<<"mig: "<<mig<<endl<<endl;
		
		//cerr<<"ik kom hier p= " << p << endl;  //OLAV 2011 2 22 TEST
		dta[p][0]=L; dta[p][1]=H; dta[p][2]=gri; dta[p][3]=mig;
	}
	
	//cerr<<"ik kom hier 1" << endl;  //OLAV 2011 2 22 TEST
	
	if (stabWrite==1) {
		ofstream outstab("out_stab.out");
		outstab.precision(16);
		for (int j=0;j<=num;j++) outstab<<dta[j][0]<<" "<<dta[j][1]<<" "<<dta[j][2]<<" "<<dta[j][3]<<endl;
		outstab.close();
		}
		else {
	    ofstream outstab("out_stab_during.out");
		outstab.precision(16);
		for (int j=0;j<=num;j++) outstab<<dta[j][0]<<" "<<dta[j][1]<<" "<<dta[j][2]<<" "<<dta[j][3]<<endl;
		outstab.close();
	}
	
    //cerr<<"ik kom hier 2" << endl;  //OLAV 2011 2 22 TEST
	
	double gr = -999;
	int row = 0;
	int iinit=0;
	
	//cerr<<"ik kom hier 3" << endl;  //OLAV 2011 2 22 TEST
	
	//OLAV 2012 09 17: commented this part away
	// if (tijd==0 && readbed1==0){ // OLAV 2011 05 12
                    
	// while (dta[iinit][2]>0) { 
          // //cerr<<"ik kom hier dta[iinit][2]= " << dta[iinit][2] << endl;
          // iinit++;
          // }
    // } // OLAV 2011 05 12
	//cerr<<"ik kom hier 4" << endl;  //OLAV 2011 2 22 TEST
	
	for(int i=iinit;i<=num;i++){
			if (dta[i][2]>gr){
					gr=dta[i][2];
					row=i;
			}
	}	
	
	L=dta[row][0];
	H=dta[row][1];
	dz=H/Npz;
	dx=L/Npx;
	setS_Av();
	dt=dtr;
	cerr<<"Finished the stability analysis"<<endl<<endl;
	cerr<<"L = "<<L<<endl;
	cerr<<"H = "<<H<<endl;
}

void doCheckQsp(vec bedflow, flow& H2O, double& q_in){
	//checken van de specifieke afvoer
	double q_sp=H2O.check_qsp();
	cerr<<"check of specific discharge: "<<q_sp<<endl;
	double q_tol=q_in/100.;//1e-3; HIER KUN JE DE WATEROPPLAATSING UITZETTEN!!!!!!!!
	//double q_tol=1e-3;//1e-3; HIER KUN JE DE WATEROPPLAATSING UITZETTEN!!!!!!!!
	double q_dif=q_sp-q_in;
	int q_tel=0;
	cerr<<"q_dif = " <<q_dif<<" (q_in = "<<q_in<<")"<<endl;
	while (fabs(q_dif)>q_tol){
		if (q_tel==0) {
			outlog<<"T="<<tijd<<" - WARNING: water depth changed to ensure constant discharge."<<endl;
			q_tel++;
		}
		double q_cor=fabs(q_dif/(H*10.));
		if (q_dif<0.) H+=q_cor;
		else if (q_dif>0.) H-=q_cor;
		dz=H/double(Npz);
		setS_Av();
		//n_it_fl=H2O.solve_gm(sand.getShape(sepflag),20);
		int n_it_fl=0;
		n_it_fl=H2O.solve_gm(bedflow,20);
		if (n_it_fl==-1) H2O.solve(bedflow);
		q_sp=H2O.check_qsp();
		q_dif=q_sp-q_in;
		cerr<<"q_dif = " <<q_dif<<" (q_in = "<<q_in<<")"<<endl;
		cerr<<"check of specific discharge: "<<endl<<q_sp<<endl;
	}
}

void setS_Av(){
	const auto ustar=pow(g*H*ii,0.5);
	S=BETA1*ustar;
	Av=BETA2*(1./6.)*kappa*H*ustar;
}
