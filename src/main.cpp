// main.cpp

#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include "bottom.h"
#include "flow.h"
#include "vecmat.h"
#include "admin.h"
#include "Config.h"
#include "BedConfig.h"
#include "FlowConfig.h"
#include "Logging.h"
#include <ctime>
#include <cstring>
#include <cassert>

using namespace std;

//ofstream outlog("out_log.txt");

double H;
double dt;
double L;
double dx;
double dz;
double tijd;
vec Avx;
double S;

void doStabAnalysis(flow& H2O, bottom& sand, const double& q_in, const Config& cfg);
void doCheckQsp(vec bedflow, flow& H2O, const bottom& sand, const double& q_in, const Config& cfg);
void setS_Av(const Config& cfg, const bottom& sand);


int main (int argc, char * const argv[]) {

cerr.precision(16);

std::string filename = (argc == 1) ? "config.cfg" : argv[1];
const Config cfg(filename);
const BedConfig bedConfig(cfg);
const FlowConfig flowConfig(cfg);

dude_log::init(cfg.FileName, cfg.FileLevel, cfg.ConsoleLevel);

admin::Npx = cfg.Npx; // still necessary for admin::o2()

// Initialize global variables
H = cfg.H0;
dt = cfg.dtr;
L = cfg.Lini;//1.0;
dx = L / cfg.Npx;
dz = H / cfg.Npz;
tijd = 0.0;
S = 0.0;

double q_in = cfg.q_in1;

flow H2O(flowConfig);
bottom sand(bedConfig);
Avx.resize(cfg.Npx);

//if (cfg.dt_write==1.) {cerr<<endl<<endl<<endl<<endl<<endl<<"         ------ NOTE!! DT_WRITE==1!! -------"<<endl<<endl<<endl<<endl<<endl;}
if (cfg.dt_write==1.)
	DUDE_LOG(warning) << SHOW_VAR(cfg.dt_write);

for (int p=1;p<=1;p++){				//superloop!!!!!!!!!!!!

	//cerr<<flowConfig.F<<endl;
	DUDE_LOG(info) << SHOW_VAR(flowConfig.F);
	const auto ampbeds = cfg.ampbeds_factor * cfg.D50;
	sand.setSin(ampbeds,1);
//	sand.setRand(ampbeds,28); //1.1*cfg.D50,(unsigned)time(0),
	//sand.setRand(1e-8,(unsigned)time(0));
	//sand.setShape(current);
	setS_Av(cfg, sand);
	H2O.resetIu();
	int iinit1=0;
	//vec inp(3,0.0); //is input if a bottom profile is read
	
	// NEW STUFF, WRITTEN BY OLAV 2014 for FLOODWAVE
	vec fw_t(2,0.0);
	vec fw_q(2,0.0);
	if (!cfg.readfw.empty()){
		double temp;
		ifstream in1(cfg.readfw);
		if (!in1) {
			DUDE_LOG(error) << "Can't open floodwave file: " << cfg.readfw;
			//std::cerr << "ERROR: Can't open floodwave file: " << cfg.readfw << std::endl;
			std::exit(1);
		}
		in1>>temp;
		int np = int(temp);
		fw_t.resize(np);
		fw_q.resize(np);
		DUDE_LOG(info) <<"Number of points in floodwave: " << np;
		//cerr<< "Number of points in floodwave: " << np<<endl;
			for (int j = 0; j < np; j++) {
				in1 >> fw_t[j] >> fw_q[j];
				DUDE_LOG(info) << SHOW_VAR(fw_t[j]) << SHOW_VAR(fw_q[j]);
				//cerr << fw_t[j] << " " << fw_q[j] << endl;
			}
	}// NEW STUFF, WRITTEN BY OLAV 2014 for FLOODWAVE
	
	if (!cfg.readbed.empty()) {
		DUDE_LOG(info) << " Bed elevation read from file " << cfg.readbed;
		//cerr<<endl<<endl<<endl
		//		<<"------>> NOTE: Bed elevation read from file " << cfg.readbed << " <<-------"
		//		<<endl<<endl<<endl;
		auto inp=sand.readBottomInp(cfg.readbed);
		H=inp[1];
				
		dz=H/cfg.Npz;
		tijd=inp[0];
		
		if (!cfg.readfw.empty()) {
			q_in=interpolate(fw_t,fw_q,tijd); // NEW STUFF, WRITTEN BY OLAV 2014 for FLOODWAVE
		}
		iinit1=int(tijd/dt);
		L=inp[2];
						
		dx=L/cfg.Npx;
		setS_Av(cfg, sand);
		DUDE_LOG(info) <<"read check:";
		DUDE_LOG(info) <<"Timeprevious: " <<tijd;
		DUDE_LOG(info) <<"H: "<<H;
		DUDE_LOG(info) <<"L: "<<L;
		DUDE_LOG(info) <<"first bed point: "<<sand.getShape(0)[0];
		DUDE_LOG(info) <<"last bed point: "<<sand.getShape(0)[cfg.Npx-1];
		DUDE_LOG(info) <<"dx: "<<dx;
		//cerr<<"read check:"<<endl;
		//cerr<<"Timeprevious: " <<tijd<<endl;
		//cerr<<"H: "<<H<<endl;
		//cerr<<"L: "<<L<<endl;
		//cerr<<"first bed point: "<<sand.getShape(0)[0]<<endl;
		//cerr<<"last bed point: "<<sand.getShape(0)[cfg.Npx-1]<<endl;
		//cerr<<"dx: "<<dx<<endl;
		//cerr<<"Av: "<<Av<<endl;
		//cerr<<"S : "<<S<<endl;
		
		//avalanching protocol
		// if(transport_eq == 2 && AllowAvalanching == 1){
		// avalanche(); 
		// }
	}
	else {
		if (!cfg.readfw.empty()){
			q_in=interpolate(fw_t,fw_q,tijd);
		}
	}// NEW STUFF, WRITTEN BY OLAV 2014 for FLOODWAVE
	
	DUDE_LOG(info) << "first discharge: " << q_in;

	ostringstream tmpbot;
	tmpbot << "out_bottom"<< p << ".txt";
	string ofname = tmpbot.str();
	ofstream outbot(ofname.c_str(),ios_base::out);

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


	auto current=sand.getShape(0);

	int write_teller = int(cfg.dt_write/dt);
	int cor=1; if (cfg.dt_write==dt) cor=0;
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

	time_t now;
	char buf[80];
	std::time(&now);
	std::strftime(buf, sizeof(buf), "%a %Y-%b-%d %H:%M:%S", std::localtime(&now));
	//outlog<<"Simulation started at "<<buf<<endl;
	DUDE_LOG(info) << "Simulation started at " << buf;

	auto myH = cfg.H0;
	auto oldL = L;
	auto Lstab = L;
	for(int i=iinit1;i<=cfg.tend/dt;i++){ // 40 minutes
		
		if (!cfg.readfw.empty()) {
			q_in=interpolate(fw_t,fw_q,tijd); // NEW STUFF, WRITTEN BY OLAV 2014
		}
		//OLAV: 2011 02 21 changed from 
        //Hdiff=(H-H0)/H0;

		// TODO: Hdiff moet eigenlijk het verschil in waterdiepte met de laatste keer dat stabanalysis is gedaan zijn, en niet de vergelijkin met de beginwaarde

		const auto Hdiff=abs((H-myH)/myH); //equals 0 when exactly the same, 1 when the difference is 100%
		const auto Hcrit = cfg.Hcrit_global;
		//cerr<<"Hdiff = " <<Hdiff <<" (Hcrit = "<<Hcrit<<")"<<endl;
		//DUDE_LOG(info) << SHOW_VAR(Hdiff) << SHOW_VAR(Hcrit);
		DUDE_LOG(info) << SHOW_2VARS(Hdiff, Hcrit);

		bool doStab=0;
		if (Hdiff>=Hcrit ) doStab=1; //OLAV: 2011 02 21 changed from 
        //if (Hdiff>=1.+Hcrit) doStab=1;
        
		auto updateMyH = false;
		/* start with initial stability analysis, or if H is sufficiently changed */
		//if ( (i==iinit1&&readbed1==0) | doStab==1 | i==1000) { // OLAV TEST 2011 2 23
		if (cfg.SimpleLength==0) { // OLAV 2012 09 06: added simple length implementation
		   if ( (i==iinit1&&cfg.readbed.empty()) || doStab==1) {
			  //outlog<<"T="<<tijd<<" - WARNING: Stability Analysis. (Hdiff="<<Hdiff<<")"<<endl;
			  DUDE_LOG(warning) << "Stability Analysis: " SHOW_VAR(Hdiff);
			  oldL = L;
			  doStabAnalysis(H2O, sand, q_in, cfg);
			  Lstab = L;
			  if ((i==iinit1&&cfg.readbed.empty())){
				  oldL =L;
			  }
			  updateMyH = true;
			  dt=cfg.dtr; // reset, stab analysis uses dt=dts
//			  const auto Hstab = H;
			  //OLAV: 2013 02 06 added doStab=0;
			  doStab=0; 
		   }
		   auto doLag= false;
		   if (!doLag ||(i == iinit1 && cfg.readbed.empty())) {
			   L = Lstab;
		   } else if (std::abs(oldL / Lstab - 1) > 0.01) {
			   L = oldL + 0.01 * (Lstab - oldL);
		   }
		   dx = L / cfg.Npx;
         }
        else if (cfg.SimpleLength==1) {
			//L=1.6;
            L = H*cfg.SimpleLengthFactor;
			//L=319.;
			//L=4.;
			
            dz=H/cfg.Npz;
	        dx=L/cfg.Npx;
			//cerr << "dx: " << dx << endl;
	        setS_Av(cfg, sand);
	        dt=cfg.dtr;
			doStab=0;
        }
        else if (cfg.SimpleLength == 2) {//don't change L
        	dz=H/cfg.Npz;
			dx=L/cfg.Npx;
			setS_Av(cfg, sand);
			dt=cfg.dtr;
			doStab=0;
        }
		assert(doStab==0);

		//q_in=qa[i];
	  	
	  	//2012 09 13 Olav 
		if(cfg.AllowFlowSep == 1){
			sand.checkFlowsep(); 
		}

		//2012 09 13 Olav
		const auto& stateFsz=sand.getFsz();
		const auto nf = stateFsz.size();
		int solve_method=stateFsz[nf-3];
		int sepflag=stateFsz[nf-2];
		int nfsz=stateFsz[nf-1];//}

		current=sand.getShape(0); // for determination of migration rate
		bedflow=sand.getShape(sepflag);

		int n_it_fl=0;
		if (i==iinit1) {
			//cerr<<"Eerste keer stroming oplossing met SOLVE"<<endl;
			//outlog<<"Eerste keer stroming oplossing met SOLVE"<<endl;
			DUDE_LOG(info) << "Eerste keer stroming oplossing met SOLVE";
			n_it_fl=H2O.solve(bedflow);
			//cerr<<"number of flow iteration required: "<<n_it_fl<<endl;
			DUDE_LOG(info) << "number of flow iteration required: " << n_it_fl;
		} else if (doStab==1) {
			//cerr<<"Stroming oplossing met SOLVE vanwege Stability Analysis"<<endl;
			//outlog<<"T="<<tijd<<" - WARNING: SOLVE vanwege Stability Analysis. (Hdiff="<<Hdiff<<")"<<endl;
			DUDE_LOG(info) << "Stroming oplossing met SOLVE vanwege Stability Analysis: " << SHOW_VAR(Hdiff);
			H2O.resetIu();
			n_it_fl=H2O.solve(bedflow);
			//cerr<<"number of flow iteration required: "<<n_it_fl<<endl;
			DUDE_LOG(info) << "number of flow iteration required: " << n_it_fl;
		} else if (solve_method>=1) {
			//cerr<<"Stroming oplossing met SOLVE vanwege sterke 'bodem' veranderingen"<<endl;
			//cerr<<"solve_method= "<<solve_method<<endl;
			DUDE_LOG(info) << "Stroming oplossing met SOLVE vanwege sterke 'bodem' veranderingen: " << SHOW_VAR(solve_method);
			H2O.resetIu();
			n_it_fl=H2O.solve(bedflow);
			//cerr<<"number of flow iteration required: "<<n_it_fl<<endl;
			DUDE_LOG(info) << "number of flow iteration required: " << n_it_fl;
		} else {
			n_it_fl=H2O.solve_gm(bedflow,20);
		}
		if (n_it_fl == -1) {
			ofstream outtemp("out_temp.txt");
			auto bedtemp = sand.getShape(0);
			for (auto temp : bedtemp)
				outtemp << temp << " ";
			outtemp << endl;
			bedtemp = sand.getShape(1);
			for (auto temp : bedtemp)
				outtemp << temp << " ";
			outtemp << endl;
			outtemp.close();
			H2O.resetIu();
			H2O.solve(bedflow);
		}

		doCheckQsp(bedflow, H2O, sand, q_in, cfg);
		vec u0_b(cfg.Npx);
		vec U0_mean(cfg.Npx);
		H2O.u_b(u0_b);
		H2O.Umean(sand.getShape(sepflag), U0_mean);
		if (updateMyH)
			myH = H;
		
		//cerr << "current: " << current[0] << " bedflow: " << bedflow[0] << endl; //OLAV 2014 03 31

		if (write_teller == cfg.dt_write / dt) {
			if (cfg.write_velocities)
				H2O.write_velocities(tijd, sand.getShape(sepflag), u0_b);
			//H2O.write_zeta(tijd);
			const auto &stateZeta = H2O.getZeta();
			outzeta << tijd << " ";
			for (auto zeta : stateZeta)
				outzeta << zeta << " ";
			outzeta << endl;
			//2014 01 27: changed q_sp to q_in
			outbot << tijd << " " << sepflag << " " << nfsz << " " << q_in << " " << H << " " << L << " ";
			for (auto b : current)
				outbot << b << " ";
			outbot << endl;
			vec bot = sand.getShape(sepflag);
			//2014 01 27: changed q_sp to q_in
			outbot << tijd + 0.1 << " " << sepflag << " " << nfsz << " " << q_in << " " << H << " " << L << " ";
			for (auto b : bot)
				outbot << b << " ";
			outbot << endl;
		}

		//cerr << "bint1: " << bint1 << " bint2: " << bint2 << endl; //OLAV 2014 03 31

		auto bint1 = sand.detint1(current);
		auto bint2 = sand.detint2(current);
		auto zetaint1 = H2O.zetaint1();
		auto zetaint2 = H2O.zetaint2();
		const auto& Dc = sand.detNd_fft(sand.getShape(0),2); // dune characteristics
		auto Nd = int(Dc[0]);
		double cr = Dc[1];
		double tr = Dc[2];
		double Hav = cr-tr;
		double Lav = L;
		
		//cerr << "Hav: " << Hav << endl; //OLAV 2014 03 31

		vec next;
		vec bss1(cfg.Npx);
		vec bss2(cfg.Npx);
		vec fluxtot(cfg.Npx);
		vec dhdx(cfg.Npx);
		if (cfg.transport_eq==0 ||cfg.transport_eq==1 || cfg.transport_eq==3 ){
			if (sepflag==0){
				next=sand.update(u0_b,U0_mean,bss1,fluxtot,dhdx);
				bss2=bss1;}
			else if (sepflag==1){
				next=sand.update_flowsep(u0_b,U0_mean,bss1,bss2,fluxtot,dhdx);}
			}
		else if (cfg.transport_eq==2) {
			if (sepflag==0){
				next=sand.update(u0_b,U0_mean,bss1,fluxtot,dhdx);
				bss2=bss1;}
			else if (sepflag==1){
				next=sand.update_flowsep(u0_b,U0_mean,bss1,bss2,fluxtot,dhdx);}
		}
		
		//cerr << "next: " << next[1] << " bint2: " << bint2 << endl; //OLAV 2014 03 31
		
		double flux_av=0.;
		for(auto f : fluxtot) flux_av+=f;
		flux_av/=cfg.Npx;

		double mig=sand.detMigr(current,next);
		//outint<<tijd<<" "<<H<<" "<<bint1<<" "<<bint2<<" "<<zetaint1<<" "<<zetaint2<<endl;

	  	if(write_teller==cfg.dt_write/dt-cor){
			const auto& stateFsz=sand.getFsz();
			const auto& stateSr=sand.getSr();
			//wegschrijven fsz:
			outfsz<<tijd<<" "; for (auto fsz : stateFsz) outfsz << fsz << " "; outfsz<<endl;
			outSround<<tijd<<" "; for(auto sr : stateSr) outSround << sr << " "; outSround<<endl;
		}

	  	if(write_teller==cfg.dt_write/dt){
			outint<<tijd<<" "<<H<<" "<<bint1<<" "<<bint2<<" "<<zetaint1<<" "<<zetaint2<<" "<<cr<<" "<<tr<<" "<<Nd<<" "<<mig<<" "<<flux_av<<endl;
			outflux<<tijd<<" ";    for(auto f : fluxtot) outflux << f << " "; outflux <<endl;
			outdhdx<<tijd<<" ";    for(auto d : dhdx) outdhdx << d <<" "; outdhdx <<endl;
			outbss<<tijd<<" ";     for(auto b : bss1) outbss << b <<" "; outbss <<endl;
			outbss<<tijd+0.1<<" "; for(auto b : bss2) outbss << b <<" "; outbss <<endl;
			write_teller = 0;
		}

	  	const auto nt = next.size();
	  	vec verschil(nt);
	  	for(auto k = 0u; k < nt ; k++) verschil[k] = next[k] - current[k];
	  	auto norm=L2(verschil);
	  	sand.setShape(next);

	  	current=next;
	  	tijd+=dt;

	  	DUDE_LOG(info) << "flowsolver " << tijd << " seconden (" << tijd / 60 <<" min)" << " onderweg";
	  	DUDE_LOG(info) << "number of flow iteration required: " << n_it_fl;
	  	DUDE_LOG(info) << SHOW_VAR(sepflag) << SHOW_VAR(nfsz) << SHOW_VAR(Nd)
	  			<< SHOW_VAR(H) << SHOW_VAR(Lav) << SHOW_VAR(Lav);
	  	DUDE_LOG(info) << "integral of bed: " << SHOW_VAR(bint1);
	  	DUDE_LOG(info) << "bodem ge update met (L2) : " << norm << " tot (L2) : " << L2(next);
	  	//cerr<<"flowsolver "<<tijd<<" seconden ("<<tijd/60.<<" min)"<<"	onderweg"<<endl;
	  	//cerr<<"number of flow iteration required: "<<n_it_fl<<endl;
	  	//cerr<<"sepflag: "<<sepflag<<"; nfsz: "<<nfsz<<endl;
	  	//cerr<<"Nd: "<<Nd;
	  	//cerr<<"; wd: "<<H;
	  	//cerr<<"; Lav: "<<Lav;
	  	//cerr<<"; Hav: "<<Hav<<endl;
	  	//cerr<<"integral of bed: "<<bint1<<endl;
	  	//cerr<<"bodem ge update met (L2) : "<<norm<<" tot (L2) : "<<L2(next)<<endl<<endl;

	  	write_teller+=1;
	  	if(Hav<.5*ampbeds) {
#if 0
	  		cerr<<"Dune height very low: " << Hav <<" , bailing out."<<endl<<endl;
	  		break;
#else
	  		sand.setSin(ampbeds,1);
	  		DUDE_LOG(warning) << "Hav very low, bed set to initial disturbance.";
	  		//cerr<<"Hav very low, bed set to initial disturbance."<<endl<<endl;
#endif
	  	}
	  	setS_Av(cfg, sand);
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
	std::strftime(buf, sizeof(buf), "%a %Y-%b-%d %H:%M:%S", std::localtime(&now));
	DUDE_LOG(info) << "Simulation ended at " << buf;
	//outlog<<"Simulation ended at "<<buf<<endl;

#if 0
	outlog.close();

	char ofname9[15]="out_log";
	char pstr [10];
	sprintf(pstr,"%d",p);
	strcat ( ofname9, pstr );
	strcat ( ofname9, ".txt" );
	//cerr<<ofname9<<endl;

  char oldname[] ="out_log.txt";
  rename( oldname , ofname9 );
	outlog.open("out_log.txt");
  outlog<<"run initialized; this file may be safely deleted"<<endl;
#endif

}
	return 0;
}

double maxval(const vec& vinp) {
	auto nm = vinp[0];
	for (auto v : vinp)
		nm = max(v, nm);
	return nm;
}

void doStabAnalysis(flow& H2O, bottom& sand, const double& q_in, const Config& cfg){
	int num=cfg.numStab; int cols=4;
	// JW vector<vector<double> > dta(num+1,cols);
	vector<vector<double> > dta(num+1,vector<double>(cols));
	double Lmin;
	double Lmax;
	if (cfg.Lrangefix){
		Lmin=cfg.Minfactor;
		Lmax=cfg.Hifactor;
	}
	else {
		Lmin=H*cfg.Minfactor;
		Lmax=H*cfg.Hifactor;
	}
	double Lstep= (Lmax-Lmin)/num;
	
	dt=cfg.dts;
	vec bedstab(cfg.Npx,0.0);
	vec newbed(cfg.Npx,0.0);
	vec ubed(cfg.Npx,0.0);
	vec U_mean(cfg.Npx,0.0);
	vec dump1(cfg.Npx,0.); vec dump2(cfg.Npx,0.); vec dump3(cfg.Npx,0.);
	const auto ampbeds = cfg.ampbeds_factor * cfg.D50;
	for(int i=0;i<cfg.Npx;i++){
		bedstab[i]=ampbeds*sin(1*2.0*M_PI/cfg.Npx*(i));
	} 
	for (int p=0;p<=num;p++){
		//L=Hi/10+Lstep*(p+1);  //2013 1 31: OLAV (was L=Hi/10+Lstep*(p);)
		//L=Hi/numStab+Lstep*(p); //2012 09 17: OLAV (was with /10., now with numStab)
		//L=Hi*5+Lstep*(p); //2012 09 17: OLAV test
		if (cfg.Lrangefix){
			L=Lmin+Lstep*(p);
		}
		else {
			L=H*cfg.Minfactor+Lstep*(p);
		}
		dx=L/cfg.Npx;

		setS_Av(cfg, sand);
		H2O.resetIu();
		H2O.solve(bedstab);
		if (p==0) {
			doCheckQsp(bedstab, H2O, sand, q_in, cfg);
			DUDE_LOG(info) << "Stab Analysys starts: " << SHOW_VAR(H);
		}
		H2O.u_b(ubed);
		H2O.Umean(bedstab, U_mean);//Umean over the bed in the linstab
#if 0
		auto& mySand = sand;
#else
		const BedConfig bcfg(cfg);
		bottom mySand(bcfg);
		mySand.setSin(ampbeds, 1);
#endif
		newbed=mySand.update(ubed,U_mean,dump1,dump2,dump3);
		//TODO LL: Kijken welke modus groeit (fourier analyse)
		double gri=(1/dt)*log(maxval(newbed)/ampbeds);
		//TODO LL: double gri = sand.detGrow()
		double mig=mySand.detMigr(bedstab,newbed);

		DUDE_LOG(debug) << SHOW_VAR(p) << SHOW_VAR(L) << SHOW_VAR(gri);
		dta[p][0]=L; dta[p][1]=H; dta[p][2]=gri; dta[p][3]=mig;
	}
	
	double gr = -999;
	int row = 0;
	int iinit=0;

	for(int i=iinit;i<=num;i++){
			if (dta[i][2]>gr){
					gr=dta[i][2];
					row=i;
			}
	}

	if (row>num-1){
		DUDE_LOG(warning) << "No real maximum found in Linstab; increase Hifactor";
	}
	else if (row<1){
		DUDE_LOG(warning) << "No real maximum found in Linstab; decrease Minfactor";
	}

	L=dta[row][0];
	H=dta[row][1];
	dz=H/cfg.Npz;
	dx=L/cfg.Npx;
	setS_Av(cfg, sand);
	DUDE_LOG(info) << "Finished the stability analysis";
	DUDE_LOG(info) << SHOW_VAR(L) << SHOW_VAR(H);

	static int seq;
	seq++;
	ostringstream oss;
	oss << "out_stab" << seq << ".out";
	ofstream outstab(oss.str());
	//outstab.precision(16);
	const auto w = 15;
	if (row == 0 || row > num-1)
		outstab << "# !! WARNING no real maximum";
	outstab << endl;
	outstab << "# t=" << tijd << " row=" << row << " L=" << L << " grow=" << gr << endl << endl;
	outstab << "#" << setw(w-1) << "L" << setw(w) << "H" << setw(w) << "grow" << setw(w) << "migr" << endl;
	for (int j = 0; j <= num; j++)
		outstab << setw(w) << dta[j][0] << setw(w) << dta[j][1]
				<< setw(w) << dta[j][2] << setw(w) << dta[j][3] << endl;
	outstab.close();
}

void doCheckQsp(vec bedflow, flow& H2O, const bottom& sand, const double& q_in, const Config& cfg){
	//checken van de specifieke afvoer
	double q_sp=H2O.check_qsp();
	DUDE_LOG(info) << "check of specific discharge: " << SHOW_VAR(q_sp);
	//cerr<<"check of specific discharge: "<<q_sp<<endl;
	double q_tol=q_in/100.;//1e-3; HIER KUN JE DE WATEROPPLAATSING UITZETTEN!!!!!!!!
	//double q_tol=1e-3;//1e-3; HIER KUN JE DE WATEROPPLAATSING UITZETTEN!!!!!!!!
	double q_dif=q_sp-q_in;
	int q_tel=0;
	//cerr<<"q_dif = " <<q_dif<<" (q_in = "<<q_in<<")"<<endl;
	while (fabs(q_dif)>q_tol){
		if (q_tel==0) {
			DUDE_LOG(warning) << "water depth changed to ensure constant discharge";
			//outlog<<"T="<<tijd<<" - WARNING: water depth changed to ensure constant discharge."<<endl;
		}
		q_tel++;
		double q_cor=fabs(q_dif/(H*10.));
		if (q_dif<0.) H+=q_cor;
		else if (q_dif>0.) H-=q_cor;
		dz=H/double(cfg.Npz);
		setS_Av(cfg, sand);
		//n_it_fl=H2O.solve_gm(sand.getShape(sepflag),20);
		int n_it_fl=0;
		n_it_fl=H2O.solve_gm(bedflow,20);
		if (n_it_fl==-1) H2O.solve(bedflow);
		q_sp=H2O.check_qsp();
		q_dif=q_sp-q_in;
		//cerr<<"q_dif = " <<q_dif<<" (q_in = "<<q_in<<")"<<endl;
	}
	if (q_tel > 0) {
		DUDE_LOG(info) << q_tel << " rechecks of specific discharge: " << SHOW_VAR(q_dif);
		//cerr<<q_tel<<" rechecks of specific discharge: q_dif="<<q_dif<<endl;
	}
}

void setS_Av(const Config& cfg, const bottom& sand){
	const auto ustar = sqrt(cfg.g * H * cfg.ii);
	S = cfg.BETA1 * ustar;//0.01;//0.0001;
	const auto Av = cfg.BETA2 * (1./6.) * cfg.kappa * H * ustar;//0.004;//

	//DUDE_LOG(warning) << SHOW_VAR(Av);
	auto dhdx = sand.get_dhdx();
	for (auto i = 0 ; i < cfg.Npx; i++) {
		Avx[i] = Av;// * (1 + dhdx[i]);
//		std::cout << std::setprecision(4) << Avx[i] << " ";
	}
	//std::cout << std::endl;
}
