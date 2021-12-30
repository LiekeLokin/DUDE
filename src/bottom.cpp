// bottom.cpp

#include "bottom.h"
#include "admin.h"
#include "fft.h"
#include <iomanip>
#include <fstream>

using namespace std;
using namespace admin;

namespace {
const double grad_2_deg = 360./(2.*M_PI);
	// help variable for grad to degress
}
bottom::bottom(const BedConfig& cfg) : cfg(cfg), Npx(cfg.Npx), nf((int(ceil(Npx/10+1)))*7+2), nf2(int(ceil(Npx/10+1))) {
	b=new vec(Npx,0.0);
	bp=new vec(Npx,0.0);
	flux=new vec(Npx,0.0);
	x=new vec(Npx,0.0);
	fsz=new vector<int>(nf,0);
	Sr=new vec(nf2,0.0);
}
bottom::~bottom(){
	delete b;
	delete bp;
	delete flux;
	delete x;
	delete fsz;
	delete Sr;
}

int bottom::o3(int i_in) const {
	/*adres vertaling voor periodieke rvw ten behoeve van parameterisatie*/
	int i_uit = i_in;

	if (i_in>=Npx) i_uit = i_in-Npx;
	if (i_in<0) i_uit = i_in+Npx;

	if (i_uit>=0 && i_uit<Npx)
	 	{return i_uit;} // alleen dan return (want dan 0<=i<Npx)
	 					// zoniet print error
	else {
		cerr<<"  ERROR: o3() buffer overflow (i_in="<<i_in<<"; i_uit="<<i_uit<<")"<<endl;
		cout<<"T="<<tijd<<" - ERROR: o3() buffer overflow (i_in="<<i_in<<"; i_uit="<<i_uit<<")"<<endl;
        return i_uit;
	}
}

/*
======================================================
======================================================

BLOCK I: BOTTOM SET ROUTINES
in this block are functions to set up the initial or
temporal bottom profiles

======================================================
======================================================
*/

void bottom::setShape(vec b_in){
	(*b)=b_in;
}

void bottom::setSin(double amp){
	/*vullen van vector met voor de bodemverstoring geschaalde dz*/
	for(int i=0;i<Npx;i++){
		(*b)[i]=amp*sin(2.0*M_PI/Npx*(i));
	}
}
void bottom::setSin(double amp,int n){
	/*vullen van vector met voor de bodemverstoring geschaalde dz*/
	for(int i=0;i<Npx;i++){
		(*b)[i]=amp*sin(n*2.0*M_PI/Npx*(i));
	}
}

#if 0
void bottom::setDistSin(double amp,int n){
	/*vullen van vector met voor de bodemverstoring geschaalde dz*/
	for(int i=0;i<Npx;i++){
		(*b)[i]+=amp*sin(n*2.0*M_PI/Npx*(i));
	}
}


void bottom::setRand(double amp){
	/*fill bottom vector with a random disturbance, without seed*/
		for(int i=0;i<Npx;i++){
		(*b)[i]=amp*2.0*(0.5-(double)rand()/(double)(RAND_MAX));
	}
}
#endif

void bottom::setWave(int xwi, int xcin){

	double H_wave=0.004;
	double deltax = dx;
	int Np=xcin-xwi;
	//cerr<<Np<<endl;
	Np=int(floor((Np+0.0001)/2.));
	//cerr<<Np<<endl;
	Np*=2;
	//cerr<<Np<<endl;
	double L_wave= Np*deltax;
	double alpha=30.; //graden

	double tanalpha = sqrt(3.)/3.;
	//int Npx_wave = int(round((L_wave/deltax)/2.)*2.)+1;
	int Npx_wave = int(round(L_wave/deltax)) +1;

	vec x_wave (Npx_wave,0.0);

	for(int i=0;i<(Npx_wave);i++){
			(x_wave)[i]=(round(((-1*L_wave/2.)+deltax*i)*100.))/100.;
		}

	double x_trough1 = H_wave/tanalpha;
	double x_trough = deltax * round(x_trough1/deltax); //round
	double x_crest = -1*x_trough;

	//redefine H_wave
	double H2 = x_trough*tanalpha;

	//for(int i=0;i<x_wave.size();i++)cerr<<x_wave[i]<<" "; cerr<<endl<<endl;
	cerr<<"Npx_wave="<<Npx_wave<<endl;
	/*
	cerr<<"xtrough="<<x_trough;
	cerr<<"xcrest="<<x_crest;
	cerr<<"H2="<<H2<<endl;
	*/

	outlog<<"T="<<tijd<<" -          Npx_wave: "<<Npx_wave<<endl;

	//slope of leeside
	double slope_lee = tanalpha;
	//slope of stoss side
	double slope_stoss = H2/((L_wave/2.)-x_trough);

	vec wavelet1 (Npx_wave,0.0);

	for(int i=0;i<(Npx_wave);i++){

		if (x_wave[i]<x_crest){
		(wavelet1)[i]=slope_stoss*(x_wave[i]+(L_wave/2.));
		}
		else if(x_wave[i]>x_trough){
		(wavelet1)[i] = slope_stoss*(x_wave[i]-(L_wave/2.));
		}
		else{
			(wavelet1)[i]=-1*slope_lee*x_wave[i];
		}

		//cerr<<x_wave[i]<<" - "<<wavelet1[i]<<endl;

	}

	for(int i=xwi;i<xwi+Npx_wave;i++){							//xwi = beginpositie van de wave

			(*b)[o3(i)]+=wavelet1[i-xwi];

}}

#if 0
void bottom::setDist(double amp_dist){

		vec rand_dist(Npx,0.0);
		srand((unsigned)time(0));
		for(int k=0;k<Npx;k++){
		rand_dist[k]=amp_dist*2.0*(0.5-(double)rand()/(double)(RAND_MAX));
		}
		double randint = 0.0;
		for(int k=0;k<Npx;k++){
			randint+=rand_dist[k];
		}
		//cerr<<randint<<endl;
		for(int k=0;k<Npx;k++){
			rand_dist[k]-=randint/Npx;
		}
		randint = 0.0;
		for(int k=0;k<Npx;k++){
			randint+=dx*rand_dist[k];
		}

		cerr<<"dit moet nul zijn na correctie: " <<randint<<endl;

		//for(int k=0;k<rand_dist.size();k++)cerr<<rand_dist[k]<<" "; cerr<<endl;


		for(int k=0;k<Npx;k++){
			(*b)[k]+=rand_dist[k];
		}
}

void bottom::setRand(double amp,int seed){
	/*fill bottom vector with a random disturbance, with seed*/
	srand(seed);
	for(int i=0;i<Npx;i++){
		(*b)[i]=amp*2.0*(0.5-(double)rand()/(double)(RAND_MAX));
	}
}
#endif

void bottom::setMidSin(double amp, double length){
	/*fill bottom vector with a sine in the middle
	length of the sine is length, domain length is L*/
	int N2=(int)(length/dx);
	int start=(int)((Npx/2.0)-(N2/2.0));
	for(int i=0;i<start;i++)(*b)[i]=0.0;
	for(int i=start;i<start+N2;i++){
		(*b)[i]=amp*sin(2.0*M_PI/N2*(i-start));
	}
	for(int i=start+N2;i<Npx;i++)(*b)[i]=0.0;
}

void bottom::setCustom(double amp, int n){
	/*vullen van vector met voor de bodemverstoring geschaalde dz*/
	/*custom gemaakt op 10/02/2006*/
	for(int i=0;i<Npx;i++){
		(*b)[i]=amp*sin(n*2.0*M_PI/Npx*(i));
	}
	vec func(Npx,1.0);
	double dist=6.0; double x0=1.0; double xx;
	for(int i=0;i<Npx/2;i++) {
		xx=i*dx;
		if (xx>=0 && xx<=x0) {
			func[i]=0.0; }
		if (xx>=x0 && xx<=dist) {
			double a=0.5/dist;
			func[i]=a*(xx-x0);}
	}
	for(int i=Npx/2;i<Npx;i++){
		func[i]=func[Npx-i];
	}
	for(int i=0;i<Npx;i++){
		(*b)[i]*=func[i];
	}
}

vec bottom::readBottomInp(const std::string& readbed){
	/*fill bottom vector from input file, after Sobek computation*/
	vec inp(3,0.0);
	double dump; int sep;	double tijd; double wd; double Lin;
	ifstream in1(readbed);
	if (!in1) {
		std::cerr << "ERROR: Can't open bed file: " << readbed << std::endl;
		std::exit(1);
	}
	
	// ofstream outdebug;
    // outdebug.open ("out_debug1.txt", ofstream::out | ofstream::app);
    // outdebug.precision(16);
	// outdebug << tijd << " ";
	
	in1>>tijd; 
	// outdebug << tijd << " ";
    in1>>sep; 
   	// outdebug << sep << " ";
    in1>>dump;
	in1>>dump; 
    // outdebug << dump << " ";
    // in1>>dump;	
	// outdebug << dump << " ";
    in1>>wd; 
   	// outdebug << wd << " ";
    in1>>Lin;
   	// outdebug << Lin << " ";
   	// outdebug << endl;
	
	//wd=0.5;
	//wd=0.361417147609724;
	//wd=0.5;
	//Lin=4.;
    	
	cerr<<tijd<<" "<<sep<<" "<<wd<<endl;
	for(int i=0;i<Npx;i++)in1>>(*b)[i];
	if (wd==0.) {
		tijd=0.;
		wd=H; }
	inp[0]=tijd;
	inp[1]=wd;
	inp[2]=Lin;

	if (sep==1) {
		in1>>dump;
		for(int i=0;i<nf;i++)in1>>(*fsz)[i];
		in1>>dump;
		for(int i=0;i<nf2;i++)in1>>(*Sr)[i];
	}
	in1.close();
	
	
	    // for(int i=0;i<nf;i++){
            
            // outdebug << " " << (*fsz)[i];
            
            // }
       // outdebug << endl;     
         // for(int i=0;i<nf2;i++){
            
            // outdebug << " " << (*Sr)[i];
            
            // }
         // outdebug << endl;
    // outdebug.close();
	return inp;
}

void bottom::writeBottom(){
	/*write bottom for restart after Sobek computation*/
	int sepflag=(*fsz)[nf-2];
	int nfsz=(*fsz)[nf-1];
	ofstream bb("out_bottom.inp");
	bb.precision(16);
	bb<<tijd<<" "<<sepflag<<" "<<nfsz<<" "<<0<<" "<<H<<" "<<L<<" "; for(int i=0;i<(*b).size();i++)bb<<(*b)[i]<<" "; bb<<endl;
	if (sepflag==1){
		bb<<tijd<<" "; for(int i=0;i<(*fsz).size();i++)bb<<(*fsz)[i]<<" "; bb<<endl;
		bb<<tijd<<" "; for(int i=0;i<(*Sr).size();i++)bb<<(*Sr)[i]<<" "; bb<<endl;
	}
	bb.close();
}

/*
======================================================
======================================================

BLOCK II: CHECK FOR FLOW SEPARATION
loop over Npx bottom points, to check for flow separation
case 1: a FSZ existed at the previous time step
case 2: a new flow separation zone is formed
        if (dhdx < sepcritangle) -> flow sep
(*fsz) is filled for later use

======================================================
======================================================
*/

void bottom::checkFlowsep(){
	vector<int> fsz_prev = (*fsz);        // previous fsz characteristics
	for(int i=0;i<Npx;i++)(*x)[i]=i*dx;
	(*bp)=(*b);                           // vector for parameterized bottom
    vector<int> dta(2,0);
	int sepflag; //=(*fsz)[nf-2];
	int sepflag1=fsz_prev[nf-2];
	int xsi=-1;	int xri=0; int xci=0; int xti=0; int xdi;
	int nfsz=0;
	int nfsz1=fsz_prev[nf-1];
	int col=4; //xdi
	int cnt=0;
	int solve_method=0;
	int nmerge=0;  // keeps track of # of fsz's that have merged (and removed)
	int skipped=0; // keeps track # of skipped fsz that are too small
	int newwave=0;  // new fsz's due to wavelet generation
	int wavelet=0;
	double sepcritangle1 = tan(cfg.sepcritangle);
	int iinit = 0;

//	ofstream outdebug;
//    outdebug.open ("out_debug1.txt", ofstream::out | ofstream::app);
//    outdebug.precision(16);
//    for(int i=0;i<nf;i++){
//            
//            outdebug << " " << (*fsz)[i];
//            
//            }
//       outdebug << endl;     
//         for(int i=0;i<nf2;i++){
//            
//            outdebug << " " << (*Sr)[i];
//            
//            }
//         outdebug << endl;
//    outdebug.close();


	cout << endl << endl << "Block II: check for flow sep" << endl << endl; //OLAV
	
	/* determination of bed gradients dhdx */
	vec dhdx(Npx,0.0);
    for(int i=0;i<Npx;i++) dhdx[i]=((*b)[i]-(*b)[o2(i-1)])/(dx);
        
	/* als het begin van het domein in de vorige tijdstap een fsz was
	 * dan moet dit gedeelte overgeslagen worden; anders wordt er een
	 * fsz op de lij-zijde van het duin gevonden */
	if (sepflag1==1){ //Olav: opens if (sepflag1==1) check if start was fsz
			
		int xsi_prev=fsz_prev[(nfsz1-1)*7+0];
		int xri_prev=fsz_prev[(nfsz1-1)*7+1];
		int xdi_prev=fsz_prev[(nfsz1-1)*7+4];
		if (xdi_prev>=xsi_prev && xri_prev<xsi_prev){ //Olav: opens if (xdi_prev>=xsi_prev && xri_prev<xsi_prev) "checkFlowsep: case 0"
			cerr<<"checkFlowsep: case 0"<<endl;
	 		iinit=xri_prev;
			/* iinit kan ook binnen een statische fsz vallen
			 * die aan het begin van het domein zit: dan moet iinit
			 * 1 fsz opschuiven naar downstream */ //OLAV 02-21-2010 added */
			int xsi_prev_next=fsz_prev[0];
			int xri_prev_next=fsz_prev[1];
			int xdi_prev_next=fsz_prev[4];
			if (iinit>=xsi_prev_next && iinit<=xri_prev_next){ //OLAV: opens
				iinit=xdi_prev_next;
				cerr<<"iinit reset since iinit is in between a static fsz"<<endl;
                } // closes if (iinit>=xsi_prev_next && iinit<=xri_prev_next)
				cerr<<"iinit set to: "<<iinit<<endl;
	 		} // closes if (xdi_prev>=xsi_prev && xri_prev<xsi_prev)
	} // closes if (sepflag1==1)

	/* here starts the loop over the complete bed, two cases, see block description */
	for(int i=iinit;i<Npx;i++){

		/* case 1: fsz op dezelfde locatie als in de vorige tijdstap */
		
		//OLAV 2014 01 09 
		//Original: if (i==fsz_prev[col] && fsz_prev[nf-2]==1){ 
		//Above Correct? Seems to point to xdi, should be xsi. 
		//if (i==fsz_prev[(nfsz1-1)*7+0] && fsz_prev[nf-2]==1){ //now points to xsi
		//Above correct? 
		if (i==fsz_prev[col] && fsz_prev[nf-2]==1){ 
			sepflag=1;
			nfsz++;
			
			//OLAV 2014 01 09 
			//Original: xsi = fsz_prev[col] 
			//Above Correct? Seems to point to xdi, should be xsi. 
			//xsi = fsz_prev[(nfsz1-1)*7+0]; 
			//Above correct? 
			
			xsi = fsz_prev[col];
			col+=7;
			//col=col+skipped*7;
			
			cerr << nfsz << " " << xsi << " " << col << " " << nfsz1 << endl;
			
			cerr<<"checkFlowsep case 1 met i_in: "<<i<<endl;
			cerr<<"checkFlowsep case 1 met xsi_in: "<<xsi<<" - nfsz_in: "<<nfsz<< " - wavelet_in: " << wavelet<<endl;
			
			dta=setFSZ(xsi,nfsz,wavelet);  // het zetten van sepzone karakteristieken
			i=dta[0];

			cerr<<"checkFlowsep case 1 met xsi: "<<xsi<<" - xri(temp): "<<i<<endl;

			if (i!=-1) { // else fsz too small

				(*fsz)[(nfsz-1)*7+6]=1;  // identifier case 1
				xri=(*fsz)[(nfsz-1)*7+1];
				//cerr<<"   xri: "<<xri<<" (i: "<<i<<")"<<" col: "<<(col-4)/7<<endl;

				// following is used to force to use SOLVE as flow routine
				// if position of xri changes abruptly
				int xsip=fsz_prev[(nfsz+nmerge-1+skipped-newwave)*7+0];
				int xrip=fsz_prev[(nfsz+nmerge-1+skipped-newwave)*7+1];
				int xsin=xsi;
				int xrin=o3(xri);
				if (xsin>xrin && xsip<xrip) {xrin+=Npx;}
				int xri_dif=xrin-xrip;
				if (xri_dif>3) {
					solve_method++;
					cerr<<"   WARNING: [j="<<nfsz<<"] solve_method set to "<<solve_method<<" (xri_dif="<<xri_dif<<")"<<endl;
					outlog<<"T="<<tijd<<" - WARNING: [j="<<nfsz<<"] solve_method set to "<<solve_method<<" (xri_dif="<<xri_dif<<")"<<endl;
					//cerr<<"   Parameters for xri_dif: "<<xsip<<" "<<xrip<<" "<<xsin<<" "<<xrin<<endl;
				}

				/* check for fszs merge */
				int merge=0; // initieer merge op 0: default geen merging
				int mfsz;
				if (nfsz1>1) {

					int kk; int col2=col-4;

					//cerr<<col2/7<<" "<<nfsz1<<endl;

					if (col2/7==nfsz1) {
							//cerr<<"ik kom hier A"<<endl;
							col2=0;} // dan domein gehad

					/* het checken werkt het beste door in upstream richting te checken
				   er kunnen namelijk meerdere fsz's overlappen */
					for (int k=col2/7-2;k>=(col2/7-1)-(nfsz1-1);k--) {
						kk=k;
						if (kk<0) {
								//cerr<<"ik kom hier B"<<endl;
								kk+=nfsz1;}
						int xsi1=fsz_prev[kk*7+0];
						int xri1=fsz_prev[kk*7+1];
						if (xsi1>xri1) xri1+=Npx;
						//cerr<<"kk :"<<kk<<" - xri: "<<xri<<" - xsi1: "<<xsi1<<" - xri1: "<<xri1<<endl;
						/* valt xri tussen downstream xsi en xri?  (gelijk->dan merge) */
						//if (xri>xsi1 && xri<xri1) {
						//op 16/04/2007 vervangen door:
						// maar op 26/04/2007 weer teruggezet!
						if (xri>xsi1 && xri<xri1) {
							//cerr<<"ik kom hier C"<<endl;
							xri=xri1; merge=1;
							/* checken of ie ook niet tussen een volgende fsz zit */
							while (xri>=fsz_prev[(kk+1)*7+0] && xri<=fsz_prev[(kk+1)*7+1]) {
								xri=fsz_prev[(kk+1)*7+1];
								kk++;
								break;
							}
							//cerr<<"kk: "<<kk<<endl;
							//cerr<<"nfsz: "<<nfsz<<endl;
							//cerr<<"nfsz1: "<<nfsz1<<endl;
							mfsz=nfsz;
							//kk-=skipped;
							if (kk<(nfsz-newwave)) {
									//cerr<<"ik kom hier D"<<endl;
									mfsz-=nfsz1; i=Npx;}
							else {
									//cerr<<"ik kom hier E"<<endl;
									i=xri;}
							//i=xri;
							//write_flowsep();
							cerr<<"   WARNING: [j="<<nfsz<<"] overlapping fsz's: "<<nfsz<<" -> ("<<mfsz+1<<"-"<<kk+1<<")"<<"; i set to: "<<i<<"; col: "<<(col-4)/7<<endl;
							outlog<<"T="<<tijd<<" - WARNING: [j="<<nfsz<<"] overlapping fsz's: "<<nfsz<<" -> ("<<mfsz+1<<"-"<<kk+1<<")"<<"; i set to: "<<i<<"; col: "<<(col-4)/7<<endl;
				 			break;
				 		}
						else if (xri>=xsi1 && xri<xri1){ // xri==xsi
							//cerr<<"ik kom hier F"<<endl;
							i--;}
					} // for

					/* bij het onderstaande overlappen 2 FSZ's
					 * de static FSZ moet dan ook nog blijven "bestaan", zonder te migreren */
					if (merge==1) {
						//int mm;
						for (int mm=mfsz+nmerge-newwave+skipped;mm<=kk;mm++) {
							cerr<<"mm: "<<mm<<endl;
							cerr<<"mfsz: "<<mfsz<<endl;
 							int cnum=4; // case number
							if (fsz_prev[mm*7+6]==5) cnum=5;
							int mmfsz=mfsz-skipped;
							if (nfsz>mfsz) mmfsz-=newwave;
							(*fsz)[mmfsz*7+6]=cnum;
							for (int j=0;j<6;j++){
								(*fsz)[mmfsz*7+j]=fsz_prev[mm*7+j];
								(*Sr)[mmfsz]=(*Sr)[mm];
							}
							//cerr<<"ik kom hier.."<<endl;
							mfsz++; col+=7;
							if (i!=Npx) {nfsz++;} // anders hebben we domein al gehad
							//cerr<<mfsz<<": xri: "<<(*fsz)[(mfsz-1)*7+1]<<" (i: "<<i<<")"<<" col: "<<(col-4)/7<<endl;
						}
					} // if (merge==1)

					/* bij het onderstaande mergen 2 FSZ's
					 * de mergende moet dan ook verwijderd worden */
					else if (merge==0 && fsz_prev[(col2/7+0)*7+6]>=4) {
						//cerr<<"ik kom hier - 1"<<endl;
						//cerr<<"col: "<<(col-4)/7<<endl;
						//int xsii=(*fsz)[(col2/7-1)*7+0];
						//int xrii=(*fsz)[(col2/7-1)*7+1];
						int xsii=xsi;
						int xrii=xri;
						if (xrii<xsii) xrii+=Npx;
						int xsi1=fsz_prev[(col2/7+0)*7+0];
						int xri1=fsz_prev[(col2/7+0)*7+1];
						if (xri1<xsii) xri1+=Npx;
						if (xsi1<xsii) xsi1+=Npx;
						//cerr<<col2/7-1<<endl;
						//cerr<<"xsii: "<<xsii<<"; xrii: "<<xrii<<"xsi1: "<<xsi1<<"; xri1: "<<xri1<<endl;
						// for case = 5 moet ie er wel tussenin liggen
						if (fsz_prev[(col2/7+0)*7+6]==5 && (xsii<xsi1 && xrii>=xri1)) {
							//cerr<<"ik kom hier - 2"<<endl;
							col+=7; col2+=7;
							nmerge++;
							while (xrii<xri1) {
								//cerr<<"ik kom hier - 3"<<endl;
								nmerge++;
								col+=7; col2+=7;
								xri1=fsz_prev[(col2/7+0)*7+1];
								if (xri1<xsii) xri1+=Npx;
								if (col2/7>=nfsz1) col2=0-nfsz*7; // dan domein gehad
							}
							cerr<<"   WARNING: static fsz's ("<<nfsz+1<<"-"<<col2/7<<") merge with fsz ("<<nfsz<<") (no need for removal, goes automatically)"<<endl;
							outlog<<"T="<<tijd<<" - WARNING: static fsz's ("<<nfsz+1<<"-"<<col2/7<<") merge with fsz ("<<nfsz<<") (no need for removal, goes automatically)"<<endl;
							solve_method++;
							cerr<<"   WARNING: solve_method set to "<<solve_method<<endl;
							outlog<<"T="<<tijd<<" - WARNING: solve_method set to "<<solve_method<<endl;
							/* Als het goed is wordt hier doorgelinkt naar "check for abandonned fsz's */
						} // if (merge==1)
						//else: do nothing

					} //if (merge==1)

				} //if (nfsz1>1)

			} // (i!=-1)
			// if fsz too small to account for
			else {
				i=findTrough(xsi,filter(3,(*bp)));
				if (i<xsi) i=Npx;
				nfsz--;
				skipped++;
				cerr<<"   WARNING: fsz too small: (xsi="<<xsi<<", i set to "<<i<<")"<<endl;
				outlog<<"T="<<tijd<<" - WARNING: fsz too small: i set to "<<i<<endl;
			} // else
		} //end case 1

		else if(dhdx[i]<sepcritangle1){

		/* case 2: fsz omdat dhdx < sepcritangle */

			sepflag=1;
			nfsz++;
			xsi = i-1;

			cerr<<"checkFlowsep case 2 met xsi: "<<xsi<<endl;

			if (dhdx[o2(xsi+1)]<-0.50) { //dan wavelet
				// dit stuk: initieel steil (dus opgelegde bodem)
				wavelet=1;
				cerr<<"ik kom hier"<<endl;
				cerr<<"xsi: "<<xsi<<endl;
				cerr<<"nfsz: "<<nfsz<<endl;
				cerr<<dhdx[o2(xsi)]<<" "<<dhdx[o2(xsi+1)]<<" "<<dhdx[o2(xsi+2)]<<endl;
				dta=setFSZ(o2(xsi),nfsz,wavelet);
				i=dta[0];
				nfsz=dta[1];
				xri=i;
				newwave++;
				if (newwave>2) {
						break;
						//return 0;
				}
				solve_method++;
				cerr<<"   WARNING: [j="<<nfsz<<"] fsz behind a newly formed wavelet (newwave="<<newwave<<" & solve_method="<<solve_method<<")"<<endl;
				outlog<<"T="<<tijd<<" - WARNING: [j="<<nfsz<<"] fsz behind a newly formed wavelet (newwave="<<newwave<<" & solve_method="<<solve_method<<")"<<endl;
				//dit stuk: initieel steil
				cerr<<"nfsz: "<<nfsz<<endl;
				cerr<<"i: "<<i<<endl;
				//write_flowsep();
				wavelet=0;
			}
			else {
			// dit stuk: development
			xci=findCrest(o2(xsi),filter(1,(*b)));
			xti=findTrough(xsi,filter(3,(*b)));
			// check if fsz needs to be skipped
			// reattachment point of small dune on lee, so double flow sep point found
			cerr<<"ik kom hierrr"<<endl;
			int skip=0;
			int xsii=xsi+1;
			int xcii=o2(xci-1);
			if (xsi==-1) xsii=Npx-1;
			for(int j=0;j<nfsz1;j++){
				int xdip=fsz_prev[j*7+4];
				int xsip=fsz_prev[j*7+0];
				int xrip=fsz_prev[j*7+1];
				cerr<<"j: "<<j+1<<" - xdip: "<<xdip<<" - xsii: "<<xsii<<" - xci: "<<xci<<endl;
				if (xdip==xsii || xdip==xcii) {
					cerr<<endl<<endl<<endl<<endl<<endl<<"         ------ PUNT WAAR HET MIS GING!! -------"<<endl<<endl<<endl<<endl<<endl;
					cerr<<"   WARNING: [j="<<nfsz<<"] sep point is equal to previous, and has to be skipped"<<endl;
					cerr<<"j: "<<j+1<<" - xdip: "<<xdip<<" - xsii: "<<xsii<<" - xcii: "<<xcii<<endl;
					outlog<<"T="<<tijd<<" - WARNING: [j="<<nfsz<<"] sep point is equal to previous, and has to be skipped"<<endl;
					skip=1;
				}
			}
			xsi=xci;
			if ((nfsz>1 && xsi==(*fsz)[(nfsz-2)*7+2]) || skip==1){
				//than sep point is equal to previous, and has to be skipped
				cerr<<"   WARNING: [j="<<nfsz<<"] sep point is equal to previous, and has to be skipped"<<endl;
				outlog<<"T="<<tijd<<" - WARNING: [j="<<nfsz<<"] sep point is equal to previous, and has to be skipped"<<endl;
				nfsz--;
			}
			else {
				xri=xti;
				cerr<<"   WARNING: [j="<<nfsz<<"] xri set to xti since first time flowsep case 2"<<endl;
				outlog<<"T="<<tijd<<" - WARNING: [j="<<nfsz<<"] xri set to xti since first time flowsep case 2"<<endl;
				(*fsz)[(nfsz-1)*7+2]=xci;
				(*fsz)[(nfsz-1)*7+3]=xti;
				(*fsz)[(nfsz-1)*7+0]=xsi;
				(*fsz)[(nfsz-1)*7+1]=xri;
				(*fsz)[(nfsz-1)*7+6]=2;

				if (xti>i) i=xti;
				else if (xti<i) i=Npx;

				cerr<<"   xri: "<<xri<<" (i: "<<i<<")"<<endl;
			}
			}
		} // end case 2

	}  // end for loop over bottom points

    /* xdi gaat niet goed: NB wat was dit ook al weer?
	 * volgens mij komen we hier niet meer.
	 * checken door weg te schrijven naar het log-bestand */
	if (fsz_prev[col]==0 && fsz_prev[col-4]>0){ //Olav: opens if (fsz_prev[col]==0 && fsz_prev[col-4]>0)
			cerr<<"   WARNING: checkFlowsep: case 1 extra"<<endl;
			outlog<<"T="<<tijd<<" - WARNING: checkFlowsep: case 1 extra"<<endl;
    	    sepflag=1;
			nfsz++;
			xsi = fsz_prev[col];
			col+=7;
			dta=setFSZ(xsi,nfsz,wavelet);
			nfsz=dta[1];
			(*fsz)[(nfsz-1)*7+6]=3;
	} //Olav: closes if (fsz_prev[col]==0 && fsz_prev[col-4]>0)

	/* check for merged fsz's at the beginning of *fsz that have to be removed from array */
	int c1=(*fsz)[6];
	int c2=fsz_prev[6];
	if (c1==4 || c2==5) {
		int xsi1=(*fsz)[0];
		int xri1=(*fsz)[1]; if(xsi1>xri1){xsi1-=Npx;}
		int xsi2=(*fsz)[(nfsz-1)*7+0];
		int xri2=(*fsz)[(nfsz-1)*7+1]; if(xsi2>xri2){xsi2-=Npx;}
		//cerr<<"xsi1 "<<xsi1<<"; xri1: "<<xri1<<"xsi2: "<<xsi2<<"; xri2: "<<xri2<<endl;
		if (xri2>=xri1 && xsi2<=xsi1) {
			cerr<<"   WARNING: static flow separation zone is merged and removed from array"<<endl;
			outlog<<"T="<<tijd<<" - WARNING: static flow separation zone is merged and removed from array"<<endl;
			for (int i=0;i<7;i++){
				(*fsz)[i]=(*fsz)[(nfsz-1)*7+i];
				(*Sr)[0]=(*Sr)[(nfsz-1)];
			}
		nfsz--;
		}
	} //Olav: closes if (c1==4 || c2==5) 

	if(nfsz==0)sepflag=0;

	/* store some parameters in *fsz array */
	(*fsz)[nf-3]=solve_method;
	(*fsz)[nf-2]=sepflag;
	(*fsz)[nf-1]=nfsz;

	/* set old fsz�s to zero in array (*fsz) */
	for(int i=nfsz*7;i<nf-3;i++){
		(*fsz)[i]=0;
	}

	/* initialize rounding errors to zero */
	for(int i=nfsz;i<nf2;i++){
		(*Sr)[i]=0;
	}

	/* write separation zone characteristics to screen */
	if (sepflag==1) write_flowsep();

	/* max of bed and param bed: is bedflow in main */
	for(int i=0;i<Npx;i++) {
		(*bp)[i]=max((*bp)[i],(*b)[i]);}

	/* smooth bed to avoid strong bed gradients at separation and reattachment */
	smooth_param(5,1); // smooth at xsi
	smooth_param(5,2); // smooth at xri

	if (sepflag==0) cerr<<"Minimum dhdx: "<<atan(minval(dhdx,Npx))*grad_2_deg<<" degrees"<<endl;
} 

/*
======================================================
======================================================

BLOCK III: DETERMINATION OF THE SHAPE OF THE FSZ
Contains functions:
-setFSZ: determine bed characteristics of FSZ
-paramSepline: determine separating streamline

======================================================
======================================================
*/

vector<int> bottom::setFSZ(int xsi, int nfsz, int wavelet){

  vector<int> dta(2,0);
	int xri=0; int xci=0;	int xti=0;
	xti=findTrough(xsi,filter(3,(*bp)));
	xci=findCrest(xsi,filter(1,(*bp)));

	/*
	int tel=0;
	int xsii=xsi; //store original xsi
	// locatie xsi corrigeren als vlak stukje door sep_migr_lee.
	// NB uitgezet, werkt niet goed!!
	while ((*b)[xsi]==(*b)[o2(xsi-1)]){// && (*b)[o2(xsi-1)]!=(*b)[o2(o2(xsi-1)-1)]){
		xsi=o2(xsi-1);
		cerr<<"xsi-- voor param. sepline"<<endl;
		tel++;
	}

	if (tel>1) {
		xsi=xsii;
		cerr<<"param sepline reset"<<endl;}
	*/

	xri=paramSepline(xsi,xti,xci,nfsz);

	/*
	if (xri-xsi<=1) {
		cerr<<"   WARNING: [j="<<nfsz<<"] FSZ too small: neglected!"<<endl;
		outlog<<"T="<<tijd<<" - WARNING: [j="<<nfsz<<"] FSZ too small: neglected!"<<endl;
		nfsz--;
		dta[0]=-1; dta[1]=nfsz;
		//cerr<<"   WARNING: [j="<<nfsz<<"] xri++, since xri-xsi=1!"<<endl;
		//outlog<<"T="<<tijd<<" - WARNING: [j="<<nfsz<<"] xri++, since xri-xsi=1!"<<endl;
		//xri++;
	}
	*/

	int m=xsi;

	int xrin=xri;
	if (xri!=-1 && xrin<xsi) xrin+=Npx; // periodic bcs
	if (xrin==-1 || xrin-xsi<=1){
		// xri -1 from paramSepline | too small FSZ: then neglect
		cerr<<"   WARNING: [j="<<nfsz<<"] xri set to -1 in paramSepline(...), since FSZ too small"<<endl;
		cerr<<"xsi: "<<xsi<<"; xri: "<<xri<<endl;
		outlog<<"T="<<tijd<<" - WARNING: [j="<<nfsz<<"] xri set to -1 in paramSepline(...), since FSZ too small"<<endl;
		nfsz--;
		dta[0]=-1; dta[1]=nfsz;
	}
	else {
		(*fsz)[(nfsz-1)*7+0]=xsi;
		(*fsz)[(nfsz-1)*7+1]=xri;
		(*fsz)[(nfsz-1)*7+2]=xci;
		(*fsz)[(nfsz-1)*7+3]=xti;
		/*
		cerr<<"xsi: "<<xsi<<"; xri: "<<xri<<endl;
		cerr<<"alpha at xsi-1: "<<atan(dhdx[o2(xsi-1)])*grad_2_deg<<endl;
		cerr<<"alpha at xsi  : "<<atan(dhdx[xsi])*grad_2_deg<<endl;
		cerr<<"alpha at xsi+1: "<<atan(dhdx[o2(xsi+1)])*grad_2_deg<<endl;
		*/
		//if ( ( xri<xsi | xti<xsi ) && wavelet==0 ){  // dan domein gehad!
		if ( xri<xsi ) {
			cerr<<"ik kom hierr"<<endl;
			m=Npx;}
		else
			{m=xri;}

		dta[0]=m;
	} // if (xri==-1)

	//cerr<<"m: "<<m<<endl;
	//cerr<<"xri: "<<xri<<endl;
	dta[1]=nfsz;
	return dta;
}

/* reimplemenation of this function after
 * 10/04/2007
 */
int bottom::paramSepline(int xsi, int xti, int xci, int nfsz){
	/* bepaling separating streamline
	 * shape of the separation streamline as determined from the data-set,
	 * described in Paarlberg et al. 2007 submitted WRR
	 * this describes a function for the dimensionless length of the FSZ in terms
	 * of bed properties at the separation point. Also relations are given for the
	 * four coefficients of a third order polynomial.
	 * the cross-point between the separation streamline and the bed is determined
	 * using a numerical algorithm
	 */

	// determination of tan alpha_s & xsi
	vec dhdx(Npx,0.0);
  for(int i=0;i<Npx;i++){
   	dhdx[i]  = ( (*b)[i] - (*b)[o2(i-1)] )/(dx);
  }

	int pt=xsi;
	while ((*b)[pt]==(*b)[o2(pt-1)]){// && (*b)[o2(xsi-1)]!=(*b)[o2(o2(xsi-1)-1)]){
		pt=o2(pt-1);
		cerr<<"pt-- voor param. sepline dhdx bepaling"<<endl;
	}

	// 11-08-2006: TODO: dit kan ook meer dan 1 punt moeten zijn!! (dus een while)
	// if (dhdx[pt]==0) pt--; pt=o2(pt);
	// double alpha_b=(dhdx[o2(o2(pt-1)-1)]+dhdx[o2(pt-1)]+dhdx[pt])/3.;
	double alpha_s=dhdx[pt]/1.;

	//override!!!
	//alpha_s=0;

	// x and z coordinate of the flow separation and reattachment point
	double xs=(*x)[xsi];
	double Hs=(*b)[xsi];
	double xt=(*x)[xti];
	double Ht=(*b)[xti];
	
	// height and total lengte of the flow separation zone
	// total length: point where separation streamline crosses the trough elevation
	double Hfsz = Hs-Ht;
	double tLs = (6.48*alpha_s + 5.17);
		
	// coefficients according to Paarlberg et al (2007, WRR)
	double s0 = 1.;
	double s1 = alpha_s;
	double s2, s3;
	s3=2./pow(tLs,3.)+(-0.53+alpha_s)/(tLs*tLs);
	if (alpha_s<0) s3=0.;
	s2=-s3*tLs-alpha_s/tLs-1./(tLs*tLs);

	// numerical implementation to find crossing point
	double tol  = 1e-10;
	int max_it  = 150;
	int dir = 1;  // 1 = left2right; 2 = right2left
	double xl = xs+dx;
	double xr = xs+Hfsz*tLs;
	double dif = 0.0;
	double x_p = 0.0;
	double x_p_bed;
	double fac=1;
	int pr_dir=dir;

	vec nb(4,0.0);
	vec rp(2,0.0);
	
	for(int k=1;k<=max_it;k++){
		// find cross-point within maximum number of iterations
		dif = xr - xl;
		// 1 = left2right; 2 = right2left
		if (dir==2){
			x_p = xr - dx/fac;}
		else if (dir==1){
			x_p = xl + dx/fac;}
		// find neighboring points of the bed
		nb=paramFindNeighbors(x_p,xsi);

		x_p_bed = x_p;
		if (x_p>=Npx*dx || x_p==L) x_p_bed = x_p-Npx*dx;
		if (x_p<0) x_p_bed = x_p+Npx*dx;
		double y1 = ((nb[1]-x_p_bed)*nb[2]+(x_p_bed-nb[0])*nb[3])/(nb[1]-nb[0]);
		double xx=(x_p-xsi*dx)/Hfsz;
		double y2 = s3*pow(xx,3) + s2*pow(xx,2) + s1*xx + s0; // height of sepline
		y2=y2*Hfsz+Ht;

		if (nb[2]==0 && nb[3]==0) {
			cerr<<"   WARNING: [j="<<nfsz<<"] neighbor points not found correctly!:"<<endl;
			outlog<<"T="<<tijd<<" - WARNING: [j="<<nfsz<<"] neighbor points not found correctly!:"<<endl;
			cerr<<"x_p="<<x_p<<"; xsi="<<xsi<<endl;
			cerr<<"Npx*dx="<<Npx*dx<<"; xr_in="<<xs+Hfsz*tLs<<endl;
			cerr<<nb[0]<<" "<<nb[1]<<" "<<nb[2]<<" "<<nb[3]<<endl;
			cerr<<k<<" "<<dir<<" "<<xs+dx<<" "<<x_p<<" "<<y1<<" "<<y2<<" "<<xx<<" "<<abs(y2-y1)<<endl;
		}

		if (k==max_it) {
			cerr<<"maximum number of iterations reached"<<endl;
			outlog<<"T="<<tijd<<" - WARNING: [j="<<nfsz<<"] maximum number of iterations reached to determine reattachment point!:"<<endl;
			cerr<<x_p<<endl;
		}

		if (abs(y2-y1) <= tol || abs(y2-y1) == 0){
		// cross-point fond or maximum # iteration reached
			int nb1=o3(int(ceil(x_p/dx)));
			rp[0]=x_p;
			rp[1]=y2;
			break;
			}
		else{
			// next iteration
			if (y2 > y1){
				dir = 1;  // left2right
				xl  = x_p;}
			else{
				dir     = 2;  // right2left
				xr = x_p;}
			if (dir!=pr_dir) {fac*=2;}
			pr_dir=dir;
		}
	} // for
	
	double xret = rp[0];
	double Hret = rp[1];
	double Ls = xret - xs;
	int sxi=int(ceil(Ls/dx));
	int xri=xsi+sxi;
	
	/* actual determination of separating streamline */
	vec sx(sxi+1,0.0);
	vec s (sxi+1,0.0);

	for(int i=0;i<=sxi;i++) sx[i]=(i*dx)/Hfsz;
	
	for(int i=0;i<=sxi;i++){
		
		s[i] = s3*pow(sx[i],3) + s2*pow(sx[i],2) + s1*sx[i] + s0;
		s[i] = s[i]*Hfsz+Ht;
	}
	for(int i=xsi;i<=xri;i++) {

		int m = i;
		if (i>=Npx) m = i-Npx;
		(*bp)[m]=s[i-xsi];
	}

	if(xri>Npx) xri-=Npx;
	
	return xri;
}

/*
======================================================
======================================================

BLOCK IV: DETERMINATION OF THE SEDIMENT TRANSPORT
- routine without critical bed shear stress
- routine including critical bed shear stress
both cases take flow separation into account, by determining
the volumetric stress at xsi (flux@xsi)

======================================================
======================================================
*/

// OLAV: function detQ is not used 2011 02 25
//void bottom::detQ(vec ub, vec &dhdx){
//	/* determine sediment fluxes without critical bed shear stress */
//	int sepflag=(*fsz)[nf-2];
//	vec tau(Npx,0.0);
//
//	for(int i=0;i<Npx;i++){
//		tau[i]=S*(ub[i]);
//		dhdx[i]=((*b)[i]-(*b)[o2(i-1)])/dx;
//	}
//
//	for(int i=0;i<Npx;i++){
//		(*flux)[i]=0.0;
//		double taui=(1./2.)*(tau[o2(i-1)]+tau[i]);
//		if (tau[i]>0.) (*flux)[i]+=alpha*pow(taui,be)*(1-l1/taui*(dhdx[i])-cfg.l2*(dhdx[i]));
//	}
//
//	/* in case of flow separation, the flux@xsi, based on local tau in that
//	 * point, is stored in (*flux)[xsi+1] */
//	if (sepflag==1){
//		int nfsz=(*fsz)[nf-1];
//		for (int j=0;j<nfsz;j++) {
//			int xsi=(*fsz)[j*7+0];
//			//cerr<<"j: "<<j<<"tau[xsi="<<xsi<<"]: "<<tau[xsi]<<endl;
//			if (tau[xsi]>0.) (*flux)[o2(xsi+1)]+=alpha*pow(tau[xsi],be)*(1-l1/tau[xsi]*(dhdx[xsi])-cfg.l2*(dhdx[xsi]));
//		}
//	}
//}

void bottom::detQcr(vec ub, vec &dhdx){
	/* determine fluxes in case of flowsep some extra code is used*/
	int sepflag=(*fsz)[nf-2];
	vec tau(Npx,0.0);
	double sepflux; // OLAV 2013 04 16
	double meanstle1;  // OLAV 2014 02 25
	double alpha_lag1;  // OLAV 2014 02 25
	
//        // ADDED 2011 2 25 (OLAV)
//        ofstream outdebug;
//	    outdebug.open ("out_debug1.txt", ofstream::out | ofstream::app);
//	    outdebug.precision(16);
//           ostringstream tmpbot10;
//           tmpbot10 << "out_debug" << 1 << ".txt";
//           string ofname10 = tmpbot10.str();
//           ofstream outdebug(ofname10.c_str(),ios_base::out);
//        // END ADDED 2011 2 25 (OLAV)
// outdebug << "xsi= " << xsi << " o2(xsi+1)= " << o2(xsi+1) << " (*flux)[o2(xsi+1)]= " << (*flux)[o2(xsi+1)] << endl; // ADDED 2011 2 25 (OLAV)
// ADDED 2011 2 25 (OLAV)	
//    outdebug.close();
// END ADDED 2011 2 25 (OLAV)
               
	for(int i=0;i<Npx;i++){
		tau[i]=S*(ub[i]);
		dhdx[i]=((*b)[i]-(*b)[o2(i-1)])/dx;
	}
    
	for(int i=0;i<Npx;i++){
		(*flux)[i]=0.0;      
		double taui=(1./2.)*(tau[o2(i-1)]+tau[i]);
		
		// Meyer-Peter M�ller (original)
        if(cfg.transport_eq == 1 || cfg.transport_eq == 3){
           double tauc=cfg.thetacr*cfg.g*cfg.delta*cfg.D50*((1.+cfg.l2*dhdx[i])/(pow(1.+dhdx[i]*dhdx[i],(1./2.)))); //OLAV 2011 02 24
		   //ORIGINAL: double tauc=cfg.thetacr*g*(2.65-1.)*cfg.D50*((1.+cfg.l2*dhdx[i])/(pow(1.+dhdx[i]*dhdx[i],(1./2.))));
		   
           if (tau[i]>tauc && taui>0.) {
			   taui=max(taui,tauc); //OLAV: 2011 2 23 --> why this operation?
			   (*flux)[i]+=cfg.alpha*pow(taui-tauc,cfg.be)*(1./(1.+cfg.l2*(dhdx[i])));
           }
            
           /* in case of flow separation, the flux@xsi, based on local tau in that
	       * point, is stored in (*flux)[xsi+1] 
		   * this is a correction afterwards    */
           if (sepflag==1 ){ //OLAV 2013 04 16 // && AllowAvalanching==1
		      int nfsz=(*fsz)[nf-1];
		      for (int j=0;j<nfsz;j++) {
			           int xsi=(*fsz)[j*7+0];
                       (*flux)[o2(xsi+1)]=0.0;
                       //ORIGINAL: tauc=cfg.thetacr*cfg.g*(2.65-1.)*cfg.D50*((1.+cfg.l2*dhdx[xsi])/(pow(1.+dhdx[xsi]*dhdx[xsi],(1./2.))));
                       double tauc=cfg.thetacr*cfg.g*cfg.delta*cfg.D50 *((1.+cfg.l2*dhdx[xsi])/(pow(1.+dhdx[xsi]*dhdx[xsi],(1./2.))));
			
			           if (tau[xsi]>tauc && tau[xsi]>0.) {
                          (*flux)[o2(xsi+1)]+=cfg.alpha*pow(tau[xsi]-tauc,cfg.be)*(1./(1.+cfg.l2*(dhdx[xsi])));
				          // ORIGINAL: (*flux)[o2(xsi+1)]+=alpha*pow(tau[xsi]-tauc,be)*(1./(1.+cfg.l2*(dhdx[xsi])));
                       }
              }
           } // CLOSES if (sepflag==1)
        } // CLOSES if(transport_eq == 1)
        
        // Nakagawa / Tsujimoto
		       
        else if(cfg.transport_eq == 2 ){
		   //non-corrected tauc
		   double tauc=cfg.thetacr*cfg.g*cfg.delta*cfg.D50; //OLAV 2013 07 08
		   //double tauc=cfg.thetacr*cfg.g*cfg.delta*cfg.D50*((1.+cfg.l2*dhdx[i])/(pow(1.+dhdx[i]*dhdx[i],(1./2.)))); //OLAV 2014 02 07
		   
		   double thetai = taui / (cfg.g*cfg.delta*cfg.D50);

           if (tau[i]>tauc && taui>0.) {
			   thetai=max(thetai,cfg.thetacr);
			   //thetai=max(thetai,cfg.thetacr*((1.+cfg.l2*dhdx[i])/(pow(1.+dhdx[i]*dhdx[i],(1./2.))))); //OLAV 2014 02 07

			   (*flux)[i]=cfg.F0_dim*thetai*pow(1.-cfg.thetacr/thetai,3.); //pickup

			   //(*flux)[i]=F0_dim*thetai*pow(1.-((1.+cfg.l2*dhdx[i])/(pow(1.+dhdx[i]*dhdx[i],(1./2.))))*cfg.thetacr/thetai,3.);
           }
		   
        } // CLOSES if(transport_eq == 2)
		/*
        else if(transport_eq == 2 ){
		   //non-corrected tauc
		   double tauc=cfg.thetacr*g*cfg.delta*cfg.D50; //OLAV 2013 07 08
		   //corrected tauc
		   //double tauc=cfg.thetacr*g*cfg.delta*cfg.D50*((1.+cfg.l2*dhdx[i])/(pow(1.+dhdx[i]*dhdx[i],(1./2.)))); //OLAV 2011 02 24
		   
		   double alpha_2 = correction*F0*pow(g*cfg.delta/cfg.D50,(1./2.));
           if (tau[i]>tauc && taui>0.) {
			   taui=max(taui,tauc); 
			   double thetai = taui / (g*cfg.delta*cfg.D50);
			   (*flux)[i]=alpha_2*thetai*pow(1.-tauc/taui,3.);
           }         
        } // CLOSES if(transport_eq == 2) */
		
    } // CLOSES for(int i=0;i<Npx;i++)
	  
    //Fluxes are determined, now optionally a lag is applied. 
    if(cfg.transport_eq == 3 ){
	
		if(cfg.alpha_varies>0){
			alpha_lag1=detAlphaLag(ub,cfg.alpha_varies,0);
			meanstle1 = alpha_lag1*cfg.D50;
		} 
		else {
			meanstle1 = cfg.meanstle;
		}

                   
		double dxtenth=dx/10;
		double flux_eq_k;
		double flux_k;
		double flux_k_prev;
				
		vec flux_eq(Npx,0.0);
		flux_eq=(*flux);
				
		double flux_temp; //changed 2012 09 07 OLAV, was: 
				
		if(cfg.moeilijkdoen==0){
			//Guess flux(0)
			(*flux)[0]=flux_eq[0]/2;
				
			for(int i=1;i<Npx;i++){ //From 1, because 0 was already guessed
				//Backwards Euler
				(*flux)[i]=((*flux)[i-1]+dx/meanstle1*flux_eq[i])/(1+dx/meanstle1);
			} //closes for(int i=1;i<Npx;i++)
						
			flux_temp = ((*flux)[Npx-1]+dx/meanstle1*flux_eq[0])/(1+dx/meanstle1);//changed 2012 09 07 OLAV, was: 
				
			for(int j=0;j<30;j++){
				//Re-guess flux(0)
				(*flux)[0]=((*flux)[0]+flux_temp)/2;
				//changed 2012 09 07 OLAV, was: (*flux)[0]=((*flux)[0]+(*flux)[Npx-1])/2;
							
				//(*flux)[0]=(*flux)[Npx-1];
				for(int i=1;i<Npx;i++){
					//Backwards Euler
					(*flux)[i]=((*flux)[i-1]+dx/meanstle1*flux_eq[i])/(1+dx/meanstle1);
				} //closes for(int i=1;i<Npx;i++)
				//outdebug << "j=" << j << " (*flux_eq)[0]=" << flux_eq[0] << " (*flux)[0]=" << (*flux)[0] << " (*flux)[L]="<< (*flux)[Npx-1];
				//outdebug << endl;
							
				flux_temp = ((*flux)[Npx-1]+dx/meanstle1*flux_eq[0])/(1+dx/meanstle1);;//changed 2012 09 07 OLAV, was: 
							
				//changed 2012 09 07 OLAV, was: if( abs( ((*flux)[0]-(*flux)[Npx-1])/(*flux)[0]) < 0.001 ) {j=30;};
				if( abs( ((*flux)[0]-flux_temp)/(*flux)[0]) < 0.001 ) {j=30;};
					} //closes for(int j=0;j<30;i++)
						   
					//changed 2012 09 07 OLAV, was: (*flux)[0]=(*flux)[Npx-1];
					(*flux)[0]=flux_temp;
							
		}
		else { //i.e. moeilijkdoen = 1
			//Guess flux(0)
			(*flux)[0]=flux_eq[0]/2;
			
			for(int j=0;j<30;j++){                            
				flux_k_prev = (*flux)[0];
				for(int i=0;i<Npx-1;i++){
							
					for(int k_dx=1;k_dx<11;k_dx++) {
						flux_eq_k = (flux_eq[o2(i+1)]-flux_eq[i])*k_dx/10+flux_eq[i];
						flux_k = (flux_k_prev+dxtenth/meanstle1*flux_eq_k)/(1+dxtenth/meanstle1);
						flux_k_prev = flux_k;
					}  //closes for(int k=1;k<11;j++){
							//outdebug.close();
					(*flux)[o2(i+1)]=flux_k;         
				} //closes for(int i=0;i<Npx;i++){
							  
				//Re-guess flux(0)
				(*flux)[0]=((*flux)[0]+(*flux)[Npx-1])/2;
					
				if( abs( ((*flux)[0]-(*flux)[Npx-1])/(*flux)[0]) < 0.001 ) {
					j=30;
					};
				} //closes for(int j=0;j<30;j++)
			(*flux)[0]=(*flux)[Npx-1];
					
		} //close else (if moeilijkdoen = 0)
    } //closes if(transport_eq == 3 )
} // CLOSES function		

/*
======================================================
======================================================

BLOCK V: BOTTOM UPDATE ROUTINES.
These are called from main.
Va routine in case without flow separation: update
Vb routine in case with flow separation   : update_flowsep

======================================================
======================================================
*/

vec bottom::update(vec ub, vec &bss1, vec &fluxtot, vec &dhdx){
	/* bottom-update without flow separation */

	for(int i=0;i<Npx;i++) {bss1[i]=S*ub[i];
							fluxtot[i]=0.; /* fluxtot is required for storage */
	}
	// JW vec newb(Npx,0.0);
	vec oldb(*b);
	vec distribute;
	double bint = 0.0;
	double ep=(1/(1-cfg.epsilonp));
	
	// Needed for N&T
	double alpha_lag1=0.;
	int Npsl = 0;
	int dxsubdiv;
	double dxGhost=0;
	int GhostProtocol=0;
	int ghostNpx=0;
	int ghostNpsl=0;
	vec depositemp(Npx,0.0);
	vec fluxtemp(Npx,0.0);
	
//	if(cfg.alpha_varies==0 && cfg.transport_eq == 2){
//		alpha_lag1=detAlphaLag(ub,cfg.alpha_varies,1);
//
//		distribute=detDistributeFunc(alpha_lag1,dx);
//	}
	
	double depositot =0.;
	double tobedepositot =0.;
	double tobedeposi=0.;
	int depopoint =0;
	vec deposi(Npx,0.0);
					
	for(int t=0;t<int(cfg.tt);t++){
		//Joris: deze loop zit er in omdat expliciet niet al te grote tijdstappen aankan, moet nog aan getweekt worden.
		detQcr(ub, dhdx);
	
	if(cfg.transport_eq == 1 || cfg.transport_eq == 3){ //OLAV: 2013 05 22
		for(int i=0;i<Npx;i++){
				(*b)[i]-=ep*dt/cfg.tt/dx*((*flux)[o2(i+1)]-(*flux)[i]);
				fluxtot[i]+=(*flux)[i]/cfg.tt*ep;
        }
	} //closes if(transport_eq == 1 || transport_eq == 3)
 
	//OLAV: 2013 05 22 START
	else if(cfg.transport_eq == 2){
	
		//Determine alpha and prepare distribution function if needed 
		if(cfg.alpha_varies==0 && t==0){
			alpha_lag1=detAlphaLag(ub,cfg.alpha_varies,1);
			distribute=detDistributeFunc(alpha_lag1,dx);
		}

		if(cfg.alpha_varies!=0){
//		if(cfg.alpha_varies==1 || cfg.alpha_varies == 2 || cfg.alpha_varies == 3){
			//Determine alpha
			alpha_lag1=detAlphaLag(ub,cfg.alpha_varies,t);
			//Determine distribution
			distribute = detDistributeFunc(alpha_lag1,dx);

		} //end of calculation for alpha and distribution
		
		Npsl = distribute.size();
		
		if(t==0){
			cerr << "NPSL: " << Npsl << " dx:" << dx << endl;
		} 
		
		GhostProtocol=0;
		if (Npsl < cfg.Npsl_min){
		
			GhostProtocol=1;

			dxsubdiv=int(ceil(dx/((cfg.stle_factor*alpha_lag1*cfg.D50)/cfg.Npsl_min)));
			dxGhost = dx/dxsubdiv;
			distribute = detDistributeFunc(alpha_lag1,dxGhost);
			ghostNpsl = distribute.size();
			
			ghostNpx = Npx*dxsubdiv;
			depositemp.resize(ghostNpx);
			fluxtemp.resize(ghostNpx);
			
			if(t==0){
				cerr << "Npsl_min: " << cfg.Npsl_min << " ghostNpsl:" << ghostNpsl << " dxGhost:" << dxGhost << endl;
				cerr << "Npx: " << Npx << " ghostNpx: "<< ghostNpx << endl;
			}
			
			for(int i=0;i<Npx-1;i++){ //for i = 1:length(testflux)
				fluxtemp[dxsubdiv*(i)]=(*flux)[i]; //original x point 1
				fluxtemp[dxsubdiv*(i)+dxsubdiv]=(*flux)[i+1]; //original x point 2
                
				for(int j=2;j<dxsubdiv+1;j++){
					fluxtemp[dxsubdiv*(i)-(1-j)]=(*flux)[i]+(j-1)*((*flux)[i+1]-(*flux)[i])/dxsubdiv; 
				}
			}

			fluxtemp[dxsubdiv*(Npx-1)]=(*flux)[Npx-1]; //original last x point 
            
			for (int j=2;j<dxsubdiv+1;j++) {
				fluxtemp[dxsubdiv*(Npx-1)-(1-j)]=(*flux)[Npx-1]+(j-1)*((*flux)[0]-(*flux)[Npx-1])/dxsubdiv; 
			}
		}
		
		//Deposition protocol
		if (GhostProtocol==0){
		for(int i=0;i<Npx;i++){
			depositot = (*flux)[i];    //to keep track of the running total of deposited material
			tobedepositot = depositot; //the goal value of the total deposited material
			
			while(depositot>0){
				for(int j=0;j<Npsl;j++){
					tobedeposi=distribute[j]*tobedepositot; //was tobedeposi=distribute[j]*dx*tobedepositot; Olav 2014 03 17
					
					if(i+j>Npx-1){
						depopoint=i+j;
						while (depopoint > Npx-1){
							depopoint=depopoint-Npx;
						}
					}
					else{depopoint=i+j;}
					
					if(tobedeposi < depositot){depositot=depositot-tobedeposi;}
					else{tobedeposi=depositot;depositot=0;}
					
					depositemp[depopoint] += tobedeposi; //correct with i+j and o3/o2? 
				}
			} // closes while
		} //closes loop over x 
		}
		else if (GhostProtocol==1){
				
		for(int i=0;i<ghostNpx;i++){
			depositot = fluxtemp[i];    //to keep track of the running total of deposited material
			tobedepositot = depositot; //the goal value of the total deposited material
			
			while(depositot>0){
				for(int j=0;j<ghostNpsl;j++){
					tobedeposi=distribute[j]*tobedepositot; //was tobedeposi=distribute[j]*dxGhost*tobedepositot; Olav 2014 03 17
					
					if(i+j>ghostNpx-1){
						depopoint=i+j;
						while (depopoint > ghostNpx-1){
							depopoint=depopoint-ghostNpx;
						}
					}
					else{depopoint=i+j;}
					
					if(tobedeposi < depositot){depositot=depositot-tobedeposi;}
					else{tobedeposi=depositot;depositot=0;}
					
					depositemp[depopoint] += tobedeposi; //correct with i+j and o3/o2? 
				}
			} // closes while
		} //closes loop over x 
		}
		
		// Make the real deposi vector
		if (GhostProtocol==0){
			deposi=depositemp;
		}
		else if (GhostProtocol==1){ 
		    for(int i=0;i<Npx;i++){
				//cerr << dxsubdiv*(i) << endl;
				deposi[i]= depositemp[dxsubdiv*(i)];
			}
		}
		
		/*
		ofstream outdebug2;
		outdebug2.open ("out_debug2.txt", ofstream::out | ofstream::app);
		outdebug2.precision(16);
		
		ofstream outdebug3;
		outdebug3.open ("out_debug3.txt", ofstream::out | ofstream::app);
		outdebug3.precision(16); 
		
		for(int i=0;i<Npx;i++){ 
			outdebug2 << (*flux)[i] << " ";
		}
		outdebug2 << endl;
		for(int i=0;i<Npx;i++){ 
			outdebug2 << deposi[i] << " ";
		}
		outdebug2 << endl;
		
		for(int i=0;i<ghostNpx;i++){ 
			outdebug3 << fluxtemp[i] << " ";
		}
		outdebug3 << endl;
		for(int i=0;i<ghostNpx;i++){ 
			outdebug3 << depositemp[i] << " ";
		}
		outdebug3 << endl;
		*/
		
		//Smoothing OLAV 2014 04 14
		int sepsmooth=0;
		for(int i=0;i<Npx;i++){
			if(dhdx[i]<tan(cfg.sepcritangle)){sepsmooth=1;}
		}
		
		int poscr=-1;
		int postr=-1;
		int pos=0;
		int nsmooth=0;
		
		if (sepsmooth == 100){
			double maxval=-1.e99;
			double minval=1.e99;
			
			for(int i=0;i<Npx;i++){
				if ((*b)[i]>maxval){
					poscr=i;
					maxval=(*b)[i];
				}
			}
			
			// for(int i=0;i<Npx;i++){
				// if((*b)[i]<minval){
					// postr=i;
					// minval=(*b)[i];
				// }
			// }
			for(int i=poscr;i<poscr+Npx;i++){
			
				if (i>Npx-1){pos=i-Npx;}

				if((*b)[pos]<minval){
					postr=pos;
					minval=(*b)[pos];
				}
			}
			
			cerr << poscr << " " << postr << endl;
			if (postr<poscr){postr=postr+Npx;}
			
			double deposmooth=0.;
			double fluxsmooth=0.;
			//poscr=poscr+1; //shift one further
			
			for(int i=poscr;i<postr+1;i++){
			
				nsmooth+=1;
				
				if (i>Npx-1){pos=i-Npx;}
				else{pos=i;}
				
				cerr << pos << " " ;
				
				deposmooth+=deposi[pos];
				fluxsmooth+=(*flux)[pos];
			}
			cerr << endl; 
			
			//int nsmooth=postr-poscr+1;
			for(int i=poscr;i<postr+1;i++){
				
				if (i>Npx-1){pos=i-Npx;}
				else{pos=i;}
				
				cerr << pos << " " ;
				
				deposi[pos]=deposmooth/nsmooth;
				(*flux)[pos]=fluxsmooth/nsmooth;
			}
			cerr << endl; 
		}
		
		fluxtemp.resize(Npx);
		depositemp.resize(Npx);
		distribute.resize(1);
		for(int i=0;i<Npx;i++){
			(*b)[i]-=ep*dt/cfg.tt*cfg.D50*((*flux)[i]-deposi[i]);
			
			fluxtot[i]+=(*flux)[i]/cfg.tt*ep; //not the real flux! Is pick-up only!
			
			deposi[i]=0; 
			depositemp[i]=0;
			fluxtemp[i]=0;
		} //closes loop over x	
		
		if (sepsmooth == 100){
			double maxdh = 1./2.*((*b)[o2(poscr-1)]-(*b)[o2(poscr-2)]);
			int nsandvault =nsmooth-1;
			pos = poscr;
			//double maxdh = (*b)[o2(poscr-2)]-(*b)[o2(poscr-3)];
			//int nsandvault =nsmooth;
			//pos = poscr-1;
			
			double sandvault = 0.;
			if (maxdh > 0){

				// while((*b)[o2(pos)]-(*b)[o2(pos-1)]>maxdh && pos<postr){
					// sandvault = (*b)[o2(pos)]-(*b)[o2(pos-1)] - maxdh;
					// (*b)[o2(pos)] = (*b)[o2(pos-1)]+maxdh;
					// for(int i=pos+1;i<postr+1;i++){
					
						// if (i>Npx-1){pos=i-Npx;}
						// else{pos=i;}
					
						// cerr << pos << " " ;
					
						// (*b)[pos]+=sandvault/nsandvault;
					// }
					
					// nsandvault-=1;
					// pos+=1;
				// }
								
				if ((*b)[o2(poscr)]-(*b)[o2(poscr-1)]>maxdh){ //then the angle towards the crest is too high
					sandvault = (*b)[o2(poscr)]-(*b)[o2(poscr-1)] - maxdh;
					(*b)[o2(poscr)] = (*b)[o2(poscr-1)]+maxdh;
					
					for(int i=poscr+1;i<postr+1;i++){
					
						if (i>Npx-1){pos=i-Npx;}
						else{pos=i;}
					
						cerr << pos << " " ;
					
						(*b)[pos]+=sandvault/nsandvault;
					}	
				}
			}
		}
		
	} //closes else if(transport_eq == 2){
	//OLAV: 2013 05 22 END
	
	} // closes for(int t=0;t<int(cfg.tt);t++){
	
	//avalanching protocol
	if(cfg.transport_eq == 2 && cfg.AllowAvalanching == 1){
		avalanche(); 
	}
			
	auto newb(*b);
	(*b)=oldb;

	return newb;
}

vec bottom::update_flowsep(vec ub, vec &bss1, vec &bss2, vec &fluxtot, vec &q_spec){
	/* bottom-update with flow separation */

	for(int i=0;i<Npx;i++) {bss1[i]=S*ub[i];}
	vec newb(Npx,0.0);
	vec oldb=(*b);
	
	//int Npsl = (int)ceil(meanstle*5/dx); //needed for x-dependent distribution
	double depositot =0.;
	double tobedepositot =0.;
	double tobedeposi=0.;
	double depocorrection;
	double n_depocorrection;
	int depopoint =0;
	int smooth1st; 
	int smoothlast;
	double smoothdiff;
	vec deposi(Npx,0.0);
	//vec distribute(Npsl,0.0);            //needed for x-dependent distribution
	vec distribute;
	double alpha_lag1 = 0.;
		
	/* determine shear stress distribution over dune */
	bss2=sep_tau_distr(bss1);

	// op het reattachment point is de bodemschuifspanning gelijk aan de lokale tau_cr.
	// dit betekent dat er iets fout gaat op het moment dat 2 fsz's mergen en voorheen
	// de lokale bodemschuifspanning 0 was op xr.
	// daarom wordt er (tijdelijk) een correctie uitgevoerd op de bodemschuifspanning
	// door te checken of twee omliggende bss-waarden 0  zijn. De tussenliggende wordt
	// dan ook op 0 gezet
	for (int i=0;i<Npx;i++){
		if (bss2[o2(i-1)]==0 && bss2[o2(i+1)]==0) bss2[i]=0;
		//cerr<<"bodemschuifspanning op nul t.b.v. merging of fsz zones."<<endl;
	}
	for(int i=0;i<Npx;i++){bss2[i]/=S;} //makes it u_bed

	// rearrangement of fsz array based on xsi
	sep_sort_fsz(0);

	/* het uitvoeren van een bodem-update zonder migratie lijzijde */

	/* fluxtot is required for storage */
	for(int i=0;i<Npx;i++){
		fluxtot[i]=0.;
	}

	write_flowsep();
	

	//if(cfg.alpha_varies==0 && transport_eq == 2) { //needed for x-dependent distribution
	//    for(int j=0;j<Npsl;j++){ //needed for x-dependent distribution
	//		distribute[j] = exp(-dx*j/meanstle)/meanstle; //needed for x-dependent distribution
	//	} //needed for x-dependent distribution
	//} //needed for x-dependent distribution

	double ep=(1/(1-cfg.epsilonp));
	for(int t=0;t<int(cfg.tt);t++){  //deze loop zit er in omdat expliciet niet al te grote tijdstappen aankan, moet nog aan getweekt worden.

		detQcr(bss2, q_spec);

		int nfsz=(*fsz)[nf-1];
		int j=0;
		int xsi=(*fsz)[j*7+0];
		
		if(cfg.transport_eq == 1 || cfg.transport_eq == 3){ //OLAV: 2013 05 22
		for(int i=0;i<Npx;i++){
			/* see notes on details on this */
			if (i==xsi){ // OLAV 2013 04 16 (was just i==xsi check)
				double fluxxsi = (*flux)[o2(xsi+1)];
				(*b)[xsi]-=ep*dt/cfg.tt/dx*(fluxxsi-(*flux)[xsi])*2.;
				(*b)[o2(xsi+1)]=oldb[o2(xsi+1)];
				
				if (j<nfsz) {j++; i++; xsi=(*fsz)[j*7+0];} }
			else {
			(*b)[i]-=ep*dt/cfg.tt/dx*((*flux)[o2(i+1)]-(*flux)[i]);}
		}
		for(int i=0;i<Npx;i++){
			fluxtot[i]+=(*flux)[i]/cfg.tt*ep;
		}
		
		} //} //OLAV: 2011 02 28 end if (transport_eq == 1)

	else if(cfg.transport_eq == 2){
	
		//START OLAV 2014 01 30
		//Determine alpha and prepare distribution function if needed
		if(cfg.alpha_varies==0 && t==0){
			alpha_lag1=detAlphaLag(ub,cfg.alpha_varies,1);
			distribute=detDistributeFunc(alpha_lag1,dx);
		}

		if(cfg.alpha_varies!=0){
			//Determine alpha
			alpha_lag1=detAlphaLag(bss2,cfg.alpha_varies,t);
			
			//Determine distribute function
			distribute=detDistributeFunc(alpha_lag1,dx);
		} //end of calculation for alpha and distribution

		int Npsl = distribute.size();
		//END OLAV 2014 01 30
	
		for(int i=0;i<Npx;i++){
			depositot = (*flux)[i];    //to keep track of the running total of deposited material
			tobedepositot = depositot; //the goal value of the total deposited material
			
			while(depositot>0){
				for(int j=0;j<Npsl;j++){
					tobedeposi=distribute[j]*tobedepositot; //was tobedeposi=distribute[j]*dx*tobedepositot; OLAV 2014 03 17
					
					// if(i+j>Npx-1){
					// depopoint=i+j-Npx;}
					// }
					// else{depopoint=i+j;}
					
					if(i+j>Npx-1){
						depopoint=i+j;
						while (depopoint > Npx-1){
							depopoint=depopoint-Npx;
						}
					}
					else{depopoint=i+j;}
					
					if(tobedeposi < depositot){depositot=depositot-tobedeposi;}
					else{tobedeposi=depositot;depositot=0;}
					
					deposi[depopoint] += tobedeposi; //correct with i+j and o3/o2? 
				}
			} // closes while
		} //closes loop over x 
		
		///*
		// if (DebugOutput==1 && tijd > 2108){
	
		// ofstream outdebug;
		// outdebug.open ("out_debug1.txt", ofstream::out | ofstream::app);
		// outdebug.precision(6);
		
		// outdebug << tijd << " ";
		// for(int i=0;i<Npx;i++){outdebug << S*bss2[i] << " ";} 
		// outdebug << endl;
		// outdebug << tijd << " ";
		// for(int i=0;i<Npx;i++){outdebug << (*flux)[i] << " ";} 
		// outdebug << endl;
		// outdebug << tijd << " ";
		// for(int i=0;i<Npx;i++){outdebug << deposi[i] << " ";} 
		// outdebug << endl;
	
		// }//*/
		
		depocorrection = 0;
		n_depocorrection = 0;
		smooth1st =0;
		smoothlast =0;
		for(int i=0;i<Npx;i++){
			//if((*flux)[i]==0 && deposi[i] > 0){
			if((*flux)[i]==0 ){
			
				if( (*flux)[o2(i-1)] > 0) {smooth1st=i;}
				else if( (*flux)[o2(i+1)] > 0){smoothlast=i;}
								
				depocorrection += deposi[i];
				n_depocorrection +=1;
			}
		} 
		
		/*
		for(int i=0;i<Npx;i++){
			//if((*flux)[i]==0 && deposi[i] > 0){
			if((*flux)[i]==0 ){
				deposi[i]=depocorrection/n_depocorrection;
			}
		} 
		*/
		
		for(int i=0;i<Npx;i++){
			(*b)[i]-=ep*dt/cfg.tt*cfg.D50*((*flux)[i]-deposi[i]);
			//fluxtot[i]+=(*flux)[i]/cfg.tt*ep; //not the real flux! Is pick-up only!
			deposi[i]=0; 
		} //closes loop over x	
		
		if (cfg.DebugOutput==1 && tijd > 2079){
			ofstream outdebugbed;
			outdebugbed.open ("out_debugbed.txt", ofstream::out | ofstream::app);
			outdebugbed.precision(6);
		
			outdebugbed << tijd << " " << smooth1st << " " << smoothlast << " ";
			for(int i=0;i<Npx;i++){outdebugbed << (*b)[i] << " ";} 
			outdebugbed << endl;
		}
		
		// smoothing protocol
		/*
		smooth_param(5,1); // smooth at xsi
		smooth_param(5,2); // smooth at xri */
		
		/*
		int np = 5;
		int num = 1;
		
		vec oldbed=(*bp);
		int nfsz=(*fsz)[nf-1];
		int xi;

		for (int j=0;j<nfsz;j++) {

			if (num==1) { //smooth1st
				xi = smooth1st; }
			else if (num==2) { //xri
				xi = o3((*fsz)[j*7+1]); }

			for(int i=xi-3;i<=xi+3;i++) {

				int ii=i;
				if (i<0) ii+=Npx;
				if (i>Npx-1) ii-=Npx;

				int k=(np+1)/2-1;
				double sum=0.0;
				for(int j=-k;j<=k;j++){
					int m=ii+j;
					if (m<0) m+=Npx;
					else if (m>=Npx) m-=Npx;
					sum+=oldbed[m];}

				(*bp)[ii]=sum/double(np);
			}
		} */
		
		/*
		smoothdiff = (1./2.)*((*b)[smooth1st]-(*b)[o2(smooth1st-1)]);
		(*b)[o2(smooth1st-6)] += (1./4.)*smoothdiff;
		(*b)[o2(smooth1st-5)] += (1./4.)*smoothdiff;
		(*b)[o2(smooth1st-4)] += (1./4.)*smoothdiff;
		(*b)[o2(smooth1st-3)] += (1./4.)*smoothdiff;
		(*b)[o2(smooth1st-2)] += (1./4.)*smoothdiff;
		(*b)[o2(smooth1st-1)] += (1./4.)*smoothdiff;
		(*b)[smooth1st]       -= (3./4.)*smoothdiff;
		(*b)[o2(smooth1st+1)] -= (1./4.)*smoothdiff;
		(*b)[o2(smooth1st+2)] -= (1./4.)*smoothdiff;
		(*b)[o2(smooth1st+3)] -= (1./4.)*smoothdiff; */
		
		/*
		(*b)[o2(smooth1st-2)] += (1./4.)*smoothdiff;
		(*b)[o2(smooth1st-1)] += (3./4.)*smoothdiff;
		(*b)[smooth1st]       -= (3./4.)*smoothdiff;
		(*b)[o2(smooth1st+1)] -= (1./4.)*smoothdiff; */
		//(*b)[o2(smooth1st+2)] -= (1./2.)*smoothdiff;
		
		smoothdiff = (1./2.)*((*b)[smoothlast]-(*b)[o2(smoothlast+1)]);
		//(*b)[o2(smoothlast-2)] -= (1./4.)*smoothdiff;
		(*b)[o2(smoothlast-1)] -= (1./4.)*smoothdiff;
		(*b)[smoothlast]       -= (3./4.)*smoothdiff;
		(*b)[o2(smoothlast+1)] += (3./4.)*smoothdiff;
		(*b)[o2(smoothlast+2)] += (1./4.)*smoothdiff; //*/
		
		
		if (cfg.DebugOutput==1 && tijd > 2108){
			ofstream outdebugbed;
			outdebugbed.open ("out_debugbed.txt", ofstream::out | ofstream::app);
			outdebugbed.precision(6);
		
			outdebugbed << tijd << " " << smooth1st << " " << smoothlast << " ";
			for(int i=0;i<Npx;i++){outdebugbed << (*b)[i] << " ";} 
			outdebugbed << endl;
		}
		
	} //closes else if(transport_eq == 2){
	} // ends for(int t=0;t<int(cfg.tt);t++)

	for(int i=0;i<Npx;i++){bss2[i]*=S;}
	
	/* uitvoeren migratie van de lij-zijde */
	if (cfg.transport_eq == 1 || cfg.transport_eq == 3){
		sep_migr_lee(fluxtot,oldb);
	}
	
	//avalanching protocol
	if(cfg.transport_eq == 2 && cfg.AllowAvalanching == 1){
		avalanche(); 
	}
	
		// npasses=0;
		// npoints=0;
		// reposefound = 1;
		// while(reposefound==1){
					
			// reposefound =0;
			// npasses+=1;
			
			// for(int i=0;i<Npx;i++){
				// if (i < Npx-1) {
					// downstream = Npx-i-1;
					// upstream   = Npx-i-2;
					// }
				// else {
					// downstream = 0;
					// upstream   = Npx-1;
				// }
				// if(((*b)[downstream]-(*b)[upstream])/(dx) < reposeangle1){
					// npoints+=1;
					// reposefound=1;
					// deltab = (1./2.)*(reposeangleplus*dx+(*b)[upstream]-(*b)[downstream]); 
					// (*b)[upstream]-=deltab;
					// (*b)[downstream]+=deltab;					
				// }
			// }
		// }
		// cerr << "Avalanched in " << npasses << " pass(es) along the bed, adjusting " << npoints << " point(s)." << endl;
	// }

	newb=(*b);
	(*b)=oldb;

	// rearrangement of fsz array based on xdi for next time step
	sep_sort_fsz(4);

	
	write_flowsep();
	

	return newb;
}

/*
======================================================
======================================================

BLOCK VI: DETERMINATION OF THE BED SHEAR STRESS DISTRIBUTION
over a dune in case of flow separation.
input is bbs1 from flow solver, output is bbs
shear stress is set to zero inside the flow separation zone
parameterization bed shear stress distribution from reattachment point
until the next maximum of the bed shear stress.
NB!! still a buggy routine!!

======================================================
======================================================
*/

vec bottom::sep_tau_distr(vec tb){

	int nfsz=(*fsz)[nf-1];
	
	cout << endl << endl << "Block VI: bed shear stress" << endl << endl; //OLAV

	/* we gaan eerst xcin opnieuw bepalen, zijnde de maximale downstream bodemschuifspanning */
	/* NB this is a buggy piece of code!! */
	for (int j=0;j<nfsz;j++) {
		int xsi=(*fsz)[j*7+0];
		int xri=(*fsz)[j*7+1];
		int xti=(*fsz)[j*7+3];
		//dit doen vanaf de volgende trough, rekening houden met periodieke randvoorwaarden:
		if (xri<xsi) xri+=Npx;
		if (xti<xsi) xti+=Npx;
		int xcin = 0;
    int xs=xri;
		int nfp=1;
		vec ftau(Npx,0.0);
		//ftau=filter(nfp,tb);
		ftau=tb;
		//ftau=(*b);
		//cerr<<"xs: "<<xs;
		xs=o3(xs);
		//cerr<<" - "<<xs;
		if (ftau[xs]>ftau[o2(xs+1)]) xs=findTrough(xs,ftau);
		//cerr<<"xri: "<<xri<<"; xs: "<<xs<<";fbed: "<<fbed[o2(o2(xs-1)-1)]<<fbed[o2(xs-1)]<<fbed[xs]<<fbed[o2(xs+1)]<<fbed[o2(o2(xs+1)+1)]<<fbed[3]<<fbed[4]<<endl;
		//cerr<<" - "<<xs;
		for(int aa=xs;aa<xs+Npx;aa++){
			int m = o3(aa);
			if(ftau[m]>max(ftau[o2(m-1)],ftau[o2(m+1)])){
				xcin=m;
				//cerr<<j+1<<": xs: "<<xs<<"; xci_next: "<<xci_next<<endl;
				//cerr<<"bed: "<<(*bp)[16]<<" : "<<(*bp)[17]<<" : "<<(*bp)[18]<<" : "<<(*bp)[19]<<" : "<<(*bp)[20]<<" : "<<(*bp)[21]<<" : "<<(*bp)[22]<<endl;
				//cerr<<"m: "<<m<<endl;
				//cerr<<"tau m-2: "<<tb[m-2]<<" : "<<"tau m-1: "<<tb[m-1]<<" : "<<"tau m: "<<tb[m]<<" : "<<"tau m+1: "<<tb[m+1]<<" : "<<"tau m+2: "<<tb[m+2]<<" : "<<endl;
				//cerr<<" - xcin: "<<xcin<<endl;
				break;
			}
		}
		/*while ( max((*bp)[o2(xcin-1)],(*bp)[o2(xcin+1)]) > (*bp)[xcin]) {
			if      ((*bp)[xcin]>(*bp)[o2(xcin-1)]) {xcin++;} // left of top
			else if ((*bp)[xcin]>(*bp)[o2(xcin+1)]) {xcin--;} // right of top
			xcin=o3(xcin);
		}
		while ( (*bp)[xcin] == (*bp)[o2(xcin+1)] ) { xcin++; xcin=o3(xcin);}*/

		//cerr<<"bed: "<<(*bp)[xcin-2]<<" : "<<(*bp)[xcin-1]<<" : "<<(*bp)[xcin]<<" : "<<(*bp)[xcin+1]<<" : "<<(*bp)[xcin+2]<<" : "<<(*bp)[xcin+3]<<" : "<<(*bp)[xcin+4]<<" : "<<(*bp)[xcin+5]<<" : "<<(*bp)[xcin+6]<<endl;
		(*fsz)[j*7+5]=xcin;
	}

	/* overige deel met polynoom */
	for (int j=0;j<nfsz;j++) {
		int xsi=(*fsz)[j*7+0];
		int xri=(*fsz)[j*7+1];
		int xci=(*fsz)[j*7+2];
		int xti=(*fsz)[j*7+3];
		int xcin=(*fsz)[j*7+5];

		/* bepalen van polynoom om tau aan te passen om
		   numerieke problemen te voorkomen
		   deze polynoom past aan over de stoss-side */

//		cerr<<"xci_next: "<<xci_next<<endl;
//		cerr<<"xti: "<<xti<<endl;
		int npts;
		int xendi = xcin;
		if (xendi<xri){npts = Npx+xendi-(xri-1)+1;}
		else {npts = xendi-(xri-1)+1;}
		double Lt = (npts-1)*dx;

		double dhdx_x1=((*b)[o2(xri-1)]-(*b)[o2(o2(xri-1)-1)])/dx;
		double tauc_x1=cfg.thetacr*cfg.g*(2.65-1.)*cfg.D50*((1.+cfg.l2*dhdx_x1)/(pow(1.+dhdx_x1*dhdx_x1,(1./2.))));
		double tau_x1 = tauc_x1;
		double tau_x2 = tb[o2(xendi)];
		double dtaudx_x1 = 2*(tau_x2-tau_x1)/(Lt);
		double dtaudx_x2 = (tb[o2(xendi)]-tb[o2(xendi-1)])/(1.*dx);
		/*
		cerr<<"xendi: "<<xendi<<endl;
		cerr<<"x1: "<<(*x)[o2(xri-1)]<<endl;
		cerr<<"x2: "<<(*x)[o2(xendi)]<<endl;
		cerr<<"tau_x1: "<<tau_x1<<endl;
		cerr<<"tau_x2 [-1]: "<<tb[o2(xendi-1)]<<endl;
		cerr<<"tau_x2: "<<tau_x2<<endl;
		cerr<<"tau_x2 [+1]: "<<tb[o2(xendi+1)]<<endl;
		cerr<<"dtaudx_x1: "<<dtaudx_x1<<endl;
		cerr<<"dtaudx_x2: "<<dtaudx_x2<<endl;
		*/
		//cerr<<"tau_b ("<<j+1<<")"<<endl;
		//cerr<<"tb[xcin-1]: "<<tb[xcin-1]<<"tb[xcin]: "<<tb[xcin]<<endl;
		//cerr<<"dtaudx_x2: "<<dtaudx_x2<<endl;
		double p0 = tau_x1;
		double p1 = dtaudx_x1;
		//double p2 = (3.*tau_x2-L*(2.*dtaudx_x1+dtaudx_x2))/pow((L),2);
		//double p3 = -(-L*(dtaudx_x1+dtaudx_x2)+2.*tau_x2)/pow((L),3);
		double p2 = 3.*(tau_x2-tau_x1)/(Lt*Lt) - (dtaudx_x2 + 2.*dtaudx_x1)/Lt;
		double p3 = (dtaudx_x2 - 2.*p2*Lt - dtaudx_x1)/(3.*Lt*Lt);
		vec p(npts,0.0);
		int k,i;
		double xx=0.0;
		int is = xri-1;
		if (xri==0) is = Npx;
		for(k=0, i=is; k<npts; k++, i++){
			p[k]=p3*pow((xx),3)+p2*pow((xx),2)+p1*(xx)+p0;
			tb[o3(i)]=p[k];
			xx+=dx;
		}
		//cerr<<"p ("<<j+1<<")"<<endl;
		//cerr<<"p[end-1]: "<<p[npts-2]<<"p[end]: "<<p[npts-1]<<endl;
		//cerr<<"dpdx_x2: "<<(p[npts-1]-p[npts-2])/(1.*dx)<<endl;
		//cerr<<"Lt: "<<Lt<<"; xx: "<<xx-dx<<endl;
		 // cerr<<"p: "<<p<<endl;
	} // einde loop over flowsep zones


	/* naar nul zetten */
	for (int j=0;j<nfsz;j++) {
		int xsi=(*fsz)[j*7+0];
		int xri=(*fsz)[j*7+1];
		int end=xri-1;
		if (xri<xsi) end+=Npx;
		for(int i=xsi+1;i<=end;i++)	tb[o3(i)]=0.0;
	}

	for(int i=0;i<Npx;i++) tb[i]=max(tb[i],0.0);

	return tb;
}

/*
======================================================
======================================================

BLOCK VII: MIGRATION OF THE LEE-SIDE.
Volumetric sediment transport over the separation point
has to "fit" at the lee-side of a dune.
Consists of:
- sep_migr_lee()
- crossPoint_migrlee()
- area2D_Polygon()
See details on paper!
Includes one TODO!

======================================================
======================================================
*/

void bottom::sep_migr_lee(vec fluxtot, vec oldb){
	/* migration of the lee-side */

	int sepflag=(*fsz)[nf-2];
	int nfsz=(*fsz)[nf-1];
	int ntolreset=1;
	
//	//OLAV: added 2011/03/31
//	ofstream outdebug;
//    outdebug.open ("out_debug1.txt", ofstream::out | ofstream::app);
//    outdebug.precision(16);
//    //OLAV: added 2011/03/31

	for (int j=0;j<nfsz;j++) {
		cerr<<j+1<<": bepaling migratie lij-zijde per loslaatzone"<<endl;
		
		cout << endl << endl << "Block VII: migration lee side" << endl << endl;
		
		int	xsi=(*fsz)[j*7+0];
		int	xti=(*fsz)[j*7+3];
		double DH=(*b)[xsi]-(*b)[xti];
		cerr<<j+1<<": DH="<<DH<<endl;
		//double Ss=fluxtot[xsi]*dt-dt/dx*(fluxtot[xsi]-fluxtot[o2(xsi-1)]);
		double Ss=dt*(fluxtot[o2(xsi+1)]);
		Ss-=(*Sr)[j];
		cerr<<"fluxtot[xsi="<<o2(xsi+1)<<"]: "<<fluxtot[o2(xsi+1)]*dt<<"; fluxtot[xsi+1="<<o2(o2(xsi+1)+1)<<"]: "<<fluxtot[o2(o2(xsi+1)+1)]*dt<<"; Sround_prev"<<(*Sr)[j]<<endl;

		vec rp(2,0.0);
		double xsep=(*x)[xsi];
		double ysep=(*b)[xsi];
		double Npxcorr=pow(2.,round(log(double(Npx)/120.)/log(2.)));

		double step=pow(2.,-2.)*Npxcorr*dx;  // has to be multiple of dx!! ,e.g. dx/(2^3). or dx*2.
		double xdown=xsep;
		double itol=1e-6;
		double tol=itol;
		// double repose=-30.;
		// recalculate to r.c.'s
		double repose1=tan(cfg.repose);
		double b1 = repose1;
		double b2 = ysep;
		int dir = 1;  // 1 = left2right; 2 = right2left
		int max_it  = 100;
		int teller = 0;
		float area;
		double Sdif=100.;

		/*
		ofstream outpolyarea;
 		outpolyarea.open ("out_polyarea.txt", ofstream::out | ofstream::app);
		outpolyarea.precision(16);
		ofstream outproglee;
 		outproglee.open ("out_proglee.txt", ofstream::out | ofstream::app);
		outproglee.precision(16);
		*/

		vec b_prev = (*b);
		(*b)[xsi] = oldb[xsi];
		double xrr=DH/-repose1;

		if (fluxtot[o2(xsi+1)]>0. && Ss>0.) {

			while(fabs(Sdif)>tol && teller<max_it){
				xdown+=step;
				teller++;
				
				rp=crossPoint_migrlee(xdown,xdown+1.2*xrr,max_it,dir,b1,b2,tol,xsi,j+1);
				
				cout << endl << "Block VII: crosspoint migrlee exited" << endl << endl; //OLAV
				
				//cerr << 1 << " " ;// OLAV 2011-03-30
				
				int npts=(int(ceil((rp[0]-xsep)/step)+1.)*2-2)*2;
				vec xarr(npts,0.0);
				vec yarr(npts,0.0);
				vec nb(4,0.0);
				xarr[0]=xsep; yarr[0]=ysep;
				xarr[npts-1]=xsep;
				if (oldb[xsi]<ysep) {yarr[npts-1]=oldb[xsi];} // dit is dus de bodem uit de vorige tijdstap
				else {yarr[npts-1]=ysep;}
				double xval;
				double xloc=xsep;
				
				//cerr << 2 << " " ;// OLAV 2011-03-30
				
				
//				if (tijd >= tend){
//                         outdebug << endl << tijd << " Volume blabla "; //OLAV: added 2011/03/31
//                         cerr << 4 << " " ;// OLAV 2011-03-30
//                         }

				for(int i=1;i<npts-2;i=i+2){

					/* setup arrays for determination of volume */
					if (i<npts/2-2){
						xloc+=step;}
					else if (i==npts/2-1){
						xloc=rp[0];}
					else if (i==npts/2+1){
						xloc=xarr[i-3];}
					else if (i>npts/2+1){
						xloc-=step;
					} // if

					xarr[i]=xloc;
					xarr[i+1]=xloc;

					if (i<npts/2){
						if (xloc<=xdown) {
							yarr[i]=ysep;}
						else {
							yarr[i]=(xloc-xdown)*b1+b2;
						} // if
						xval=xarr[i+1];
   					    nb=paramFindNeighbors(xval,xsi);
						if (nb[1]<xval) {
							nb[0]+=Npx*dx;
							nb[1]+=Npx*dx;
						} // if
						yarr[i+1] = ((nb[1]-xval)*nb[2]+(xval-nb[0])*nb[3])/(nb[1]-nb[0]);
					} // if (i<npts/2)

					if (yarr[i]<yarr[i+1]) yarr[i]=yarr[i+1];
					else{
						if (xloc<=xdown) {
							yarr[i+1]=ysep;}
						else {
							yarr[i+1]=(xloc-xdown)*b1+b2;
						} // if
						xval=xarr[i];
   					nb=paramFindNeighbors(xval,xsi);
						if (nb[1]<xval) {
							nb[0]+=Npx*dx;
							nb[1]+=Npx*dx;
							}
						yarr[i] = ((nb[1]-xval)*nb[2]+(xval-nb[0])*nb[3])/(nb[1]-nb[0]);
					}
					if (yarr[i+1]<yarr[i]) yarr[i+1]=yarr[i];
				} // for
				
				//cerr << 3 << " " ;// OLAV 2011-03-30
				
				cout << endl << "Block VII: migrlee for finished--" << endl << endl; //OLAV

               /* determination of volume */
				area = area2D_Polygon( npts, xarr, yarr );
				
				cout << endl << "Block VII: migrlee area finished" << endl << endl; //OLAV
				
				Sdif=(area-Ss);
				//cerr<<"teller: "<<teller<<"; xdown = "<<xdown<<"; rp = "<<rp<<endl;
				cout << endl << "Block VII: migrlee Sdif finished" << endl << endl; //OLAV

				if (teller>max_it) {
					/* no appropriate volume found */
					tol=itol+ntolreset*itol;
					cerr<<"   WARNING: [j="<<j+1<<"] tolerance reset to "<<tol<<" (area="<<area<<")"<<endl;
					outlog<<"T="<<tijd<<" - WARNING: [j="<<j+1<<"] tolerance reset to "<<tol<<" (area="<<area<<")"<<endl;
					ntolreset+=1;
					xdown=xsep; step=dx; teller=1; Sdif=100.;
				}

				if (fabs(Sdif)<tol || teller==max_it){
					/* solution (appropriate volume) found */
					/*
					outpolyarea<<tijd<<" "<<npts<<" ";
					for(int i=0;i<xarr.size();i++) outpolyarea<<xarr[i]<<" ";
					for(int i=0;i<yarr.size();i++) outpolyarea<<yarr[i]<<" ";
					outpolyarea<<endl;
					outproglee<<tijd<<" "<<sepflag<<" "<<nfsz<<" "<<j<<" "<<teller<<" "<<Ss<<" "<<area<<" "<<Sdif<<" ";//endl;
					*/
				}
				else if (area>Ss){
					/* sediment volume to large, go step back, making "step" smaller */
					xdown=xdown-step;
					step/=2;
				}

                cout << endl << "Block VII: migrlee ifs are done" << endl << endl; //OLAV

				//if (tijd==690. && j+1==1){
				//	cerr<<"Ss = "<<Ss<<"; area = "<<area<<"; Sdif = "<<Sdif<<"; step="<<step<<"; npts="<<npts<<endl;
				//}

				/* sediment volume does not fit yet... continue */
			}   //while(fabs(Sdif)>tol)
			
			cout << endl << "Block VII: migrlee while is done" << endl << endl; //OLAV

			/* when we come here, the distance to migrate the lee-side is found */

			cerr<<"teller = "<<teller<<"; Ss = "<<Ss<<"; area = "<<area<<"; Sdif = "<<Sdif;
			(*b) = b_prev;

			/* eerst worden er een aantal karakteristieken van xdown bepaald */
			int xdi = 0;
			vec nb(4,0.0);
			double xdown_grid;
			nb = paramFindNeighbors(xdown,xsi);
			xdown_grid = nb[0];
			for (int i=0;i<Npx;i++){
				if ( (*x)[i]==xdown_grid ) {
					xdi = i;}
			}

			/* en vervolgens van het punt waar lij-zijde bodem raakt */
			int leerpi = 0;
			double leerp_grid;
			double leerp = rp[0];
			if (leerp>Npx*dx) leerp-=Npx*dx;
			nb = paramFindNeighbors(leerp,xdi);
			leerp_grid = nb[0];//dit moet altijd de linker buurman zijn
			for (int i=0;i<Npx;i++){
				if ( (*x)[i]==leerp_grid ) {
					leerpi = i;}
			}

			/* setting of the lee-side, interpolated to grid */
			double xx=xsep;
			if (leerpi<xsi) leerpi+=Npx;
			for (int i=xsi;i<=leerpi;i++){
				if (xx<=xdown) {(*b)[o3(i)]=ysep;}
				else {
					double val = (xx-xdown)*b1+b2;
					(*b)[o3(i)]=val;}
					//cerr<<leerpi<<" "<<i<<" "<<xx<<" "<<xdown<<" "<<(xx-xdown)<<" "<<(*b)[o3(i)]<<endl;
				xx+=dx;
			}

			/* due to numerical grid, rounding errors are made, which can explicitely be computed
			   these are separated in 3 areas: A,B&C. See notes on this
			   errors are corrected for in the next time step */
			double Sround = 0.0; double chiA = 0.0; double chiB = 0.0; double chiC = 0.0; double A = 0.0; double B = 0.0; double C = 0.0;
			double xd = xdown; if (xdown>=L) xd-=L;
			double rp0 = rp[0]; if (rp[0]>=L) rp0-=L;
			A -= (1./2.)*(dx+xd-((*x)[o3(xdi)]+dx))*(ysep-(*b)[o2(o3(xdi+1))]);
			chiB = ((*b)[o3(leerpi)]-(*b)[o3(leerpi+1)])*((*x)[o3(leerpi)]+dx-rp0)/dx - (rp[1]-(*b)[o3(leerpi+1)]);
			B = (1./2.)*chiB*dx;
			double xtemp = (*x)[o3(leerpi)];
			chiC = rp[1]- ( (xtemp+dx-rp0)*oldb[o3(leerpi)] + (rp0-xtemp)*oldb[o3(leerpi+1)] )/dx;
			C = (1./2.)*chiC*dx;
			Sround=B+A+C;

			/*
	 		xtemp = (*x)[o3(xdi)];
			double ytemp = (*b)[xdi];
			outproglee<<xtemp<<" "<<xd<<" "<<xtemp+dx<<" "<<ytemp<<" "<<ytemp<<" "<<(*b)[o3(xdi+1)]<<" "<<A<<" ";
			xtemp = (*x)[o3(leerpi)];
			ytemp = (*b)[o3(leerpi)];
			outproglee<<xtemp<<" "<<xtemp+dx<<" "<<rp0<<" "<<xtemp<<" "<<ytemp<<" "<<(*b)[o3(leerpi+1)]<<" "<<rp[1]<<" "<<ytemp<<" "<<B<<" ";
			ytemp = oldb[o3(leerpi)];
			outproglee<<xtemp<<" "<<rp0<<" "<<xtemp+dx<<" "<<xtemp<<" "<<ytemp<<" "<<rp[1]<<" "<<(*b)[o3(leerpi+1)]<<" "<<ytemp<<" "<<C<<" ";
			*/

			cerr<<"; Sround + Sdif = "<<Sround+Sdif<<endl;

			/* store some parameters */
			//outproglee<<Sround<<endl;	outproglee.close();
			(*fsz)[j*7+4]=xdi;
			(*Sr)[j]=Sround+Sdif;

		} //if (fluxtot[o2(xsi+1)]>0.)

		else {
			/* this is the case when tau is below critical over entire stoss-side */
			(*fsz)[j*7+4]=(*fsz)[j*7+0]; // op xsi zetten, niets doen
			(*fsz)[j*7+6]=5; // case 5 identifier
			if (fluxtot[o2(xsi+1)]<=0.) {
				cerr<<"   WARNING: [j="<<j+1<<"] tau below critical over entire stoss-side; nothing done."<<endl;
				outlog<<"T="<<tijd<<" - WARNING: [j="<<j+1<<"] tau below critical over entire stoss-side; nothing done."<<endl;
			}
			else if (Ss < 0.) {
				cerr<<"   WARNING: [j="<<j+1<<"] Ss < 0 due to low sediment transport; nothing done."<<endl;
				outlog<<"T="<<tijd<<" - WARNING: [j="<<j+1<<"] Ss < 0 due to low sediment transport; nothing done."<<endl;
			}
		}

	} //loop over flowsep zones

    //outdebug.close(); //OLAV: added 2011/03/31


	/* dit doen we liever in main.cpp.
	ofstream outSround;
	outSround.open ("out_Sround.txt", ofstream::out | ofstream::app);
	outSround.precision(16);
	outSround<<tijd<<" "; for(int i=0;i<(*Sr).size();i++)outSround<<(*Sr)[i]<<" "; outSround<<endl;
	outSround.close();
	*/
}

vec bottom::crossPoint_migrlee(double xl_in, double xr_in, int max_it, int dir, double b1, double b2, double tol, int xi, int j){
  vec rp(2,0.0);
	vec nb(4,0.0);
	dir=1; max_it=150;
	double x_p = 0.0;
	double xl=xl_in;
	double xr=xr_in;
	double x_p_bed;
	double fac=1;
	int pr_dir=dir;

	//cout << endl << endl << "Block VII: crosspoint migrlee start" << endl << endl; //OLAV


  for(int k=1;k<=max_it;k++){
		// 1 = left2right; 2 = right2left

		/* toegevoegd op 23/4/2007 */
    if (k>1 && x_p<=xl_in) {
				//cerr<<"aap"<<endl;
				x_p=xr_in;
				dir=1;
				pr_dir=dir;
				fac*=2;
				}

		if (dir==2){
     	x_p = xr - dx/fac;}
		else{
     	x_p = xl + dx/fac;
    } // if
    // find neighboring points of the bed
    nb=paramFindNeighbors(x_p,xi);
		x_p_bed = x_p;
		if (x_p>=Npx*dx || x_p==L) x_p_bed = x_p-Npx*dx;
		double y1 = ((nb[1]-x_p_bed)*nb[2]+(x_p_bed-nb[0])*nb[3])/(nb[1]-nb[0]);
    double xx=x_p-xl_in;
    double y2 = b1*xx + b2; // height of -30 degrees line

		if (nb[2]==0 && nb[3]==0) {
			cerr<<"   WARNING: [j="<<j<<"] neighbor points not found correctly!:"<<endl;
			outlog<<"T="<<tijd<<" - WARNING: [j="<<j<<"] neighbor points not found correctly!:"<<endl;
			cerr<<"x_p="<<x_p<<"; xi="<<xi<<endl;
			cerr<<"Npx*dx="<<Npx*dx<<endl;
			cerr<<nb[0]<<" "<<nb[1]<<" "<<nb[2]<<" "<<nb[3]<<endl;
			cerr<<k<<" "<<dir<<" "<<xl_in<<" "<<x_p<<" "<<y1<<" "<<y2<<" "<<xx<<" "<<abs(y2-y1)<<endl;
		}


   	if (abs(y2-y1) <= tol || abs(y2-y1) == 0 || k==max_it){

     	rp[0]=x_p;
     	rp[1]=y2;

     	if (k==max_it) {
				/* no point found */
				cerr<<"   WARNING: [j="<<j<<"] no cross-point found to determine migration lee-side! restarting by exiting with xl...!"<<endl;
				outlog<<"T="<<tijd<<" - WARNING: [j="<<j<<"] no cross-point found to determine migration lee-side! restarting by exiting with xl...!"<<endl;
				rp[0]=xl_in;
				rp[1]=-999.;
				cerr<<"xl_in: "<<xl_in<<"xr_in: "<<xr_in<<endl;
				cerr<<nb[0]<<" "<<nb[1]<<" "<<nb[2]<<" "<<nb[3]<<endl;
				cerr<<k<<" "<<dir<<" "<<xl_in<<" "<<x_p<<" "<<y1<<" "<<y2<<" "<<xx<<" "<<abs(y2-y1)<<endl;
				break;
      }
     	else {
				break;
			}

		}
   	else{
			/* next iteration */
      if (y2 > y1){
       	dir     = 1;  // left2right
       	xl  = x_p;}
      else{
       	dir     = 2;  // right2left
       	xr = x_p;}
      if (dir!=pr_dir) fac*=2; // for interation to solution
			pr_dir=dir;
		}
		//if ((tijd==690. && j==1)){
		//		cerr<<k<<" "<<dir<<" "<<xl_in<<" "<<x_p<<" "<<y1<<" "<<y2<<" "<<xx<<" "<<abs(y2-y1)<<" "<<nb[0]<<" "<<nb[1]<<" "<<nb[2]<<" "<<nb[3]<<endl;
		//}
	} // for

  	cout << endl << "Block VII: crosspoint migrlee end" << endl << endl; //OLAV

	return rp;
}

double bottom::area2D_Polygon(int n, vec xarr, vec yarr ){
// area2D_Polygon(): computes the area of a 2D polygon
// see: http://softsurfer.com/Archive/algorithm_0101/algorithm_0101.htm#area2D_polygon()
    double area = 0;
    int   i, j, k;     // indices

    cout << endl << "Block VII: area started" << endl << endl; //OLAV

    for (int i=1, j=2, k=0; i<n; i=i+2, j=j+2, k=k+2) {
		if (j==n) j=0;
        area += fabs(xarr[i]-xarr[k]) * fabs(yarr[j] - yarr[i]);
    }
    return area / 2.0;
}

/*
======================================================
======================================================

BLOCK VIII: DETERMINATION OF DUNE CHARACTERISTICS:
number of dunes, average height and average length
we have to methods available: one with fft
NB both don't work correctly yet.

======================================================
======================================================
*/

vec bottom::detNd(vec bot){
	/* without fast fourier transform */

	vec Dc(3,0.0);
	// create temporary arrays to store positions
	vector<int> tpos_temp(Npx/2,0); //NB: Npx moet even zijn
	vector<int> cpos_temp(Npx/2,0); //NB: Npx moet even zijn
	vec botf(Npx,0.0);
	botf=filter(9,bot);

	//ofstream outtemp("out_temp.txt");
	//outtemp.precision(10);
	//outtemp<<botf<<endl;
	
	cout << endl << endl << "Block VIII: dune char" << endl << endl; //OLAV

	int Nd = 0; int stop=0; int pos=Npx-1;
	int pos_t1 = 0; // position of first trough from right.
	for (int m=Npx-1; m>=0; m--){
		if (botf[m]<min(botf[o2(m-1)],botf[o2(m+1)])) {pos_t1=m; break;}
	}
	//cerr<<"pos first trough: "<<pos_t1<<" - "<<endl;
	while(stop==0) {
		//find Trough
		for (int i=pos; i>pos-Npx; i--){
			int m=i; if (i<0) m=i+Npx;
			//if (botf[m]<min(botf[o2(m-1)],botf[o2(m+1)])) {pos=m;
			if (botf[m]<min(botf[o2(m-1)],botf[o2(m+1)]) && botf[m]!=botf[o2(m-1)]) {
				pos=m;
				break;
			}
		}
		//cerr<<"pos trough: "<<pos;
		tpos_temp[Nd]=pos;
		if (pos==pos_t1 && Nd>0) {stop = 1; break; cerr<<endl;}
		//find crest
		for (int i=pos; i>pos-Npx; i--){
			int m=i; if (i<0) m=i+Npx;
			//if (botf[m]>max(botf[o2(m-1)],botf[o2(m+1)]) | (botf[m]==botf[o2(m-1)] && botf[m]==botf[o2(m+1)]) ) {
			if (botf[m]>max(botf[o2(m-1)],botf[o2(m+1)]) && botf[m]!=botf[o2(m-1)] ) {
				pos=m;
				break;
			}
		}
		//cerr<<"pos crest: "<<pos<<endl;
		cpos_temp[Nd]=pos;
		Nd++;
		if (Nd>Npx) {
			cerr<<"   WARNING: detNd gaat nog niet helemaal ok"<<endl; Nd=Npx;
			outlog<<"T="<<tijd<<" - WARNING: detNd gaat nog niet helemaal ok"<<endl; Nd=Npx;
			break;
		}
	}

	// store real positions in array of correct dimension
	vector<int> tpos(Nd,0);
	vector<int> cpos(Nd,0);
	for (int i=0; i<Nd; i++){
		tpos[i]=tpos_temp[i];
		cpos[i]=cpos_temp[i];
	}
	//cerr<<"troughs: "; for (int i=0; i<Nd; i++){cerr<<tpos[i]<<" ";} cerr<<endl;
	//cerr<<"crests:  "; for (int i=0; i<Nd; i++){cerr<<cpos[i]<<" ";} cerr<<endl;

	//store La and Ha arrays
	vec Ha(Nd,0.0);
	vec La(Nd,0.0);
	for (int i=0; i<Nd; i++){
		int tposi=tpos[i];
		int cposi=cpos[i];
		Ha[i]=(*b)[cposi]-(*b)[tposi];
		int il=i+1; if(i==Nd-1) il=0;
		int tposil=tpos[il];
		if (tposi<tposil) tposi+=Npx;
		//cerr<<tposil<<" "<<tposi<<endl;
		double val=double(tposi-tposil)*dx;
		if (val>0) {La[i]=val;} else {La[i]=L;}
	}
	//cerr<<"La: "; for (int i=0; i<Nd; i++){cerr<<La[i]<<" ";} cerr<<endl;
	//cerr<<"Ha: "; for (int i=0; i<Nd; i++){cerr<<Ha[i]<<" ";} cerr<<endl;
	double Lav=meanval(La,Nd);
	double Hav=meanval(Ha,Nd);
	//cerr<<"Lav: "<<Lav<<endl;
	//cerr<<"detNd Hav: "<<Hav<<endl; // OLAV 2014 03 31
	Dc[0]=Nd; Dc[1]=Lav; Dc[2]=Hav;

	return Dc;
}

vec bottom::detNd_fft(vec bot, int fftnum){
	/* dune characteristics using FFT */

	// create temporary arrays to store positions
	vec Dc(3,0.0);
	vec botf(Npx,0.0);
	vector<int> tpos_temp(Npx/2,0); //NB: Npx moet even zijn
	vector<int> cpos_temp(Npx/2,0); //NB: Npx moet even zijn
	
	cout << endl << endl << "Block VIII: dune char FFT" << endl << endl; //OLAV

	botf=fftBed(bot,2);
	//botf=filter(31,botf);

	/*
	ofstream outtemp;
	outtemp.open ("out_temp.txt", ofstream::out | ofstream::app);
	outtemp.precision(6);
	for(int i=0;i<Npx;i++) outtemp<< bot[i]<<" "; outtemp<<endl;
	for(int i=0;i<Npx;i++) outtemp<<botf[i]<<" "; outtemp<<endl;
	//outtemp.close; 
	*/

	int Nd = 0; int stop=0; int pos=Npx-1;
	int pos_t1=-1; // position of first trough from right.
	
	//BEGIN OLAV 2014 03 31
	double testval=1.e99;
	int testm=-1;
	//END OLAV 2014 03 31
	
	//BEGIN OLAV 2014 03 31
	if (cfg.nd==1){
		for (int m=0; m<Npx; m++){
			if (botf[m]<testval){ //niet met botf[m]<=testval? Nu wordt het meest linkse laagste punt geselecteerd
				testval=botf[m];
				testm=m;
			}
		}
		pos_t1=testm; 
	}
	else{ //OLAV 2014 03 31 this is the original part
		for (int m=Npx-1; m>=0; m--){ 
			if (botf[m]<min(botf[o2(m-1)],botf[o2(m+1)])) {pos_t1=m; break;}
		}
	}
	//END OLAV 2014 03 31 //WAS: for (int m=Npx-1; m>=0; m--){if (botf[m]<min(botf[o2(m-1)],botf[o2(m+1)])) {pos_t1=m; break;}}
	
	//cerr<<"pos first trough: "<<pos_t1<<" - "<<endl;
	while(stop==0) { 
		if (cfg.nd==1){ //OLAV 2014 03 31
			stop =1;
			tpos_temp[Nd]=pos_t1;
			cerr<<endl;
			
			testval=-1.e99;
			testm=-1;
			
			//find crest
			for (int m=0; m<Npx; m++){
				if (botf[m]>testval){ //niet met botf[m]=>testval? Nu wordt het meest linkse hoogste punt geselecteerd
					testval=botf[m];
					testm=m;
				}
			}
			pos=testm;
		}//OLAV 2014 03 31
		else{//OLAV 2014 03 31 this is the original part
		//find Trough
			for (int i=pos; i>pos-Npx; i--){
				int m=i; if (i<0) m=i+Npx;
				//if (botf[m]<min(botf[o2(m-1)],botf[o2(m+1)])) {pos=m;
				if (botf[m]<min(botf[o2(m-1)],botf[o2(m+1)]) && botf[m]!=botf[o2(m-1)]) {
					pos=m;
					break;
				}
			}
			//cerr<<"pos trough: "<<pos;
			tpos_temp[Nd]=pos;
			if (pos==pos_t1 && Nd>0) {stop = 1; break; cerr<<endl;}
			//find crest
			for (int i=pos; i>pos-Npx; i--){
				int m=i; if (i<0) m=i+Npx;
				//if (botf[m]>max(botf[o2(m-1)],botf[o2(m+1)]) | (botf[m]==botf[o2(m-1)] && botf[m]==botf[o2(m+1)]) ) {
				if (botf[m]>max(botf[o2(m-1)],botf[o2(m+1)]) && botf[m]!=botf[o2(m-1)] ) {
					pos=m;
					break;
				}
			}
		}//OLAV 2014 03 31
		//cerr<<"pos crest: "<<pos<<endl;
		cpos_temp[Nd]=pos;
		Nd++;
		if (Nd>Npx) {
			cerr<<"   WARNING: detNd gaat nog niet helemaal ok"<<endl; Nd=Npx;
			outlog<<"T="<<tijd<<" - WARNING: detNd gaat nog niet helemaal ok"<<endl; Nd=Npx;
			break;
		}
	}

	// store real positions in array of correct dimension
	vector<int> tpos(Nd,0);
	vector<int> cpos(Nd,0);
	for (int i=0; i<Nd; i++){
		tpos[i]=tpos_temp[i];
		cpos[i]=cpos_temp[i];
	}
	//cerr<<"troughs: "; for (int i=0; i<Nd; i++){cerr<<tpos[i]<<" ";} cerr<<endl;
	//cerr<<"crests:  "; for (int i=0; i<Nd; i++){cerr<<cpos[i]<<" ";} cerr<<endl;

	/* The real tops and troughs are between the previously found estimates
	 * using fft transforms
	 */

	for (int i=0; i<Nd; i++){
		int tposi=tpos[i];
		int cposi=cpos[i];
		if (tposi<cposi) tposi+=Npx;
		double minval = 100.; double maxval = -100.;
		for (int j=cposi; j<=tposi; j++) {
			int m=o3(j);
			if ( (*b)[m]>maxval ) {
				maxval=(*b)[m];
				cpos[i]=m; }
			else if ( (*b)[m]<minval ) {
				minval=(*b)[m];
				tpos[i]=m; }
		}
	}

	//store La and Ha arrays
	vec Ha(Nd,0.0);
	vec La(Nd,0.0);

	for (int i=0; i<Nd; i++){
		int tposi=tpos[i];
		int cposi=cpos[i];
		Ha[i]=(*b)[cposi]-(*b)[tposi];
		int il=i+1; if(i==Nd-1) il=0;
		int tposil=tpos[il];
		if (tposi<tposil) tposi+=Npx;
		//cerr<<tposil<<" "<<tposi<<endl;
		double val=double(tposi-tposil)*dx;
		if (val>0) {La[i]=val;} else {La[i]=L;}
	}

	//cerr<<"La: "; for (int i=0; i<Nd; i++){cerr<<La[i]<<" ";} cerr<<endl;
	//cerr<<"Ha: "; for (int i=0; i<Nd; i++){cerr<<Ha[i]<<" ";} cerr<<endl;
	double Lav=meanval(La,Nd);
	double Hav=meanval(Ha,Nd); 
	
	//cerr<<"detNd_fft Hav: "<<Hav<<endl; // OLAV 2014 03 31
	
	//cerr<<"Lav: "<<Lav<<endl;
	//cerr<<"Hav: "<<Hav<<endl;
	Dc[0]=Nd; //Dc[1]=Lav; Dc[2]=Hav;
	
	//cerr<<"detNd_fft Dc0: "<<Dc[0]<<" detNd_fft Dc1: "<<Dc[1]<<" detNd_fft Dc2: "<<Dc[2]<<endl; // OLAV 2014 03 31

	Dc[1]=(*b)[cpos[0]];
	Dc[2]=(*b)[tpos[0]];
	
	//cerr<<"detNd_fft Dc0: "<<Dc[0]<<" detNd_fft Dc1: "<<Dc[1]<<" detNd_fft Dc2: "<<Dc[2]<<endl; // OLAV 2014 03 31

	return Dc;

}

vec bottom::fftBed(vec bed, int fftnum) {

	vec bedfft(Npx,0.0);

	fft dft(Npx);


cout << endl << endl << "Block VIII: dune char fft bed" << endl << endl; //OLAV

	//for(int i=0;i<Npx;i++) outlog<<bed[i]<<" "; outlog<<endl;

	vector<complex<double> > d(Npx/2+1);
	for(int i=0;i<Npx/2+1;i++)d[i]=i*2.*M_PI/L*complex<double>(0,1);
	vector<complex<double> > du(Npx/2+1);

	dft.heen(bed,du);

	vector<int> max_comp(fftnum,-1);
	double max_val = 0.0;
	int max_loc = -1;
	double max_val_prev = 0.0;
	for (int k=0;k<fftnum;k++) {
		for(int i=0;i<Npx/2+1;i++) {
			double val = abs(du[i]);
			if (k==0 && val > max_val) {
				max_loc = i; max_val = val;
			}
			else if (k>0 && val > max_val && val < max_val_prev) {
				max_loc = i; max_val = val;
			}
		}
		max_val_prev=max_val;
		max_val=0.0;
		max_comp[k]=max_loc;
	}

	//for (int k=0;k<fftnum;k++) cerr<<max_comp[k]<<" "; cerr<<endl;

	vector<complex<double> > du_new(Npx/2+1);

	for (int k=0;k<fftnum;k++) {
		int loc = max_comp[k];
		du_new[loc]=du[loc];
	}

	dft.terug(du_new,bedfft);

	return bedfft;

}

/*
======================================================
======================================================

BLOCK IX: DETERMINATION OF DUNE MIGRATION AND GROWTH.
using FFT

======================================================
======================================================
*/

//TODO maak bottom::detGrow(vec current, vec next), om daarmee de groei te bepalen in plaats van met het zoeken van het lokale maximum
//TODO schrijf bij de laatste tijdstap van het model ook de de abs (amplitude) en het argument (fase) van de fourier weg.

/*
double bottom::detGrow(vec current, vec next){
	double grow = 0;
	vector<complex<double> > du0 = fftbot(current);

}
*/

double bottom::detMigr(vec current,vec next) {

	double mig = 0;

	vector<complex<double> > du0 = fftbot(current);
	int loc0 = maxloc_complex(du0);
	double k = 2.*M_PI*double(loc0)/L;
	double a0 = real(du0[loc0]);
	double b0 = imag(du0[loc0]);
	double phi0 = atan2(-b0,a0);
//#define SHOW_SPECTRUM
#ifdef SHOW_SPECTRUM
	vec amp(du0.size());
	for (auto i = 0; i < du0.size(); ++i) {
		amp[i] = abs(du0[i]);
		//cout << i << " " << amp[i] << " " << arg(du0[i]) <<endl;
	}
#endif

	vector<complex<double> > du1 = fftbot(next);
	int loc1 = maxloc_complex(du1);
	double a1 = real(du1[loc1]);
	double b1 = imag(du1[loc1]);
	double phi1 = atan2(-b1,a1);
#ifdef SHOW_SPECTRUM
	for (auto i = 0; i < next.size(); ++i){
		cout << i << " "  << next[i] << endl;
	}
	for (auto i = 0; i < du0.size(); ++i) {
		amp[i] = abs(du1[i]);
		if (amp[i] > 1e-8)
			cout << i << " " << amp[i] << " " << arg(du1[i]) <<endl;
	}
#endif

	double phidif = phi1-phi0;
    if (phi1<0 && phi0>0) {phidif = phidif + 2.*M_PI;}
	//cerr<<"phidif: "<<phidif<<endl;

	mig = phidif/k/dt;

	return mig;
}

vector<complex<double> > bottom::fftbot(vec bed) {

	fft dft(Npx);

	vector<complex<double> > d(Npx/2+1);
	for(int i=0;i<Npx/2+1;i++)d[i]=i*2.*M_PI/L*complex<double>(0,1);
	vector<complex<double> > du(Npx/2+1);

	dft.heen(bed,du);

	return du;

}

int bottom::maxloc_complex(vector<complex<double> > du) {
	double max_val = 0.0;
	int max_loc = -1;
	for(int i=0;i<Npx/2+1;i++) {
		double val = abs(du[i]);
		if (val > max_val) {
			max_loc = i;
			max_val = val;
		}
	}
	return max_loc;
}

/*
======================================================
======================================================

BLOCK X: REMAINING ROUTINES

======================================================
======================================================
*/ 

double bottom::detAlphaLag(vec ub, int method,int suppressoutput){

	double alpha_lag1=0.;
	double alpha_lag_temp=0.;
	double u_star = 0.;
	double u_star_temp = 0.;
	double theta  = 0.;
	double theta_temp  = 0.; 

	int skipped = 0;
	
	//determine alpha
	if(method == 0){ // constant alpha
		alpha_lag1=cfg.alpha_lag;
	}
	else if(method == 1){ // Sekine & Kikkawa 
		for(int i=0;i<Npx;i++) {
			//S*ub=volumetric bed shear stress; 
			//u_star=(tau_vol)^(1/2)
			u_star_temp= pow(S*(1./2.)*(ub[o2(i-1)]+ub[i]),(1./2.)); 
			if (u_star_temp>cfg.u_star_cr){
				u_star = u_star+u_star_temp; //u_star = u_star+(1/double(Npx))*u_star_temp;
			}
			else {
				skipped+=1;
			}
		}	
		
		u_star/=double(Npx-skipped);
		
		alpha_lag1 = (cfg.alpha_2*pow(u_star/cfg.w_s,(3./2.))*(1-(cfg.u_star_cr/cfg.w_s)/(u_star/cfg.w_s)));
		alpha_lag1 = max(alpha_lag1,cfg.alpha_min_SK);
		alpha_lag1 = min(alpha_lag1,cfg.alpha_max_SK);
		
		if(suppressoutput==0){
			cerr << "alpha_lag1=" << alpha_lag1<< " D_star=" << cfg.D_star<< " u_star=" << u_star << " u_star_cr=" << cfg.u_star_cr << " w_s=" << cfg.w_s << endl;
		}
	}			
	else if(method == 2){ // Shimizu et al. original model
		for(int i=0;i<Npx;i++) {
			// S*ub=volumetric bed shear stress;
			// S*ub*ro = bed shear stress
			theta_temp = S*(1./2.)*(ub[o2(i-1)]+ub[i])/(cfg.delta*cfg.g*cfg.D50);
			
			if(theta_temp>cfg.thetacr){ //2014 02 04 theta_temp>0
				theta = theta+theta_temp; //theta = theta+(1/double(Npx))*theta_temp;
			}
			else{
				skipped+=1;
			}
		}
		
		theta/=double(Npx-skipped);
		
		//theta = H*ii / (cfg.D50*cfg.delta);	//Change Olav 2015 02 28
		//theta_temp = theta; //Change Olav 2015 02 28 
		//theta=0.06+0.3*pow(theta_temp,(3/2)); //Change Olav 2015 02 28 
		
		if(theta>cfg.theta_max_S){
			alpha_lag1=cfg.alpha_max_S;

		}
		else if(theta<cfg.theta_min_S)	{
			alpha_lag1=cfg.alpha_min_S;
		}
		else{
			alpha_lag1=cfg.alpha_min_S+(theta-cfg.theta_min_S)*(cfg.alpha_max_S-cfg.alpha_min_S)/(cfg.theta_max_S-cfg.theta_min_S);
		} 	 
		
		if(suppressoutput==0){cerr << "alpha IS ADJUSTED from " << alpha_lag1 << endl;} //Change Olav 2015 02 28 
		
//		alpha_lag1 = alpha_lag1 / cfg.H_ref * H; //not in orignial Shimizu model only in changed model Change Olav 2015 02 28
		
		if(suppressoutput==0){cerr << "alpha IS ADJUSTED to " << alpha_lag1 << endl;} //Change Olav 2015 02 28 
		
		if(suppressoutput==0){
			//cerr << "THETA IS ADJUSTED from " << theta_temp << " --> " << theta << endl; //Change Olav 2015 02 28 
			cerr << "alpha_lag1=" << alpha_lag1<< " theta=" << theta << endl;
		}
		
	}
	
	else if(method == 3){ // Shimizu et al. adjusted to van Duin 2021
		for(int i=0;i<Npx;i++) {
			theta_temp = S*(1./2.)*(ub[o2(i-1)]+ub[i])/(cfg.delta*cfg.g*cfg.D50);

			if(theta_temp>cfg.thetacr){ //2014 02 04 theta_temp>0
				theta = theta+theta_temp; //theta = theta+(1/double(Npx))*theta_temp;
			}
			else{
				skipped+=1;
			}
		}

		theta/=double(Npx-skipped);

		if(theta<cfg.theta_min_S)	{
			alpha_lag1=cfg.alpha_min_S;
		}
		else{
			alpha_lag1=cfg.alpha_min_S+(theta-cfg.theta_min_S)*(cfg.alpha_max_S-cfg.alpha_min_S)/(cfg.theta_max_S-cfg.theta_min_S);
		}

		if(suppressoutput==0){cerr << "alpha IS ADJUSTED from " << alpha_lag1 << endl;} //Change Olav 2015 02 28

		alpha_lag1 = alpha_lag1 / cfg.H_ref * H; //Change Olav 2015 02 28

		if(suppressoutput==0){cerr << "alpha IS ADJUSTED to " << alpha_lag1 << endl;} //Change Olav 2015 02 28

		if(suppressoutput==0){
			cerr << "alpha_lag1=" << alpha_lag1<< " theta=" << theta << endl;
		}
	}

	// double reductionfactor = 0.99;
	// if (5.*alpha_lag1*cfg.D50 > L) {
		// if(suppressoutput==0) {cerr << "WARNING: alpha too large with: " << alpha_lag1 << ". Set to: " << reductionfactor*L/5./cfg.D50 << endl;}
		
		// alpha_lag1 = reductionfactor*L/5./cfg.D50;
	// }
	
	if(suppressoutput==0){
		ofstream outdebug;
		outdebug.open ("out_debug1.txt", ofstream::out | ofstream::app);
		outdebug.precision(6);
		
		outdebug << tijd << " " << theta << " " << u_star << " " << alpha_lag1 << endl;
	}
	
	return alpha_lag1; 
}

vec bottom::detDistributeFunc(double alpha_lag1,double deltax){	
	double meanstle1  = alpha_lag1*cfg.D50;
	int Npsl = (int)ceil(meanstle1*cfg.stle_factor/deltax);
	vec distribute(Npsl,0.0);
	for(int j=0;j<Npsl;j++){ 
		distribute[j]= -exp(-deltax*(j+1)/meanstle1) + exp(-deltax*j/meanstle1);
	} 
	
	return distribute;
}

// vec bottom::detDistributeFunc(double alpha_lag1,double deltax){	
	// double meanstle1  = alpha_lag1*cfg.D50;
	// int Npsl = (int)ceil(meanstle1*stle_factor/deltax); 
	// vec distribute(Npsl,0.0);
	// for(int j=0;j<Npsl;j++){ 
		// distribute[j] = deltax*exp(-deltax*j/meanstle1)/meanstle1; //was distribute[j] = exp(-deltax*j/meanstle1)/meanstle1; OLAV 2014 03 17
	// } 
	
	// return distribute;
// }

void bottom::avalanche(){
	double reposeangle1 = tan(cfg.repose); //double sepcritangle1 = tan(sepcritangle); //sepcritangle=-10.;
	double reposeangleplus = 0.99*reposeangle1;
	double deltab; 
	int reposefound;
	int npasses;
	int npoints;
	int downstream;
	int upstream;
	
	npasses=0;
	npoints=0;
	reposefound = 1;
	while(reposefound==1){
				
		reposefound =0;
		npasses+=1;
		
		for(int i=0;i<Npx;i++){
			if (i < Npx-1) {
				downstream = Npx-i-1;
				upstream   = Npx-i-2;
				}
			else {
				downstream = 0;
				upstream   = Npx-1;
			}
			if(((*b)[downstream]-(*b)[upstream])/(dx) < reposeangle1){
				npoints+=1;
				reposefound=1;
				deltab = (1./2.)*(reposeangleplus*dx+(*b)[upstream]-(*b)[downstream]); 
				(*b)[upstream]-=deltab;
				(*b)[downstream]+=deltab;					
				}
		}
	}
	cerr << "Avalanched in " << npasses << " pass(es) along the bed, adjusting " << npoints << " point(s)." << endl;

}

vec bottom::smooth(vec bed_in){
	vec bed_out(Npx,0.0);
	for(int i=0;i<Npx;i++) {
		bed_out[i]=( bed_in[o2(i-1)]+bed_in[i]+bed_in[o2(i+1)] )/3.;
	}
	return bed_out;
}

vec bottom::filter(int np, vec inp_arr){
	vec out_arr(Npx,0.0);
	for(int i=0;i<Npx;i++) {
		int k=(np+1)/2-1;
		double sum=0.0;
		for(int j=-k;j<=k;j++){
			int m=i+j;
			if (m<0) m+=Npx;
			else if (m>=Npx) m-=Npx;
			sum+=inp_arr[m];}
		out_arr[i]=sum/double(np);
	}
	return out_arr;
}

void bottom::smooth_param(int np, int num){
	/* smoothen bij het reattachment point
	 * np geeft het aantal punten waarover gemiddeld wordt
	 * num=1: smooth at xsi; num=2: smooth at xri */

	vec oldbed=(*bp);
	int nfsz=(*fsz)[nf-1];
	int xi = 0;

	for (int j=0;j<nfsz;j++) {

		if (num==1) { //xsi
			xi = (*fsz)[j*7+0]; }
		else if (num==2) { //xri
			xi = o3((*fsz)[j*7+1]); }

		for(int i=xi-3;i<=xi+3;i++) {

			int ii=i;
			if (i<0) ii+=Npx;
			if (i>Npx-1) ii-=Npx;

			int k=(np+1)/2-1;
			double sum=0.0;
			for(int j=-k;j<=k;j++){
				int m=ii+j;
				if (m<0) m+=Npx;
				else if (m>=Npx) m-=Npx;
				sum+=oldbed[m];}

			(*bp)[ii]=sum/double(np);
		}
	}
}

void bottom::sep_sort_fsz(int num){
	/* sort the fsz�s in *fsz array
	 * 2 cases: 0 uses xsi; 4 uses xdi */

	int nfsz=(*fsz)[nf-1];
	vector<int> fsztemp = (*fsz);
	vec Srtemp = (*Sr);
	int t=0;
	while (t<nfsz){
		for(int i=0;i<Npx;i++){
			for(int j=0;j<nfsz;j++){
				if(fsztemp[j*7+num]==i){
					(*fsz)[t*7+0]=fsztemp[j*7+0];
					(*fsz)[t*7+1]=fsztemp[j*7+1];
					(*fsz)[t*7+2]=fsztemp[j*7+2];
					(*fsz)[t*7+3]=fsztemp[j*7+3];
					(*fsz)[t*7+4]=fsztemp[j*7+4];
					(*fsz)[t*7+5]=fsztemp[j*7+5];
					(*fsz)[t*7+6]=fsztemp[j*7+6];
					(*Sr)[t]=Srtemp[j];
					t+=1;
					break;}
			}
		}
	}
}

vec bottom::paramFindNeighbors(double x_p, int xi){
  /* find neighboring points of the bed */
  vec Neighbors(4,0.0);
	double a1=0.0; double a2=0.0; double b1=0.0; double b2=0.0;
	//int mm = int( floor( (x_p + 1e-10)/dx ) + 1e-10 );
	//if (tijd==3737. && xi==118-1){
	//		cerr<<"x_p="<<x_p<<"; x_p/dx="<<x_p/dx<<endl;
	//}
	int mm = int( floor( (x_p+1e-16)/dx ) );
  int m=o3(mm);
  a1 = m*dx;
  a2 = (m+1)*dx;
  b1 = (*b)[m];
  b2 = (*b)[o2(m+1)];
/*
	double dif=x_p-Npx*dx;
	if (x_p<0) x_p+=Npx*dx;
	if (x_p>=Npx*dx | x_p==L) x_p-=Npx*dx;
	//if (tijd==3782. && xi==262-1){
	//		cerr<<"x_p="<<x_p<<"; dif="<<dif<<endl;
	//}
	for(int i=xi;i<=xi+Npx;i++){
		int m = o3(i);
        a1 = m*dx;
        a2 = (m+1)*dx;
        if ( (x_p>=a1 | fabs(x_p-a1)<1e-10) && x_p<a2 && fabs(x_p-a2)>1e-10){
        //if ( x_p >= a1 && x_p < a2){
						if (tijd==3782. && xi==262-1){
							cerr<<"x_p="<<x_p<<"; a1="<<a1<<"; a2="<<a2<<endl;
						}
            b1 = (*b)[m];
            b2 = (*b)[o2(m+1)];
            break;}
	}
*/
	Neighbors[0]=a1; Neighbors[1]=a2; Neighbors[2]=b1; Neighbors[3]=b2;
	return Neighbors;
}

int bottom::findTrough(int xs, vec bed){
	int xti = 0;

  vec dhdx(Npx,0.0);
  for(int i=0;i<Npx;i++) dhdx[i]=(bed[i]-bed[o2(i-1)])/(dx);

	for(int j=xs+1;j<Npx+(xs+1);j++){

		int m = o3(j);
		if(fabs(bed[m]-bed[o2(m-1)])<1e-10 && fabs(atan(dhdx[o2(o2(m-1)-1)])*grad_2_deg+30)<2){
			 xti=m-1;
			 break;
		}
		else if(bed[m]<min(bed[o2(m-1)],bed[o2(m+1)])){
			 xti=m;
			 break;
		}
	} // for

	return o2(xti);
}

int bottom::findCrest(int xsi, vec bed){
	int xci = 0;

 	vec dhdx(Npx,0.0);
  for(int i=0;i<Npx;i++) dhdx[i]=(bed[i]-bed[o2(i-1)])/(dx);

  for(int j=xsi;j>xsi-Npx;j--){
		int m = j;
		if (j<0) m = j+Npx;
		if(bed[m]>=bed[o2(m-1)]){
			xci=m; //which is xsi
			break;
		}
	} // for

	return o2(xci);
}

double bottom::meanval(vec vr, int Np) {
	double nm = vr[0];
	for (int i = 1; i < Np; i++) {
		nm+=vr[i];
	}
	return nm/Np;
}

double bottom::minval(vec vr, int Np) {
	double nm = vr[0];
	for (int i = 1; i < Np; i++) {
		if (vr[i]<nm) nm=vr[i];
	}
	return nm;
}

vec bottom::getShape(int sepflag){
	if (sepflag==0){
		return (*b);}
	if (sepflag==1){
		return (*bp);}
    return (*b); // JW avoid warning
}

vector<int> bottom::getFsz(){
	return (*fsz);
}

vector<double> bottom::getSr(){
	return (*Sr);
}

void bottom::write_flowsep(){
	// Write flow separation zone characteristics to screen
	cerr<<"Flowsep characteristics:"<<endl;
	cerr<<setw(7)<<"xsi"<<setw(5)<<"xri"<<setw(5)<<"xci"<<setw(5)<<"xti"<<setw(5)<<"xdi"<<setw(6)<<"xcin"<<setw(6)<<"case"<<endl;
	int nfsz=(*fsz)[nf-1];
	for (int j=0;j<nfsz;j++) {
		int xsi=(*fsz)[j*7+0];
		int xri=(*fsz)[j*7+1];
		int xci=(*fsz)[j*7+2];
		int xti=(*fsz)[j*7+3];
		int xdi=(*fsz)[j*7+4];
		int xcin=(*fsz)[j*7+5];
		int fscase=(*fsz)[j*7+6];
		cerr<<j+1<<":"<<setw(5)<<xsi<<setw(5)<<xri<<setw(5)<<xci<<setw(5)<<xti<<setw(5)<<xdi<<setw(6)<<xcin<<setw(6)<<fscase<<endl;
	}
}
double bottom::detint1(vec bed){
	double bint = 0.0;
	for(int i=0;i<Npx;i++){
		bint+=dx*(bed[i]+bed[o2(i+1)])/2.;
	}
	//cerr<<"bint: "<<bint<<endl;
	return bint;
}

double bottom::detint2(vec bed){
	double bint = 0.0;
	double minbed = bed[0];
	for(int i=1;i<Npx;i++){
		if (bed[i]<minbed) minbed=bed[i];}
	for(int i=0;i<Npx;i++){
		bint+=dx*((bed[i]+bed[o2(i+1)])-minbed)/2.;
	}
	//cerr<<"bint: "<<bint<<endl;
	return bint;
}

