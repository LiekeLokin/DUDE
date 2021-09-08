// flow.h

#ifndef _FLOW_H
#define _FLOW_H
#include "linalg.h"
#include "admin.h"
#include <fstream>
#include <cstdlib>
using namespace std;
using namespace admin;

class flow{
	private:
		spMat *A;
		crLU *prLU;
		int prec;
		vec *b;
		vec *iu;
		vec *beta;
		vec *alpha;
		vec *u0;
		vec *Avx;
		vec *Sx;
		/*spMat *Am;*/
		void zl(int r,int i);
		void zr(int r,int i);
		void wl(int r,int j_ex,int i);
		void wr(int r,int j_ex,int i);
		void w_from_uA(int r,int j_ex,int i_ex,double co);
		double w_from_u(int j_ex,int i_ex);
		void a0l(int r,int j,int i);
		void a0r(int r,int j,int i);
		void a1l(int r,int j,int i);
		void a1r(int r,int j,int i);
		void a1l_opperv(int r,int i);
		void a1r_opperv(int r,int i);
		void a1l_bodem(int r,int i);
		void a1r_bodem(int r,int i);
		void gzl(int r,int i);
		void gzr(int r,int i);
		void vl(int r,int j,int i);
		void vr(int r,int j,int i);
		void vl_bodem(int r,int i);
		void vr_bodem(int r,int i);
		void vl_opperv(int r,int i);
		void vr_opperv(int r,int i);
		void last_zetal();
		void last_zetar();
		void vulb();
		void vulA();
		void vulu();
		void initIu();
		void dzs_init(vec bottom_state);
		void det_AvS(vec bottom_state);
	public:
		flow();
		~flow();
		int solve(vec bottom_state);
		int solve_gm(vec bottom_state,int gmn);
		void reprec();
		vec getiu();
		double check_qsp();
		void u_b(vec &u0);
		void resetIu();
		void resetIu(vec u);
		void u_anal(double eps);
		/*void testNewton(vec bottom_state,double eps);*/
		void write_velocities(double tijd, vec bottom_state, vec u0_b);
		//void write_zeta(double tijd);
		vec getZeta();
		double zetaint1();
		double zetaint2();
};




#endif
