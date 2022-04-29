// flow.h

#ifndef _FLOW_H
#define _FLOW_H

#include "FlowConfig.h"
#include "linalg.h"

class flow {
private:
	const FlowConfig cfg;
	const int Npx;
	const int Npz;
	const int nt;
	spMat A;
	crLU *prLU;
	int prec;
	vec b;
	vec iu;
	vec beta;
	vec alpha;
	vec u0;
	/*spMat *Am;*/
	int o(int j_ex,int i_ex,int v) const;
	void zl(int r,int i);
	void zr(int r,int i);
	//void wl(int r,int j_ex,int i);
	//void wr(int r,int j_ex,int i);
	void w_from_uA(int r,int j_ex,int i_ex,double co);
	double w_from_u(int j_ex,int i_ex) const;
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
	void dzs_init(const vec& bottom_state);
public:
	flow() = delete;
	flow(const FlowConfig& cfg);
	~flow();
	int solve(const vec& bottom_state);
	int solve_gm(const vec& bottom_state,int gmn);
	void reprec();
	vec getiu() const;
	double check_qsp() const;
	void u_b(vec &u0) const;
	void resetIu();
	void resetIu(const vec& u);
	//void u_anal(double eps);
	/*void testNewton(vec bottom_state,double eps);*/
	void write_velocities(double tijd, const vec& bottom_state, const vec& u0_b) const;
	//void write_zeta(double tijd);
	vec getZeta() const;
	double zetaint1() const;
	double zetaint2() const;
};

#endif
