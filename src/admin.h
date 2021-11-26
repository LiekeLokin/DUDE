// admin.h

#ifndef _ADMIN_H
#define _ADMIN_H

#include <iosfwd>

extern double tijd; // JW wordt eigenlijk alleen gebruikt voor logging (en q_in interpolatie)
extern double dz; // JW per definitie H/Npz
extern double Av;
extern double S;
extern double H;
extern double L;
extern double dx;
extern double dt;
extern std::ofstream outlog;

namespace admin{
	extern int Npx;

	/*
	 * Used in both flow.cpp and bottom.cpp
	 */
	void o2_abort(int i_ex);
#define INLINE_O2
#ifdef INLINE_O2
	inline int o2(int i_ex) {
		/*adres vertaling voor periodieke rvw ten behoeve van dz*/
		int uit = 0;
		if (i_ex >= 0 && i_ex < Npx)
			uit = i_ex;
		else if (i_ex == -1)
			uit = Npx-1;
		else if (i_ex == Npx)
			uit = 0;
		else
			o2_abort(i_ex);
		return uit;
	}
#else
	int o2(int i_ex);
#endif
}

#endif
