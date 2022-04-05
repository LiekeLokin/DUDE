// admin.cpp

#include "admin.h"
#include "Logging.h"
//#include <iostream>

namespace admin {
int Npx = 120;
}

void admin::o2_abort(int i_ex) {
	DUDE_LOG(fatal) << "o2() buffer overflow (i_in=" << i_ex << ")";
	//std::cerr << "T=" << tijd << " - ERROR: o2() buffer overflow (i_in=" << i_ex << ")" << std::endl;
	std::abort();
}
#ifndef INLINE_O2
int admin::o2(int i_ex) {
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
#endif
