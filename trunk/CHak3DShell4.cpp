/*
 * CHak3DShell4.cpp
 *
 *  Created on: 23 Feb 2015
 *      Author: knai20-admin
 */

#include "CHak3DShell4.h"

CHak3DShell_4::CHak3DShell_4() {
	// TODO Auto-generated constructor stub
	// constructor sets num nodes
	numNode=4;
	nodeSpace();
    t = -1.0;
    volume = -1.0; // indicates volume not computed
	Ke=0;
	Me=0;
    Kstress=0;
    // create 4 sub elems
    int i;
    for(i=0;i<4;i++){ sub_shells[i] = new Shell3(); }
}

CHak3DShell_4::~CHak3DShell_4() {
	// TODO Auto-generated destructor stub
	// destructor (delete memory)
	delMat();
    int i;
    for(i=0;i<4;i++) { delete sub_shells[i]; }
}

