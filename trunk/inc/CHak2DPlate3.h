/*
 * CHak2DPlate3.h
 *
 *  Created on: 23 Feb 2015
 *      Author: knai20-admin
 */

#ifndef CHAK2DPLATE3_H_
#define CHAK2DPLATE3_H_

#include "inc/CHak2DElement.h"

// DKT plate element
// formed of one tirangle
class CHak2DPlate_3: public CHak2DElement {
public:
	CHak2DPlate_3();
	virtual ~CHak2DPlate_3();

public:
	double t;	// element thickness (restricted to const value over element)
public:

	// set material and thickness
	void setVals(double thk){t=thk;}
	// create stiffness matrix
	void makeKe();
	// create mass matrix
	void makeMe();
    // compute element volume
    void makeVol();
    // get element displacements
	void getDisp(double*, double*);
	// compute DKT strain tensor at a point in the element
	void strain(Coord, double*, double*);
	// compute DKT stress tensor at a point in the element
	void stress(Coord, double*, double*);
	// output
	void KeOut() {
		int i;
		for(i=0;i<45;i++) {
            std::cout << "\n" << i << " - " << Ke[i];
		}
	}
	void MeOut() {
		int i;
		for(i=0;i<45;i++) {
            std::cout << "\n" << i << " - " << Me[i];
		}
	}

};

#endif /* CHAK2DPLATE3_H_ */
