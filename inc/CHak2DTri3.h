/*
 * CHak2DTri3.h
 *
 *  Created on: 23 Feb 2015
 *      Author: knai20-admin
 */

#ifndef CHAK2DTRI3_H_
#define CHAK2DTRI3_H_

#include "inc/CHak2DElement.h"

class CHak2DTri_3: public CHak2DElement {
public:
	CHak2DTri_3();
	virtual ~CHak2DTri_3();
	double t;	// element thickness (restricted to const value over element)

public:
	// set material and thickness
	void setVals(double thk){t=thk;}
	// create stiffness matrix
	void makeKe();
	// create mass matrix
	void makeMe();
    void makeMe2();
    // compute element volume
    void makeVol();
    // get element displacements
	void getDisp(double*, double*);
	// compute DKT strain tensor at a point in the element
	void strain(Coord, double*, double*);
	// compute DKT stress tensor at a point in the element
	void stress(Coord, double*, double*);
    // compute Tri3 element thickness sensitivity
    void elemSens(int numDual, int numCase, double **disp_prim, double **disp_dual, double **fg,
                        double *wgt_fact, double *wgt_case, double *elSens, int eNum, double *Nf);
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

#endif /* CHAK2DTRI3_H_ */
