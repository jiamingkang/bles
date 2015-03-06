/*
 * CHak2DPlate4.h
 *
 *  Created on: 23 Feb 2015
 *      Author: knai20-admin
 */

#ifndef CHAK2DPLATE4_H_
#define CHAK2DPLATE4_H_

#include "inc/CHak2DElement.h"

class CHak2DPlate_4: public CHak2DElement {
public:
	CHak2DPlate_4();
	virtual ~CHak2DPlate_4();

public:
	double t;	// element thickness (restricted to const value over element)
	int diag;	// store which diagonal the two triangles are joined 1 = node 1-3, 2 = node 2-4
public:
	// set material and thickness
	void setVals(double thk, int d){t=thk;diag=d;}
	// get diag
	int getDiag() {return(diag);}
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
    // fucntion to work out which triangle a point lies in
    int findTri(Coord);
	// output
	void KeOut() {
		int i;
		for(i=0;i<78;i++) {
            std::cout << "\n" << i << " - " << Ke[i];
		}
	}
	void MeOut() {
		int i;
		for(i=0;i<78;i++) {
            std::cout << "\n" << i << " - " << Me[i];
		}
	}

};

#endif /* CHAK2DPLATE4_H_ */
