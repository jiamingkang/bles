/*
 * CHak3DContIncomp8.h
 *
 *  Created on: 23 Feb 2015
 *      Author: knai20-admin
 */

#ifndef CHAK3DCONTINCOMP8_H_
#define CHAK3DCONTINCOMP8_H_

#include "inc/CHak3DElement.h"

class CHak3DContIncomp_8: public CHak3DElement {
public:
	CHak3DContIncomp_8();
	virtual ~CHak3DContIncomp_8();

public:
	double *Km;
public:

	// create stiffness matrix - also addtional matrix
	void makeKe();
	// create mass matrix
	void makeMe();
	// get element displacements
	void getDisp(double*, double*);
    void getDispC(dcompLPK *elDisp, dcompLPK *globDisp); // complex version
	// compute C8 strain tensor at a point in the element
	void strain(Coord3D, double*, double*);
	// compute C8 stress tensor at a point in the element
	void stress(Coord3D, double*, double*);
    // compute element volume
    void makeVol();
    // add self-weight to a global load vector
    void addWeight(Coord3D, double, int*, double*);
    // compute C8M sensiticies at all gauss points
    void shpSens(int, int, double**, double**, double*, double*, double*, int, double*, double, int);
    // compute C8M sensitivities at all gauss points, for eigenvalue problems
    void shpSens_Eig(int numDual, int numEig, double **disp_prim, double **disp_dual,
                     double *eig, double *wgt_case, double *gSens, int pNum, double alpha);
    // compute C8M sensitivities at all gauss points, for complex eigenvalue problems
    void shpSens_CEig(int numDual, int numEig, dcompLPK **disp_prim, dcompLPK **disp_dual,
                      dcompLPK *eig, dcompLPK *wgt_case, double *gSens, int pNum, double K_fact, double alpha);
    // compute C8M sensitivities at all gauss points, for complex eigenvalue problems
    void shpSens_CEig2(int numDual, int numEig, dcompLPK **eig_flut, double **eig_vec,
                       dcompLPK *eig, dcompLPK *wgt_case, double *gSens, int pNum, double alpha);
    // compute C8M element density sensitivity
    void elemSens(int numDual, int numCase, double **disp_prim, double **disp_dual,
                  double *wgt_fact, double *wgt_case, double *elSens, int eNum, double *Nf, int mode);
    // compute C8M element density sensitivity - complex Eigenvalue
    void elemSens_CEig(int numDual, int numEig, dcompLPK **disp_prim, dcompLPK **disp_dual,
                       dcompLPK *eig, dcompLPK *wgt_case, double *elSens, double K_fact, int eNum);
    void elemSens_CEig2(int numDual, int numEig, double **eig_vecs, dcompLPK **right, dcompLPK **left,
                             dcompLPK eig, dcompLPK *wgt_case, double *elSens, double K_fact, int eNum);
    //output
	void KeOut() {
		int i;
		for(i=0;i<300;i++) {
			std::cout << "\n" << i << " - " << Ke[i];
		}
	}
	void MeOut() {
		int i;
		for(i=0;i<300;i++) {
			std::cout << "\n" << i << " - " << Me[i];
		}
	}


};

#endif /* CHAK3DCONTINCOMP8_H_ */
