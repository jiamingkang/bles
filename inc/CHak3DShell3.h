/*
 * CHak3DShell3.h
 *
 *  Created on: 23 Feb 2015
 *      Author: knai20-admin
 */

#ifndef CHAK3DSHELL3_H_
#define CHAK3DSHELL3_H_

#include "inc/CHak3DElement.h"

// 3 node 3D shell element (6 dof / node)
// DKT + Tri3 elements
class CHak3DShell_3: public CHak3DElement {
public:
	CHak3DShell_3();
	virtual ~CHak3DShell_3();

    double t; // thickness
    double Kbend[171]; // bending stiffness matrix
    double Kmem[171];  // membrane stiffness matrix
    double *Kstress;   // stress stiffness matrix for buckling problems
    double Trans[9];   // coordinate transformation matrix
    Coord3D local[3];  // local 2D coordinate system (z coords only for pointer referencing)
    DKT *plate;      // associated plate element
    Tri3 *mem;       // associated membrane element
    bool tr_yes;
public:

    // delete stress stiffness matrix
    void delKs(){if(Kstress){delete [] Kstress; Kstress=0;}}
    // return stess stiffness matrix
    double* getKs(){return Kstress;}
    // set thickness value
    void setVals(double thk){t=thk; plate->setVals(thk); mem->setVals(thk);}
    // compute transformation matrix and local coords
    void makeTrans();
	// create stiffness matrix
	void makeKe();
	// create mass matrix
	void makeMe();
    // create stress stiffness matrix
    void makeKs(double*);
	// get element displacements
	void getDisp(double*, double*);
    // compute Shell3 element thickness sensitivity
    void elemSens(int numDual, int numCase, double **disp_prim, double **disp_dual,
                          double *wgt_fact, double *wgt_case, double *elSens, int eNum, double *Nf, int mode);
    // compute Shell3 element thickness sensitivity w.r.t an eigenvalue
    void EigSens(int numDual, int numEig, double **disp_prim, double **disp_dual,
                 double *eig, double *wgt_case, double *elSens, int eNum, int mode);
    // compute Shell3 flutter sensitivity - complex Eigenvalue
    void CEigSens(int numDual, int numEig, dcompLPK **eig_flut, double **eig_vec,
                  dcompLPK *eig, dcompLPK *wgt_case, double *elSens, int eNum);
    // function to compute the contibution of an element to dks_du (buckling)
    void DispSens(int numModes, double **disp_prim, double *fact, double *dk_du);
    // function to compute the contibution of an element to dvm_du (stress)
    void DispSensVM(int numModes, double **disp_prim, double *fact, double *dvm_du);
	// compute C8 strain tensor at a point in the element
	void strain(Coord3D, double*, double*);
	// compute C8 stress tensor at a point in the element
	void stress(Coord3D, double*, double*);
    // compute von Mises stress (at an input point)
    double vonMises(Coord3D, double*);
    // compute element volume
    void makeVol();
    // function to update data due to a thickness change
    void change_thk(double tin);
    // return thickness
    double getThk(){return t;}
    // add self-weight loading to a global loading vector
    void addWeight(Coord3D vec, double mag, int *dof_ind, double *gload);
	// output
	void KeOut() {
		int i;
		for(i=0;i<171;i++) {
			std::cout << "\n" << i << " - " << Ke[i];
		}
	}
	void MeOut() {
		int i;
		for(i=0;i<171;i++) {
			std::cout << "\n" << i << " - " << Me[i];
		}
	}
};

#endif /* CHAK3DSHELL3_H_ */
