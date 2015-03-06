/*
 * CHak3DMass1.h
 *
 *  Created on: 23 Feb 2015
 *      Author: knai20-admin
 */

#ifndef CHAK3DMASS1_H_
#define CHAK3DMASS1_H_

#include "inc/CHak3DElement.h"

// Mass element
// either at single node, or connected using distribued coupling
class CHak3DMass_1: public CHak3DElement {
public:
	CHak3DMass_1();
	virtual ~CHak3DMass_1();

public:
    double mass; // element mass
    mset *connect; // set of nodes that mass is connected to
    double *Trans; // pointer to the transfer matrixts
    // data for rigid connections
    int numR; // number
    double *massR; // masses
    int *nodeR; // node numbers
    double **T_loc; // local transfer matrices for rigid element connections
public:

	// set material and thickness
	void setVals(double m, mset *ms){mass=m; connect=ms;}
    // function to add rigid mass connections
    void setRigid(int nR, double *mR, int *ndR);
	// create stiffness matrix
	void makeKe(){}; // no stiffness matrix
	// create mass matrix
	void makeMe();
    // compute element volume
    void makeVol(){volume=mass;
        int i; for(i=0;i<numR;i++){volume+=massR[i];} }; // default to total mass
    // get element displacements
	void getDisp(double*, double*);
    // return rigid attachment node nums
    int getRigid(int *ptr)
    {
        if(numR>0)
        {
            int i;
            for(i=0;i<numR;i++){ptr[i] = nodeR[i];}
        }
        return numR;
    }
    // compute transformation matrix for the distributed coupling
    void makeTrans();
    void makeTrans2();
    // return number of connected nodes
    int getNnodes(){return (connect->getPos());}
    // return pointer to connected node numbers
    int* getCNodes(){return (connect->nums);}
    // add self-weight loading to a global loading vector
    void addWeight(Coord3D vec, double mag, int *dof_ind, double *gload);
	// output
	void KeOut(){}
	void MeOut() {
		int i; int end=3*connect->getPos(); // dof
		for(i=0;i<end;i++) {
            std::cout << "\n" << i << " - " << Me[i];
		}
	}
};

#endif /* CHAK3DMASS1_H_ */
