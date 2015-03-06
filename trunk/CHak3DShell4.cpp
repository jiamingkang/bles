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

// compute Shell4 element thickness sensitivity
void CHak3DShell_4::elemSens(int numDual, int numCase, double **disp_prim, double **disp_dual,
              double *wgt_fact, double *wgt_case, double *elSens, int eNum, double *Nf, int mode)
{
    int i;
    double *sens = new double [numDual](); // sum element contributions from sub elems (for each dual)

    for(i=0;i<4;i++)
    {
        sub_shells[i]->elemSens(numDual, numCase, disp_prim, disp_dual, wgt_fact, wgt_case, sens, 0, Nf, mode);
    }

    // for each dual response
    for(i=0;i<numDual;i++)
    {
        // add to overall sensitivity
        elSens[index2(eNum,i,numDual)] += 0.5*sens[i]; // half as each sub-elem only contributes half stiffness and mass
    }

    delete [] sens;
}

// compute Shell4 element thickness sensitivity w.r.t an eigenvalue
void CHak3DShell_4::EigSens(int numDual, int numEig, double **disp_prim, double **disp_dual,
             double *eig, double *wgt_case, double *elSens, int eNum, int mode)
{
    int i;
    double *sens = new double [numDual](); // sum element contributions from sub elems (for each dual)

    for(i=0;i<4;i++)
    {
        sub_shells[i]->EigSens(numDual, numEig, disp_prim, disp_dual, eig, wgt_case, sens, 0, mode);
    }

    // for each dual response
    for(i=0;i<numDual;i++)
    {
        // add to overall sensitivity
        elSens[index2(eNum,i,numDual)] += 0.5*sens[i]; // half as each sub-elem only contributes half stiffness and mass
    }

    delete [] sens;
}

// compute Shell4 flutter sensitivity - complex Eigenvalue
void CHak3DShell_4::CEigSens(int numDual, int numEig, dcompLPK **eig_flut, double **eig_vec,
              dcompLPK *eig, dcompLPK *wgt_case, double *elSens, int eNum)
{
    int i;
    double *sens = new double [numDual](); // sum element contributions from sub elems (for each dual)

    for(i=0;i<4;i++)
    {
        sub_shells[i]->CEigSens(numDual, numEig, eig_flut, eig_vec, eig, wgt_case, sens, 0);
    }

    // for each dual response
    for(i=0;i<numDual;i++)
    {
        // add to overall sensitivity
        elSens[index2(eNum,i,numDual)] += 0.5*sens[i]; // half as each sub-elem only contributes half stiffness and mass
    }

    delete [] sens;
}

// function to compute the contibution of an element to dks_du (buckling)
void CHak3DShell_4::DispSens(int numModes, double **disp_prim, double *fact, double *dk_du)
{
    int i;
    double *fact_sub = new double [numModes];
    for(i=0;i<numModes;i++){ fact_sub[i] = fact[i]*0.5; } // half from sub elems
    for(i=0;i<4;i++)
    {
        sub_shells[i]->DispSens(numModes, disp_prim, fact_sub, dk_du);
    }

    delete [] fact_sub;
}

// function to compute the contibution of an element to dvm_du (stress)
void CHak3DShell_4::DispSensVM(int numModes, double **disp_prim, double *fact, double *dvm_du)
{
    int i;
    double *fact_sub = new double [numModes];
    for(i=0;i<numModes;i++){ fact_sub[i] = fact[i]*0.25; } // average of von mises stresses from sub elems

    for(i=0;i<4;i++)
    {
        sub_shells[i]->DispSensVM(numModes, disp_prim, fact_sub, dvm_du);
    }

    delete [] fact_sub;
}

// add self-weight loading to a global loading vector
void CHak3DShell_4::addWeight(Coord3D vec, double mag, int *dof_ind, double *gload)
{
    int i;
    double mag_sub = mag*0.5; // half for sub elems
    for(i=0;i<4;i++)
    {
        sub_shells[i]->addWeight(vec, mag_sub, dof_ind, gload);
    }
}

// function to update data due to a thickness change
void CHak3DShell_4::change_thk(double tin)
{
    int i, nt[3];
    double fact = tin / t; // new / old

    t = tin; // update thickness
    volume *= fact; // update volume

    for(i=0;i<4;i++)
    {
        sub_shells[i]->change_thk(tin); // chane thk of sub elems
    }

    // update stiffness matrix - from sub-matrices
    if(Ke)
    {
        delete [] Ke;
        Ke = new double [300]();

        // assemble sub matrices into total elem matrix
        // first triangle is 0,1,2
		nt[0] = 0; nt[1] = 1; nt[2] = 2;
		addMatrix(3, 4, sub_shells[0]->getKe(), Ke, 6, nt);

		// second triangle is 0,2,3
		nt[1] = 2; nt[2] = 3;
		addMatrix(3, 4, sub_shells[1]->getKe(), Ke, 6, nt);

        // third triangle is 0,1,3
		nt[1] = 1;
		addMatrix(3, 4, sub_shells[2]->getKe(), Ke, 6, nt);

        // forth triangle is 1,2,3
		nt[0] = 1; nt[1] = 2;
		addMatrix(3, 4, sub_shells[3]->getKe(), Ke, 6, nt);

        // now divide by 2
        for(i=0;i<300;i++){Ke[i] *= 0.5;}
    }

    // update mass matrix - from sub-matrices
    if(Me)
    {
        delete [] Me;
        Me = new double [300]();

        // assemble sub matrices into total elem matrix
        // first triangle is 0,1,2
		nt[0] = 0; nt[1] = 1; nt[2] = 2;
		addMatrix(3, 4, sub_shells[0]->getMe(), Me, 6, nt);

		// second triangle is 0,2,3
		nt[1] = 2; nt[2] = 3;
		addMatrix(3, 4, sub_shells[1]->getMe(), Me, 6, nt);

        // third triangle is 0,1,3
		nt[1] = 1;
		addMatrix(3, 4, sub_shells[2]->getMe(), Me, 6, nt);

        // forth triangle is 1,2,3
		nt[0] = 1; nt[1] = 2;
		addMatrix(3, 4, sub_shells[3]->getMe(), Me, 6, nt);

        // now divide by 2
        for(i=0;i<300;i++){Me[i] *= 0.5;}
    }
}
