/*
	CHakSensitivity.cpp

	Created on: Nov 24, 2014
	Author: Peter Dunning, JeeHang Lee

 	-- This file is part of Topology Optimisation Opensource Project,
 	owned by MSO (Multidisciplinary and Structural Optimisation) Research Group
 	(http://people.bath.ac.uk/ens/MSORG/index.html) at University of Bath.

 	The project is led by Dr. Hyunsun Alicia Kim
 	(http://www.bath.ac.uk/mech-eng/people/kim/).

	-- This is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    -- This is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this. If not, see <http://www.gnu.org/licenses/>.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "CommonTypes.h"

#include "CHakFiniteElement.h"	// jeehang.lee@gmail.com, solving required.
#include "CHakSolver.h"
#include "CHakSensitivity.h"

//
// Constructor & Destructor
//
CHakSensitivity::CHakSensitivity() {
	// TODO Auto-generated constructor stub

}

CHakSensitivity::~CHakSensitivity() {
	// TODO Auto-generated destructor stub
}

//
// Implementation - Interfaces
//

// calculate sensitivies using least squares of integration points for AFG method
void CHakSensitivity::AFG_Sens(CHakMesh *inMesh, CHakBoundary *bound_in, double *alpha, CHakMaterial *inMat,  double *Nsens,
				double **prim, double **dual, int numDual, int numCase, double *wgt, Coord *gCoord,
					double aMin, int mode, double *fact, bool sw, Coord *acc)
{
	// read in mesh data
	double h = inMesh->m_lenEdge;
	double thk = inMesh->m_thinkness;
	int NumNodes = inMesh->m_numNodes;
	int elemX = inMesh->m_elemX;
	int elemY = inMesh->m_elemY;
	int NumElem = inMesh->m_numElem;
	Elem **Number = inMesh->m_pNumber;
	Coord *NodeCoord = inMesh->m_pNodeCoord;
	double rad = 2.0 * h; // hard coded for two elements around a node

	int irad = 2; //ceil(rad/h); // search space
	rad *=rad; // input squared dist to Lsens
	Coord *AuxNodes = bound_in->m_pAuxNodes;
	int numBound = bound_in->m_numBound;
	Bseg *Boundary = bound_in->Bound;
	int Ntot = NumNodes + bound_in->m_numAux;
	double alMin = (aMin < 1.0e-6) ? 1.0e-6 : aMin; // numerical tollerance for sens calc
	if(mode==1){ alMin = (aMin < 1.0e-2) ? 1.0e-2 : aMin; } // increase for eigenvalue sensitivities
	else if(mode==2){alMin = 0.0;} // for compliant mechanisms

	int i,n,m,o,num,ind,node,ex,ey,xMin,xMax,yMin,yMax;
	double atemp;
	//if(mode==2){w1 = -1.0/wgt[0]; w2 = wgt[1]*w1*w1;} // extra weights for complinat mech disp constraint
	int Gcount = 0;		// varible to track number of points evaluated
	int tnodes[4];
	int gpoints = NumElem * 4; // number of sensitivity values (4 per element)
	double *gSens = (double *) calloc(gpoints*numDual,sizeof(double));	// define memory for gauss point sensitivities
	CHakMaterial *mptr; // pointer for material

	int mat_count = 0; // count for designable material
	double mtemp; // variables for designable material modification
	CHakMaterial fgm;
	CHakMultiMaterial multiMaterial;


	// Step 1. Compute gauss point senstivities

	// For All Elements
	for(m=0;m<elemY;m++)
	{
		for(n=0;n<elemX;n++)
		{
			num = Number[n][m].n; // Element number
			atemp = alpha[num];
			if(atemp > alMin) // If Element isn't OUT
			{
				tnodes[0] = Number[n][m].a;
				tnodes[1] = Number[n][m].b;
				tnodes[2] = Number[n][m].c;
				tnodes[3] = Number[n][m].d;

				o = inMesh->m_pMaterialType[num]; // material number
				mptr = &inMat[o]; // point to correct material

                // modify matrices for designable materials
                if(inMesh->m_bDesignableMaterial && num==inMesh->m_pMaterialElems[mat_count])
                {
                    // material mix variable
                    mtemp = inMesh->m_pMaterialVars[mat_count++];

                    // modulus
                    if(inMesh->m_bMaterialLin)
                    {
                    	fgm.m_e =  (1.0-mtemp)*inMat[inMesh->m_materialOne].m_e + mtemp*inMat[inMesh->m_materialTwo].m_e;
                    }
                    else
                    {
                    	fgm.m_e = multiMaterial.HS_mat(mtemp, 0.5, &inMat[inMesh->m_materialOne], &inMat[inMesh->m_materialTwo]);
                    }

                    // density
                    fgm.m_rho =  mtemp*inMat[inMesh->m_materialTwo].m_rho + (1.0-mtemp)*inMat[inMesh->m_materialOne].m_rho;

                    // material prop matrix
                    mtemp = fgm.m_e / inMat[inMesh->m_materialOne].m_e;
                    for(i=0;i<9;i++) {
                        //fgm.mat[i] = mtemp*inMat[inMesh->m_materialOne].mat[i] + (1.0-mtemp)*inMat[inMesh->m_materialTwo].mat[i];
                        fgm.m_mat[i] = mtemp*inMat[inMesh->m_materialOne].m_mat[i];
                    }

                    fgm.m_v = inMat[inMesh->m_materialOne].m_v;
                    mptr = &fgm;
                }
                else
                {
                    mptr = &inMat[o]; // point to correct material (not FGM)
                }

				//calcualte sensitivity at Gauss points
				if(mode==1) // eigenvalue sensitivites
				{
					GaEigSens_Q4(tnodes, prim, dual, atemp, h, thk, mptr, Gcount, gSens, numCase, wgt);
				}
				else if(mode==2) // compliant mechanism sensitivites
				{
					double *gSens2 = (double *) calloc(8, sizeof(double));
					double *gSens3 = (double *) calloc(4, sizeof(double));	// temp memory for gauss point sensitivities

					GaSens_Q4(tnodes, prim, dual, atemp, h, thk, mptr, 0, gSens2, 1, 2, wgt, sw, acc);
					GaSens_Q4(tnodes, &dual[1], &dual[1], atemp, h, thk, mptr, 0, gSens3, 1, 1, wgt, sw, acc);

					// process sensitivities for compliant mech
					for(i=0;i<4;i++)
					{
						ind = (Gcount+i)*3; // point to place in gSens
						gSens[ind] = -fact[0]*gSens2[(2*i)+1] - fact[1]*gSens3[i];
						gSens[ind+1] = gSens2[2*i]*fact[4] + fact[2]*gSens2[(2*i)+1] + fact[3]*gSens3[i];
						gSens[ind+2] = gSens3[i]*fact[5];
					}

					free(gSens2);
					free(gSens3);
				}
				else
				{ GaSens_Q4(tnodes, prim, dual, atemp, h, thk, mptr, Gcount, gSens, numCase, numDual, wgt, sw, acc); }
			}
            else if(inMesh->m_bDesignableMaterial && num==inMesh->m_pMaterialElems[mat_count]){mat_count++;} // keep counting
			Gcount += 4; // update point count
		}
	}

	// Step 2. Compute node sensitivities using least squares method

	double *sens_temp = (double *) malloc(numDual * sizeof(double)); // array for sens of each dual state

	// calculate smoothed sensitivities for all boundary (or non-OUT) nodes
	bool *done = (bool *) calloc(Ntot, sizeof(bool));
	double ftemp;

	// first search through each boundary segment
	for(n=0;n<numBound;n++)
	{
		ey = Boundary[n].e / elemX;
		ex = Boundary[n].e - elemX*ey; // element indices
		xMin = ex-irad; xMax=ex+irad+1;
		yMin = ey-irad; yMax=ey+irad+1; // search limits
		xMin = (xMin < 0) ? 0 : xMin;
		yMin = (yMin < 0) ? 0 : yMin;
		xMax = (xMax > elemX) ? elemX : xMax;
		yMax = (yMax > elemY) ? elemY : yMax; // adjust search limits

		node=Boundary[n].n1; // 1st node number
		if(!done[node]) // if sens not already computed
		{
			ftemp = 1.0;
			if(node < NumNodes) // If node is a grid node
			{
				num = Lsens(&NodeCoord[node], xMax, xMin, yMax, yMin, alMin, alpha, rad, gCoord, gSens, Number, 4, numDual, sens_temp);
			}
			else // Otherwise node is an auxillary node
			{
				m = node - NumNodes; // position in AuxNodes array
				num = Lsens(&AuxNodes[m], xMax, xMin, yMax, yMin, alMin, alpha, rad, gCoord, gSens, Number, 4, numDual, sens_temp);
			}
			for(i=0;i<numDual;i++)
			{
				m = Ntot*i + node; // point to correct place in Nsens
				Nsens[m] = ftemp * sens_temp[i]; // multiply smoothed sensitivities by weight
			}
			done[node] = true;
		}
		node=Boundary[n].n2; // 2nd node number
		if(!done[node]) // if sens not already computed
		{
			ftemp = 1.0;
			if(node < NumNodes) // If node is a grid node
			{
				num = Lsens(&NodeCoord[node], xMax, xMin, yMax, yMin, alMin, alpha, rad, gCoord, gSens, Number, 4, numDual, sens_temp);
			}
			else // Otherwise node is an auxillary node
			{
				m = node - NumNodes; // position in AuxNodes array
				num = Lsens(&AuxNodes[m], xMax, xMin, yMax, yMin, alMin, alpha, rad, gCoord, gSens, Number, 4, numDual, sens_temp);
			}
			for(i=0;i<numDual;i++)
			{
				m = Ntot*i + node; // point to correct place in Nsens
				Nsens[m] = ftemp * sens_temp[i]; // multiply smoothed sensitivities by weight
			}
			done[node] = true;
		}
	}

	// clear memory
	free(done);
	free(gSens);
	free(sens_temp);
}

// Function that calculates the sensitivity of a node by a least squares (2nd order) filter of near-by gauss points
int CHakSensitivity::Lsens(Coord *pt, int xMax, int xMin, int yMax, int yMin, double aMin, double *alpha, double r2,
					Coord *gCoord, double *gSens, Elem **Number, int wFlag, int numDual, double *out)
{
	int i,j,n,m,num,ind;  // Incrementor
	double atemp,ftemp,x2,y2,xb,yb;  // variables for squared co-ordinates
	double wgt;	 // varible for weight function
	int pts[150]; // array of global gPoint numbers
	double dist[150]; // array of distances
	double aPts[150]; // array of area ratios (of the points)

	int count = 0;	 // initialize count of number of points used to zero

	// only look at nearby elements
	for(m=yMin;m<yMax;m++)
	{
		for(n=xMin;n<xMax;n++)
		{
			num = Number[n][m].n;	// Element number
			atemp = alpha[num];		// area ratio
			// check element is not OUT
			if(atemp > aMin)
			{
				ind = 4*num; // point to correct place in gCoord
				for(i=0;i<4;i++) // for all gauss points
				{
					ftemp = pt->x - gCoord[ind].x;
					x2 = ftemp * ftemp;	 // squared x-distance from current gauss point to node of interest
					ftemp = pt->y - gCoord[ind].y;
					y2 = ftemp * ftemp;	 // squared y-distance from current gauss point to node of interest

					ftemp = x2 + y2;  // squared distance between gauss point and node

					if(ftemp < r2) // if distance less than the radius squared, then add data to arrays
					{
						dist[count] = sqrt(ftemp); // save distance
						aPts[count] = atemp; // save area ratio
						pts[count++] = ind;  // save gPoint number
					}
					ind++; // next gPoint
				}
			}
		}
	}

	if(count < 10)  // if not enough points within radius then point must be on an island!
	{
		printf("\nWarning island node found: ");
		for(j=0;j<numDual;j++)
		{
			out[j] = 0.0;  // return zero
		}
		return(1);
	}

	double *A = (double *) malloc( 6 * count * sizeof(double)); // array to store gpoint values
	double *B = (double *) malloc(numDual * count * sizeof(double)); // array to store rhs' of equation
	double *Bmax =  (double *) calloc(numDual, sizeof(double)); // array to store maximum sens values
	double *Bmin =  (double *) calloc(numDual, sizeof(double)); // array to store minimum sens values

	for(i=0;i<count;i++)
	{
		ind = pts[i]; // gPoint number

		// calculate the weight function
		switch (wFlag)
		{
			case 1:
					wgt = 1.0; // no weighting
					break;
			case 2:
					wgt = 1.0 / sqrt(dist[i]); // weighted by inverse distance
					break;
			case 3:
					wgt = sqrt(aPts[i]); // weighted by area
					break;
			case 4:
					wgt = sqrt(aPts[i] / dist[i]); // weighted by area & inverse distance
					break;
			default:
					printf("\nERROR!! wFlag out of range!");
		}

		xb = gCoord[ind].x - pt->x; // relative x-coord
		yb = gCoord[ind].y - pt->y; // relative y-coord

		A[i] = wgt;
		A[count + i] = xb * wgt;
		A[(2 * count) + i] = yb * wgt;
		A[(3 * count) + i] = xb * yb * wgt;
		A[(4 * count) + i] = xb * xb * wgt;
		A[(5 * count) + i] = yb * yb * wgt;

		m = numDual*pts[i]; // point to correct place in gSens
		for(j=0;j<numDual;j++)
		{
			n = j * count;
			num = m + j;
			B[n + i] = gSens[num] * wgt; // save senstivity for each dual case
			if(gSens[num] > Bmax[j]){ Bmax[j] = gSens[num]; }
			else if(gSens[num] < Bmin[j]){ Bmin[j] = gSens[num]; }
		}
	}

	// solve least squares problem using LAPACK routine
	// jeehang.lee@gmail.com -- located in solve.c, resolving required.
	CHakSolver cSolver;
	cSolver.d_lsLPK(6, count, numDual, A, B);

	// Finally evaluate sensitivity at the node using the co-efficients
	for(j=0;j<numDual;j++)
	{
		num = count*j;
		if(B[num] > Bmax[j]){ B[num] = Bmax[j]; }
		else if(B[num] < Bmin[j]){ B[num] = Bmin[j]; }
		out[j] = B[num];
		//if(out[j] > 1.0e3)
		//{ printf("sens %i = %12.4e",j+1,out[j]); }
	}

	free(A);
	free(B);
	free(Bmax);
	free(Bmin);
	return(0);
}

// Function to calculate sensitivity values for an element at 4 gauss points
// for for a plane 4-node element (Q4)
void CHakSensitivity::GaSens_Q4(int *tnodes, double **prim, double **dual, double alpha, double h, double t,
				CHakMaterial *inMat, int Gcount, double *gSens, int numCase, int numDual, double *wgt, bool sw, Coord *acc)
{
	int i,j,k,n,p,temp,temp2,ind,ind2;	// incrementors etc
	int totDof = 4*NUM_DOF; // total dof per element
	double *Eprim = (double *) malloc(totDof*sizeof(double));	// Primary Element displacement array
	double *Edual = (double *) malloc(totDof*numDual*sizeof(double));	// Adjoint Element displacement array
	double stress[3], strain[3], sens;	 // Variabes for strain tensors
    Coord ug;
    Coord *pg = (Coord *) malloc(numDual*sizeof(Coord)); // primal & dual displacements at gauss point

	// copy material property matrix
	double Emat[9];
	for(i=0;i<9;i++){Emat[i] = inMat->m_mat[i];}

	double h_fact = (t*alpha) / (4.0*h*h); // factor to re-dimension senstivity (inc area ratio & thickness)
    double sw_fact;
    //if(sw){ sw_fact = alpha * inMat->rho * t; } // const factor for self-weight
    if(sw){ sw_fact = inMat->m_rho * t; } // const factor for self-weight

	// Evaluate sensitivities at each gauss point
	for(p=0;p<4;p++)
	{
		// for each load case (should be 1, unless multi-load compliance is obj or constr)
		for(n=0;n<numCase;n++)
		{
			// get element displacement arrays
			for(i=0;i<4;i++)
			{
				temp = i * NUM_DOF; // local
				temp2 = tnodes[i] * NUM_DOF; // global

				for(j=0;j<NUM_DOF;j++)
				{
					ind=temp+j; // local
					ind2=temp2+j; // global

					Eprim[ind] = prim[n][ind2]; // primary disp

					for(k=0;k<numDual;k++)
					{
						// NB: **dual is really a matix of pointers
						// row num = load case num (n)
						// col num = dual state num (k)
						if(dual[n*numDual + k]) // if dual state exists
						{
							Edual[totDof*k + ind] = dual[n*numDual + k][ind2]; // dual disp
						}
					}
				}
			}

            // get displacements at the gauss point (if self-weight)
            if(sw)
            {
                // re-set to zero
                ug.x=0.0; ug.y=0.0;
                for(i=0;i<numDual;i++){ pg[i].x=0.0; pg[i].y=0.0; }

                for(i=0;i<4;i++)
                {
                    temp = i * NUM_DOF; // local
                    ug.x += Q4_inter[p][i]*Eprim[temp]; // sum Ni * Uxi
                    ug.y += Q4_inter[p][i]*Eprim[temp+1]; // sum Ni * Uyi

                    // dual displacements
                    for(k=0;k<numDual;k++)
					{
						if(dual[n*numDual + k]) // if dual state exists
						{
                            ind = totDof*k + temp;
							pg[k].x += Q4_inter[p][i]*Edual[ind];
                            pg[k].y += Q4_inter[p][i]*Edual[ind + 1];
						}
					}
                }
            }

			ind = (Gcount+p)*numDual; // point to place in gSens

			// compute stresses using Eprim ( stress = [E][B]{u} )
			cblas_dgemv(CblasRowMajor, CblasNoTrans, 3, totDof, 1.0, BMat[p], totDof, Eprim, 1, 0.0, strain, 1);
			cblas_dgemv(CblasRowMajor, CblasNoTrans, 3, 3, 1.0, Emat, 3, strain, 1, 0.0, stress, 1);

			// for each dual state
			for(j=0;j<numDual;j++)
			{
				if(dual[n*numDual + j]) // if dual state exists
				{
					// compute strains using Edual ( strain = [B]{p} )
					cblas_dgemv(CblasRowMajor, CblasNoTrans, 3, totDof, 1.0, BMat[p], totDof, &Edual[j*totDof], 1, 0.0, strain, 1);

					// compute sensitivity ( stress(u) * strain(p) )
					sens = stress[0]*strain[0] + stress[1]*strain[1] + stress[2]*strain[2];
					sens *= h_fact;

					// add senstivity to main data array
					gSens[ind + j] += sens * wgt[n];	// store sensitivity value for later

                    // add sensitivity for self-weight (if required)
                    if(sw)
                    {
                        sens = sw_fact * ( (ug.x + pg[j].x)*acc[n].x + (ug.y + pg[j].y)*acc[n].y );
                        gSens[ind + j] -= sens * wgt[n];
                    }
				}
			}
		}
	}

	// clear memory
	free(Eprim);
	free(Edual);
    free(pg);
}

// Function to calculate additional sensitivity part for eigenvalues
// for for a plane 4-node element (Q4)
void CHakSensitivity::GaEigSens_Q4(int *tnodes, double **prim, double **dual, double alpha, double h, double t,
		CHakMaterial *inMat, int Gcount, double *gSens, int num_eig, double *eig)
{
	int i,j,k,p,temp,temp2,ind,ind2;	// incrementors etc
	int totDof = 4*NUM_DOF; // total dof per element
	double *Eprim = (double *) malloc(totDof*sizeof(double));	// Primary Element displacement array
	double *Edual = (double *) malloc(totDof*sizeof(double));	// Adjoint Element displacement array
	double stress[3], strain[3], sens, ftemp;	 // Variabes for strain tensors
	double Iprim, Idual; // variables for interpolated displacements

	// copy material property matrix
	double Emat[9];
	for(i=0;i<9;i++){Emat[i] = inMat->m_mat[i];}

	double h_fact = (t*alpha) / (4.0*h*h); // factor to re-dimension senstivity (inc area ratio & thickness)

	// get element displacement arrays
	for(k=0;k<num_eig;k++) // for each eigenvector
	{
		for(i=0;i<4;i++) // each node
		{
			temp = i * NUM_DOF; // local
			temp2 = tnodes[i] * NUM_DOF; // global

			for(j=0;j<NUM_DOF;j++) // each dof
			{
				ind=temp+j; // local
				ind2=temp2+j; // global

				Eprim[ind] = prim[k][ind2]; // primary disp
				Edual[ind] = dual[k][ind2]; // dual disp
			}
		}

		// Evaluate sensitivities at each gauss point
		for(p=0;p<4;p++)
		{
			// compute stresses using Eprim ( stress = [E][B]{u} )
			cblas_dgemv(CblasRowMajor, CblasNoTrans, 3, totDof, 1.0, BMat[p], totDof, Eprim, 1, 0.0, strain, 1);
			cblas_dgemv(CblasRowMajor, CblasNoTrans, 3, 3, 1.0, Emat, 3, strain, 1, 0.0, stress, 1);

			// compute strains using Edual ( strain = [B]{p} )
			cblas_dgemv(CblasRowMajor, CblasNoTrans, 3, totDof, 1.0, BMat[p], totDof, Edual, 1, 0.0, strain, 1);

			// compute sensitivity ( stress(u) * strain(p) )
			sens = stress[0]*strain[0] + stress[1]*strain[1] + stress[2]*strain[2];
			sens *= h_fact;

			// add mass senstivity ( lambda * rho * {u}.{p} )
			// interpolate displacements at the gauss point and sum squares
			ftemp = 0.0;
			for(j=0;j<NUM_DOF;j++)
			{
				Iprim = 0.0; Idual = 0.0;
				for(i=0;i<4;i++) // position in Edisp
				{
					// interpolate
					ind = (i*NUM_DOF)+j;
					Iprim += Q4_inter[p][i] * Eprim[ind];
					Idual += Q4_inter[p][i] * Edual[ind];
				}
				ftemp += Iprim*Idual; // add u.p for each direction
			}

			sens -= eig[k]*ftemp*inMat->m_rho*alpha*t; // extra factor

			// add senstivity to main data array
			ind = (Gcount+p)*num_eig + k; // point to place in gSens
			gSens[ind] = sens;	// store sensitivity value for later (minus for max eig value??)
		}
	}

	// clear memory
	free(Eprim);
	free(Edual);
}

// function that computes bar senstivites for compliance (possible multi-load case)
void CHakSensitivity::barSens(CHakMesh *inMesh, double *bar_sens, double **prim, int numCase, double *wgt)
{
	// read data
	int numBars = inMesh->m_numBars;
	double fact = inMesh->m_pBarMaterialType->m_e / inMesh->m_lenEdge;
	Bseg *bar = inMesh->m_pBarNums;

	int b,n,d1,d2;
	double udiff;

	// for all bars
	for(b=0;b<numBars;b++)
	{
		// get the dof
		d1 = bar[b].n1; d2 = bar[b].n2;

		// for all cases
		for(n=0;n<numCase;n++)
		{
			// compute uT (E/L) u for this load case
			udiff = prim[n][d1] - prim[n][d2];
			udiff *= udiff*fact; // senstivity

			bar_sens[b] -= wgt[n] * udiff; // add sensitivity to overall
		}
	}
}

// function to compute designable bc sensitvities for compliance (possible multi-load case)
void CHakSensitivity::bcSens(CHakMesh *inMesh, double *bc_sens, double **prim, int numCase, double *wgt)
{
	int n,m,o,c,i,j,temp,ind;
	int elemX = inMesh->m_elemX;
	int numBC = inMesh->m_numBc;
	Elem **Number = inMesh->m_pNumber;
	int tnodes[4];
	double ftemp,disp;

	// for all
	for(o=0;o<numBC;o++)
	{
		// compute element indices
		m = floor((double)inMesh->m_pBcNums[o] / (double)elemX);
		n = inMesh->m_pBcNums[o] - m*elemX;

		// Read in grid node numbers of current element
		tnodes[0] = Number[n][m].a;
		tnodes[1] = Number[n][m].b;
		tnodes[2] = Number[n][m].c;
		tnodes[3] = Number[n][m].d;

		for(c=0;c<numCase;c++) // each load case
		{
			// compute displacement vector dot product
			disp = 0.0;
			for(i=0;i<4;i++) // each node
			{
				temp = tnodes[i] * NUM_DOF; // node location in global dof

				for(j=0;j<NUM_DOF;j++) // each dof
				{
					ind=temp+j; // global dof
					ftemp = prim[c][ind]; // primary disp
					disp += ftemp*ftemp; // dot product
				}
			}

			// sensitivity (including penalization power)
			ftemp = inMesh->m_K_bc[o]; // current design variable
			bc_sens[o] -= wgt[c] * disp * 3.0*ftemp*ftemp * inMesh->m_K0_bc;
		}
	}
}

// function to compute designable material design varibles for compliance
void CHakSensitivity::matSens_comp(CHakMesh *inMesh, CHakMaterial *inMat, double *KE, double *mat_sens, int numCase, double *wgt,
                    double **prim, double **dual, double *alpha, double aMin, bool sw, Coord *acc)
{
    int i,j,n,m,ind,ind2,temp,temp2;
	int tnodes[4];
	double uKu, sx, sy;
	int m_elemX = inMesh->m_elemX;
	int numVar = inMesh->m_numDesignableMaterialt;
	Elem **Number = inMesh->m_pNumber;

	int totDof = 4*NUM_DOF; // total dof per element
	double *Eprim = (double *) malloc(totDof*sizeof(double));	// primary element disp array
    double *Edual = (double *) malloc(totDof*sizeof(double));	// primary element disp array
	double *aux = (double *) malloc(totDof*sizeof(double));

	// realtive material properties
	double E1 = inMat[inMesh->m_materialOne].m_e;
	double E2 = inMat[inMesh->m_materialTwo].m_e;
    double rho = inMat[inMesh->m_materialTwo].m_rho - inMat[inMesh->m_materialOne].m_rho;
    rho *= 0.25 * inMesh->m_thinkness * inMesh->m_lenEdge * inMesh->m_lenEdge; // multiply by volume to get mass / node

	// sensitivity factor
	E1 = (E2 - E1)/E1;

    // for each design variable
	for(i=0;i<numVar;i++)
	{
        if(alpha[inMesh->m_pMaterialElems[i]] > aMin)
        {
            // compute element indices
            m = floor((double)inMesh->m_pMaterialElems[i] / (double)m_elemX);
            n = inMesh->m_pMaterialElems[i] - m*m_elemX;

            // Read in grid node numbers of current element
            tnodes[0] = Number[n][m].a;
            tnodes[1] = Number[n][m].b;
            tnodes[2] = Number[n][m].c;
            tnodes[3] = Number[n][m].d;

            // for each loadcase
            for(j=0;j<numCase;j++)
            {
                // get prim and dual for current element
                for(n=0;n<4;n++) // each node
                {
                    temp = n * NUM_DOF; // local
                    temp2 = tnodes[n] * NUM_DOF; // global

                    for(m=0;m<NUM_DOF;m++) // each dof
                    {
                        ind=temp+m; // local
                        ind2=temp2+m; // global
                        Eprim[ind] = prim[j][ind2];
                        Edual[ind] = dual[j][ind2];
                    }
                }

                // multiply dk = p^T K u
                cblas_dgemv(CblasRowMajor, CblasNoTrans, totDof, totDof, 1.0, KE, totDof, Eprim, 1, 0.0, aux, 1);
                uKu = E1 * cblas_ddot(totDof, aux, 1, Edual, 1);

                // sensitivity: wgt*alpha*dk
                ind = j*numVar + i; // place in mat_sens
                mat_sens[ind] = -alpha[inMesh->m_pMaterialElems[i]]*wgt[j]*uKu;

                // add sensitivty for self-weight
                if(sw)
                {
                    // sum ux+px, uy+py
                    sx=0.0; sy=0.0;
                    for(n=0;n<4;n++)
                    {
                        temp = n * NUM_DOF; // local x_dof
                        sx += Eprim[temp] + Edual[temp];
                        temp++; // local y_dof
                        sy += Eprim[temp] + Edual[temp];
                    }
                    mat_sens[ind] += alpha[inMesh->m_pMaterialElems[i]]*wgt[j]*rho*( (sx*acc[j].x) + (sy*acc[j].y) );
                }
            }
        }
        // if element is out use zero sensitivity
        else
        {
            // for each case
            for(j=0;j<numCase;j++)
            {
                ind = j*numVar + i; // place in mat_sens
                mat_sens[ind] = 0.0; // zero
            }
        }
	}

	free(aux);
	free(Eprim);
    free(Edual);
}

// function to compute designable material design varibles for eigenvalues
void CHakSensitivity::matSens_eig(CHakMesh *inMesh, CHakMaterial *inMat, double *KE, double *ME, double *mat_sens,
					int numEig, double *eig_vals, double **eig_vecs, double *alpha, double aMin)
{
	int i,j,n,m,ind,ind2,temp,temp2;
	int tnodes[4];
	double uKu, uMu;
	int m_elemX = inMesh->m_elemX;
	int numVar = inMesh->m_numDesignableMaterialt;
	Elem **Number = inMesh->m_pNumber;

	int totDof = 4*NUM_DOF; // total dof per element
	double *Eprim = (double *) malloc(totDof*sizeof(double));	// element eigenvector array
	double *aux = (double *) malloc(totDof*sizeof(double));

	// realtive material properties
	double E1 = inMat[inMesh->m_materialOne].m_e;
	double E2 = inMat[inMesh->m_materialTwo].m_e;
	double rho1 = inMat[inMesh->m_materialOne].m_rho;
	double rho2 = inMat[inMesh->m_materialTwo].m_rho;

	// sensitivity factors
	E1 = (E2 - E1)/E1;
	rho1 = (rho2 - rho1)/rho1;

	// for each design variable
	for(i=0;i<numVar;i++)
	{
        if(alpha[inMesh->m_pMaterialElems[i]] > aMin)
        {
            // compute element indices
            m = floor((double)inMesh->m_pMaterialElems[i] / (double)m_elemX);
            n = inMesh->m_pMaterialElems[i] - m*m_elemX;

            // Read in grid node numbers of current element
            tnodes[0] = Number[n][m].a;
            tnodes[1] = Number[n][m].b;
            tnodes[2] = Number[n][m].c;
            tnodes[3] = Number[n][m].d;

            // for each eigenvalue
            for(j=0;j<numEig;j++)
            {
                // get eigenvector for current element
                for(n=0;n<4;n++) // each node
                {
                    temp = n * NUM_DOF; // local
                    temp2 = tnodes[n] * NUM_DOF; // global

                    for(m=0;m<NUM_DOF;m++) // each dof
                    {
                        ind=temp+m; // local
                        ind2=temp2+m; // global
                        Eprim[ind] = eig_vecs[j][ind2];
                    }
                }

                // multiply dk = u^T K u
                cblas_dgemv(CblasRowMajor, CblasNoTrans, totDof, totDof, E1, KE, totDof, Eprim, 1, 0.0, aux, 1);
                uKu = cblas_ddot(totDof, aux, 1, Eprim, 1);

                // multiply dm = u^T M u
                cblas_dgemv(CblasRowMajor, CblasNoTrans, totDof, totDof, rho1, ME, totDof, Eprim, 1, 0.0, aux, 1);
                uMu = cblas_ddot(totDof, aux, 1, Eprim, 1);

                // sensitivity: alpha*(dk - eig x dm)
                ind = j*numVar + i; // place in mat_sens
                mat_sens[ind] = -alpha[inMesh->m_pMaterialElems[i]]*(uKu -eig_vals[j]*uMu); // reverse sign for maximize
            }
        }
        // if element is out use zero sensitivity
        else
        {
            // for each eigenvalue
            for(j=0;j<numEig;j++)
            {
                ind = j*numVar + i; // place in mat_sens
                mat_sens[ind] = 0.0; // zero
            }
        }
	}

	free(aux);
	free(Eprim);
}

// function to compute designable material design H-S varibles for eigenvalues
void CHakSensitivity::HS_Sens_eig(CHakMesh *inMesh, CHakMaterial *inMat, double *KE, double *ME, double *mat_sens,
                 int numEig, double *eig_vals, double **eig_vecs, double *alpha, double aMin)
{
	int i,j,n,m,ind,ind2,temp,temp2;
	int tnodes[4];
	double uKu, uMu;
	int elemX = inMesh->m_elemX;
	int numVar = inMesh->m_numDesignableMaterialt;
	Elem **Number = inMesh->m_pNumber;

	int totDof = 4*NUM_DOF; // total dof per element
	double *Eprim = (double *) malloc(totDof*sizeof(double));	// element eigenvector array
	double *aux = (double *) malloc(totDof*sizeof(double));

	// realtive material properties
	double E1 = inMat[inMesh->m_materialOne].m_e;
	double rho1 = inMat[inMesh->m_materialOne].m_rho;
	double rho2 = inMat[inMesh->m_materialTwo].m_rho;

	// mass sensitivity factor
	// rho1 = (rho1 - rho2)/rho1;
    rho1 = (rho2 - rho1)/rho1;

    double Ea; // current elastic modulus

	// for each design variable
	for(i=0;i<numVar;i++)
	{
        if(alpha[inMesh->m_pMaterialElems[i]] > aMin)
        {
            // compute element indices
            m = floor((double)inMesh->m_pMaterialElems[i] / (double)elemX);
            n = inMesh->m_pMaterialElems[i] - m*elemX;

            // Read in grid node numbers of current element
            tnodes[0] = Number[n][m].a;
            tnodes[1] = Number[n][m].b;
            tnodes[2] = Number[n][m].c;
            tnodes[3] = Number[n][m].d;

            // elastic modulus factor
            Ea = dE_dalpha(inMesh->m_pMaterialVars[i], 0.5, &inMat[inMesh->m_materialOne], &inMat[inMesh->m_materialTwo]);
            Ea /= E1;

            // for each eigenvalue
            for(j=0;j<numEig;j++)
            {
                // get eigenvector for current element
                for(n=0;n<4;n++) // each node
                {
                    temp = n * NUM_DOF; // local
                    temp2 = tnodes[n] * NUM_DOF; // global

                    for(m=0;m<NUM_DOF;m++) // each dof
                    {
                        ind=temp+m; // local
                        ind2=temp2+m; // global
                        Eprim[ind] = eig_vecs[j][ind2];
                    }
                }

                // multiply dk = u^T K u
                cblas_dgemv(CblasRowMajor, CblasNoTrans, totDof, totDof, Ea, KE, totDof, Eprim, 1, 0.0, aux, 1);
                uKu = cblas_ddot(totDof, aux, 1, Eprim, 1);

                // multiply dm = u^T M u
                cblas_dgemv(CblasRowMajor, CblasNoTrans, totDof, totDof, rho1, ME, totDof, Eprim, 1, 0.0, aux, 1);
                uMu = cblas_ddot(totDof, aux, 1, Eprim, 1);

                // sensitivity: alpha*(dk - eig x dm)
                ind = j*numVar + i; // place in mat_sens
                mat_sens[ind] = -alpha[inMesh->m_pMaterialElems[i]]*(uKu -eig_vals[j]*uMu); // reverse sign for maximize
            }
        }
        // if element is out use zero sensitivity
        else
        {
            // for each eigenvalue
            for(j=0;j<numEig;j++)
            {
                ind = j*numVar + i; // place in mat_sens
                mat_sens[ind] = 0.0; // zero
            }
        }
	}

	free(aux);
	free(Eprim);
}

// derivative of E for H-S bound material model
double CHakSensitivity::dE_dalpha(double alpha, double hs_int, CHakMaterial *mat1, CHakMaterial *mat2)
{
    double dKmax, dKmin, dK, dGmax, dGmin, dG;
    double kmax, kmin, gmax, gmin, khs, ghs;
    double fact, ftemp, ftemp2, ftemp3;
    double amin = 1.0-alpha;

    // bulk modulus
    fact = amin*mat1->m_k + alpha*mat2->m_k;
    ftemp = mat2->m_k - mat1->m_k;
    ftemp *= ftemp;
    ftemp *= amin*alpha;

    kmax = fact - (ftemp / (amin*mat2->m_k + alpha*mat1->m_k + mat2->m_g));
    kmin = fact - (ftemp / (amin*mat2->m_k + alpha*mat1->m_k + mat1->m_g));

    khs = hs_int*kmax + (1.0-hs_int)*kmin;

    // shear modulus
    fact = amin*mat1->m_g + alpha*mat2->m_g;
    ftemp = mat2->m_g - mat1->m_g;
    ftemp *= ftemp;
    ftemp *= amin*alpha;

    ftemp2 = (mat2->m_k * mat2->m_g) / (mat2->m_k + 2.0*mat2->m_g);
    gmax = fact - (ftemp / (amin*mat2->m_g + alpha*mat1->m_g + ftemp2));

    ftemp2 = (mat1->m_k * mat1->m_g) / (mat1->m_k + 2.0*mat1->m_g);
    gmin = fact - (ftemp / (amin*mat2->m_g + alpha*mat1->m_g + ftemp2));

    ghs = hs_int*gmax + (1.0-hs_int)*gmin;

    // bulk modulus derivatives
    fact = mat2->m_k - mat1->m_k;
    ftemp = amin*mat2->m_k + alpha*mat1->m_k + mat2->m_g; // denom
    ftemp2 = (1.0-2.0*alpha)*fact;
    ftemp3 = (alpha-alpha*alpha)*fact*fact;

    dKmax = 1.0 - (ftemp2/ftemp) - (ftemp3/(ftemp*ftemp));
    dKmax *= fact;

    ftemp = amin*mat2->m_k + alpha*mat1->m_k + mat1->m_g; // denom
    dKmin = 1.0 - (ftemp2/ftemp) - (ftemp3/(ftemp*ftemp));
    dKmin *= fact;

    // shear modulus derivatives
    fact = mat2->m_g - mat1->m_g;
    ftemp = amin*mat2->m_g + alpha*mat1->m_g;
    ftemp += (mat2->m_k * mat2->m_g) / (mat2->m_k + 2.0*mat2->m_g);

    ftemp2 = (1.0-2.0*alpha)*fact;
    ftemp3 = (alpha-alpha*alpha)*fact*fact;

    dGmax = 1.0 - (ftemp2/ftemp) - (ftemp3/(ftemp*ftemp));
    dGmax *= fact;

    ftemp = amin*mat2->m_g + alpha*mat1->m_g;
    ftemp += (mat1->m_k * mat1->m_g) / (mat1->m_k + 2.0*mat1->m_g);

    dGmin = 1.0 - (ftemp2/ftemp) - (ftemp3/(ftemp*ftemp));
    dGmin *= fact;

    // total derivatives
    dK = hs_int*dKmax + (1.0-hs_int)*dKmin;
    dG = hs_int*dGmax + (1.0-hs_int)*dGmin;

    ftemp2 = 3.0*khs + ghs;
    ghs /= ftemp2;
    khs /= ftemp2;
    ftemp = (9.0*ghs*ghs*dK) + (27.0*khs*khs*dG);

    return (ftemp);
}

//for hole insertion method
// calculate sensitivies using least squares of integration points for AFG method
void CHakSensitivity::AFG_Sens_hole(mesh *inMesh, boundary *bound_in, double *alpha, isoMat *inMat,  double *Nsens,
              double **prim, double **dual, int numDual, int numCase, double *wgt, Coord *gCoord,
              double aMin, int mode, double *fact, bool sw, Coord *acc, int h_count, int *h_index, int *h_EmapX, int *h_EmapY, int *h_posN, int *h_posE, double *h_Esens, double *h_Nsens)
{
    // read in mesh data
    double h = inMesh->h;
    double thk = inMesh->t;
    int NumNodes = inMesh->NumNodes;
    int elemX = inMesh->elemX;
    int elemY = inMesh->elemY;
    int NumElem = inMesh->NumElem;
    Elem **Number = inMesh->Number;
    Coord *NodeCoord = inMesh->NodeCoord;
    double rad = 2.0 * h; // hard coded for two elements around a node

    int irad = 2; //ceil(rad/h); // search space
    rad *=rad; // input squared dist to Lsens
    Coord *AuxNodes = bound_in->AuxNodes;
    int numBound = bound_in->NumBound;
    Bseg *Boundary = bound_in->Bound;
    int Ntot = NumNodes + bound_in->NumAux;
    double alMin = (aMin < 1.0e-6) ? 1.0e-6 : aMin; // numerical tollerance for sens calc
    if(mode==1){ alMin = (aMin < 1.0e-2) ? 1.0e-2 : aMin; } // increase for eigenvalue sensitivities
    else if(mode==2){alMin = 0.0;} // for compliant mechanisms

    int i,n,m,o,num,ind,node,ex,ey,xMin,xMax,yMin,yMax;
    double atemp;
    //if(mode==2){w1 = -1.0/wgt[0]; w2 = wgt[1]*w1*w1;} // extra weights for complinat mech disp constraint
    int Gcount = 0;		// varible to track number of points evaluated
    int tnodes[4];
    int gpoints = NumElem * 4; // number of sensitivity values (4 per element)
    double *gSens = calloc(gpoints*numDual,sizeof(double));	// define memory for gauss point sensitivities
    isoMat *mptr; // pointer for material

    int mat_count = 0; // count for designable material
    double mtemp; // variables for designable material modification
    isoMat fgm;

    // Step 1. Compute gauss point senstivities

    // For All Elements
    for(m=0;m<elemY;m++)
    {
        for(n=0;n<elemX;n++)
        {
            num = Number[n][m].n; // Element number
            atemp = alpha[num];
            if(atemp > alMin) // If Element isn't OUT
            {
                tnodes[0] = Number[n][m].a;
                tnodes[1] = Number[n][m].b;
                tnodes[2] = Number[n][m].c;
                tnodes[3] = Number[n][m].d;

                o = inMesh->mat_type[num]; // material number
                mptr = &inMat[o]; // point to correct material

                // modify matrices for designable materials
                if(inMesh->des_mat && num==inMesh->mat_elems[mat_count])
                {
                    // material mix variable
                    mtemp = inMesh->mat_vars[mat_count++];

                    // modulus
                    if(inMesh->mat_lin){ fgm.e =  (1.0-mtemp)*inMat[inMesh->mat1].e + mtemp*inMat[inMesh->mat2].e; }
                    else{ fgm.e = HS_mat(mtemp, 0.5, &inMat[inMesh->mat1], &inMat[inMesh->mat2]); }

                    // density
                    fgm.rho =  mtemp*inMat[inMesh->mat2].rho + (1.0-mtemp)*inMat[inMesh->mat1].rho;

                    // material prop matrix
                    mtemp = fgm.e / inMat[inMesh->mat1].e;
                    for(i=0;i<9;i++) {
                        //fgm.mat[i] = mtemp*inMat[inMesh->mat1].mat[i] + (1.0-mtemp)*inMat[inMesh->mat2].mat[i];
                        fgm.mat[i] = mtemp*inMat[inMesh->mat1].mat[i];
                    }

                    fgm.v = inMat[inMesh->mat1].v;
                    mptr = &fgm;
                }
                else
                {
                    mptr = &inMat[o]; // point to correct material (not FGM)
                }

                //calcualte sensitivity at Gauss points
                if(mode==1) // eigenvalue sensitivites
                {
                    GaEigSens_Q4(tnodes, prim, dual, atemp, h, thk, mptr, Gcount, gSens, numCase, wgt);
                }
                else if(mode==2) // compliant mechanism sensitivites
                {
                    double *gSens2 = calloc(8, sizeof(double));
                    double *gSens3 = calloc(4, sizeof(double));	// temp memory for gauss point sensitivities

                    GaSens_Q4(tnodes, prim, dual, atemp, h, thk, mptr, 0, gSens2, 1, 2, wgt, sw, acc);
                    GaSens_Q4(tnodes, &dual[1], &dual[1], atemp, h, thk, mptr, 0, gSens3, 1, 1, wgt, sw, acc);

                    // process sensitivities for compliant mech
                    for(i=0;i<4;i++)
                    {
                        ind = (Gcount+i)*3; // point to place in gSens
                        gSens[ind] = -fact[0]*gSens2[(2*i)+1] - fact[1]*gSens3[i];
                        gSens[ind+1] = gSens2[2*i]*fact[4] + fact[2]*gSens2[(2*i)+1] + fact[3]*gSens3[i];
                        gSens[ind+2] = gSens3[i]*fact[5];
                    }

                    free(gSens2);
                    free(gSens3);
                }
                else
                { GaSens_Q4(tnodes, prim, dual, atemp, h, thk, mptr, Gcount, gSens, numCase, numDual, wgt, sw, acc); }
            }
            else if(inMesh->des_mat && num==inMesh->mat_elems[mat_count]){mat_count++;} // keep counting
            Gcount += 4; // update point count
        }
    }

    // Step 2. Compute node sensitivities using least squares method

    double *sens_temp = malloc(numDual * sizeof(double)); // array for sens of each dual state

    // calculate smoothed sensitivities for all boundary (or non-OUT) nodes
    bool *done = calloc(Ntot, sizeof(bool));
    double ftemp;

    // first search through each boundary segment
    for(n=0;n<numBound;n++)
    {
        ey = Boundary[n].e / elemX;
        ex = Boundary[n].e - elemX*ey; // element indices
        xMin = ex-irad; xMax=ex+irad+1;
        yMin = ey-irad; yMax=ey+irad+1; // search limits
        xMin = (xMin < 0) ? 0 : xMin;
        yMin = (yMin < 0) ? 0 : yMin;
        xMax = (xMax > elemX) ? elemX : xMax;
        yMax = (yMax > elemY) ? elemY : yMax; // adjust search limits

        node=Boundary[n].n1; // 1st node number
        if(!done[node]) // if sens not already computed
        {
            ftemp = 1.0;
            if(node < NumNodes) // If node is a grid node
            {
                num = Lsens(&NodeCoord[node], xMax, xMin, yMax, yMin, alMin, alpha, rad, gCoord, gSens, Number, 4, numDual, sens_temp);
            }
            else // Otherwise node is an auxillary node
            {
                m = node - NumNodes; // position in AuxNodes array
                num = Lsens(&AuxNodes[m], xMax, xMin, yMax, yMin, alMin, alpha, rad, gCoord, gSens, Number, 4, numDual, sens_temp);
            }
            for(i=0;i<numDual;i++)
            {
                m = Ntot*i + node; // point to correct place in Nsens
                Nsens[m] = ftemp * sens_temp[i]; // multiply smoothed sensitivities by weight
            }
            done[node] = true;
        }
        node=Boundary[n].n2; // 2nd node number
        if(!done[node]) // if sens not already computed
        {
            ftemp = 1.0;
            if(node < NumNodes) // If node is a grid node
            {
                num = Lsens(&NodeCoord[node], xMax, xMin, yMax, yMin, alMin, alpha, rad, gCoord, gSens, Number, 4, numDual, sens_temp);
            }
            else // Otherwise node is an auxillary node
            {
                m = node - NumNodes; // position in AuxNodes array
                num = Lsens(&AuxNodes[m], xMax, xMin, yMax, yMin, alMin, alpha, rad, gCoord, gSens, Number, 4, numDual, sens_temp);
            }
            for(i=0;i<numDual;i++)
            {
                m = Ntot*i + node; // point to correct place in Nsens
                Nsens[m] = ftemp * sens_temp[i]; // multiply smoothed sensitivities by weight
            }
            done[node] = true;
        }
    }

    //Get the Hole node sensitvities
    int A,B,C,D,j,l,k;
    int ind1,ind2,ind3,ind4;
    double *max_sens = malloc(numDual*sizeof(double));
    double *min_sens = malloc(numDual*sizeof(double));

    //for(i=0;i<(NumNodes*numDual);i++){h_Nsens[i] = 0.0;}

    for(i=0;i<h_count;i++)
    {
        n = h_EmapX[i];
        m = h_EmapY[i];

        //printf("\n%i: n: %i, m: %i", i, n, m);

        num = Number[n][m].n;
        A  = Number[n][m].a;
        B  = Number[n][m].b;
        C  = Number[n][m].c;
        D  = Number[n][m].d;

        for(j=0;j<numDual;j++)
        {
            ind1 = (4*num)*numDual+j; // point to place in gSens
            ind2 = (4*num+1)*numDual+j;
            ind3 = (4*num+2)*numDual+j;
            ind4 = (4*num+3)*numDual+j;
            k = h_posE[j];
            l = h_posN[j];
            h_Esens[k+i] = 0.25*(gSens[ind1] + gSens[ind2] + gSens[ind3] + gSens[ind4]);
            if(h_index[A]==1){h_Nsens[l+A] += 0.25*h_Esens[k+i];}
            if(h_index[B]==1){h_Nsens[l+B] += 0.25*h_Esens[k+i];}
            if(h_index[C]==1){h_Nsens[l+C] += 0.25*h_Esens[k+i];}
            if(h_index[D]==1){h_Nsens[l+D] += 0.25*h_Esens[k+i];}
        }
    }

    for(j=0;j<numDual;j++)
    {
        k = h_posN[j];
        max_sens[j] = 0;
        min_sens[j] = 0;

        for(i=0;i<NumNodes;i++)
        {
            if(h_index[i]==1)
            {
                max_sens[j] = ((h_Nsens[k+i])>max_sens[j]) ? (h_Nsens[k+i]):max_sens[j];
                min_sens[j] = ((h_Nsens[k+i])<min_sens[j]) ? (h_Nsens[k+i]):min_sens[j];
            }
        }
        max_sens[j] = max_sens[j]-min_sens[j];
        //printf("\n Maxsens[%i] = %f", j, max_sens[j]);
    }

    for(j=0;j<numDual;j++)
    {
        k = h_posN[j];
        for(i=0;i<NumNodes;i++)
        {
            if(h_index[i]==1)
            {
                h_Nsens[k+i] = h_Nsens[k+i]/max_sens[j];
            }
            else{h_Nsens[k+i] =1;}
        }
    }

    // clear memory
    free(max_sens);
    free(min_sens);
    free(done);
    free(gSens);
    free(sens_temp);
}
