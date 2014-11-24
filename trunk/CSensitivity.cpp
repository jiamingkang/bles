/*
	CSensitivity.cpp

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

#include "CSensitivity.h"

//
// Constructor & Destructor
//
CSensitivity::CSensitivity() {
	// TODO Auto-generated constructor stub

}

CSensitivity::~CSensitivity() {
	// TODO Auto-generated destructor stub
}

//
// Implementation - Interfaces
//

// Function that calculates the sensitivity of a node by a least squares (2nd order) filter of near-by gauss points
int CSensitivity::Lsens(Coord *pt, int xMax, int xMin, int yMax, int yMin, double aMin, double *alpha, double r2,
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

	double *A = malloc( 6 * count * sizeof(double)); // array to store gpoint values
	double *B = malloc(numDual * count * sizeof(double)); // array to store rhs' of equation
	double *Bmax =  calloc(numDual, sizeof(double)); // array to store maximum sens values
	double *Bmin =  calloc(numDual, sizeof(double)); // array to store minimum sens values

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
	d_lsLPK(6, count, numDual, A, B);

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
void CSensitivity::GaSens_Q4(int *tnodes, double **prim, double **dual, double alpha, double h, double t,
				isoMat *inMat, int Gcount, double *gSens, int numCase, int numDual, double *wgt, bool sw, Coord *acc)
{
	int i,j,k,n,p,temp,temp2,ind,ind2;	// incrementors etc
	int totDof = 4*NUM_DOF; // total dof per element
	double *Eprim = malloc(totDof*sizeof(double));	// Primary Element displacement array
	double *Edual = malloc(totDof*numDual*sizeof(double));	// Adjoint Element displacement array
	double stress[3], strain[3], sens;	 // Variabes for strain tensors
    Coord ug;
    Coord *pg = malloc(numDual*sizeof(Coord)); // primal & dual displacements at gauss point

	// copy material property matrix
	double Emat[9];
	for(i=0;i<9;i++){Emat[i] = inMat->mat[i];}

	double h_fact = (t*alpha) / (4.0*h*h); // factor to re-dimension senstivity (inc area ratio & thickness)
    double sw_fact;
    //if(sw){ sw_fact = alpha * inMat->rho * t; } // const factor for self-weight
    if(sw){ sw_fact = inMat->rho * t; } // const factor for self-weight

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
void CSensitivity::GaEigSens_Q4(int *tnodes, double **prim, double **dual, double alpha, double h, double t,
					isoMat *inMat, int Gcount, double *gSens, int num_eig, double *eig)
{
	int i,j,k,p,temp,temp2,ind,ind2;	// incrementors etc
	int totDof = 4*NUM_DOF; // total dof per element
	double *Eprim = malloc(totDof*sizeof(double));	// Primary Element displacement array
	double *Edual = malloc(totDof*sizeof(double));	// Adjoint Element displacement array
	double stress[3], strain[3], sens, ftemp;	 // Variabes for strain tensors
	double Iprim, Idual; // variables for interpolated displacements

	// copy material property matrix
	double Emat[9];
	for(i=0;i<9;i++){Emat[i] = inMat->mat[i];}

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

			sens -= eig[k]*ftemp*inMat->rho*alpha*t; // extra factor

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
void CSensitivity::barSens(mesh *inMesh, double *bar_sens, double **prim, int numCase, double *wgt)
{
	// read data
	int numBars = inMesh->NumBars;
	double fact = inMesh->bar_mat->e / inMesh->h;
	Bseg *bar = inMesh->bar_nums;

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
void CSensitivity::bcSens(mesh *inMesh, double *bc_sens, double **prim, int numCase, double *wgt)
{
	int n,m,o,c,i,j,temp,ind;
	int elemX = inMesh->elemX;
	int numBC = inMesh->NumBC;
	Elem **Number = inMesh->Number;
	int tnodes[4];
	double ftemp,disp;

	// for all
	for(o=0;o<numBC;o++)
	{
		// compute element indices
		m = floor((double)inMesh->BC_nums[o] / (double)elemX);
		n = inMesh->BC_nums[o] - m*elemX;

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
			ftemp = inMesh->K_bc[o]; // current design variable
			bc_sens[o] -= wgt[c] * disp * 3.0*ftemp*ftemp * inMesh->K0_bc;
		}
	}
}

// function to compute designable material design varibles for compliance
void CSensitivity::matSens_comp(mesh *inMesh, isoMat *inMat, double *KE, double *mat_sens, int numCase, double *wgt,
                    double **prim, double **dual, double *alpha, double aMin, bool sw, Coord *acc)
{
    int i,j,n,m,ind,ind2,temp,temp2;
	int tnodes[4];
	double uKu, sx, sy;
	int elemX = inMesh->elemX;
	int numVar = inMesh->NumDesMat;
	Elem **Number = inMesh->Number;

	int totDof = 4*NUM_DOF; // total dof per element
	double *Eprim = malloc(totDof*sizeof(double));	// primary element disp array
    double *Edual = malloc(totDof*sizeof(double));	// primary element disp array
	double *aux = malloc(totDof*sizeof(double));

	// realtive material properties
	double E1 = inMat[inMesh->mat1].e;
	double E2 = inMat[inMesh->mat2].e;
    double rho = inMat[inMesh->mat2].rho - inMat[inMesh->mat1].rho;
    rho *= 0.25 * inMesh->t * inMesh->h * inMesh->h; // multiply by volume to get mass / node

	// sensitivity factor
	E1 = (E2 - E1)/E1;

    // for each design variable
	for(i=0;i<numVar;i++)
	{
        if(alpha[inMesh->mat_elems[i]] > aMin)
        {
            // compute element indices
            m = floor((double)inMesh->mat_elems[i] / (double)elemX);
            n = inMesh->mat_elems[i] - m*elemX;

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
                mat_sens[ind] = -alpha[inMesh->mat_elems[i]]*wgt[j]*uKu;

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
                    mat_sens[ind] += alpha[inMesh->mat_elems[i]]*wgt[j]*rho*( (sx*acc[j].x) + (sy*acc[j].y) );
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
void CSensitivity::matSens_eig(mesh *inMesh, isoMat *inMat, double *KE, double *ME, double *mat_sens,
					int numEig, double *eig_vals, double **eig_vecs, double *alpha, double aMin)
{
	int i,j,n,m,ind,ind2,temp,temp2;
	int tnodes[4];
	double uKu, uMu;
	int elemX = inMesh->elemX;
	int numVar = inMesh->NumDesMat;
	Elem **Number = inMesh->Number;

	int totDof = 4*NUM_DOF; // total dof per element
	double *Eprim = malloc(totDof*sizeof(double));	// element eigenvector array
	double *aux = malloc(totDof*sizeof(double));

	// realtive material properties
	double E1 = inMat[inMesh->mat1].e;
	double E2 = inMat[inMesh->mat2].e;
	double rho1 = inMat[inMesh->mat1].rho;
	double rho2 = inMat[inMesh->mat2].rho;

	// sensitivity factors
	E1 = (E2 - E1)/E1;
	rho1 = (rho2 - rho1)/rho1;

	// for each design variable
	for(i=0;i<numVar;i++)
	{
        if(alpha[inMesh->mat_elems[i]] > aMin)
        {
            // compute element indices
            m = floor((double)inMesh->mat_elems[i] / (double)elemX);
            n = inMesh->mat_elems[i] - m*elemX;

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
                mat_sens[ind] = -alpha[inMesh->mat_elems[i]]*(uKu -eig_vals[j]*uMu); // reverse sign for maximize
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
void CSensitivity::HS_Sens_eig(mesh *inMesh, isoMat *inMat, double *KE, double *ME, double *mat_sens,
                 int numEig, double *eig_vals, double **eig_vecs, double *alpha, double aMin)
{
	int i,j,n,m,ind,ind2,temp,temp2;
	int tnodes[4];
	double uKu, uMu;
	int elemX = inMesh->elemX;
	int numVar = inMesh->NumDesMat;
	Elem **Number = inMesh->Number;

	int totDof = 4*NUM_DOF; // total dof per element
	double *Eprim = malloc(totDof*sizeof(double));	// element eigenvector array
	double *aux = malloc(totDof*sizeof(double));

	// realtive material properties
	double E1 = inMat[inMesh->mat1].e;
	double rho1 = inMat[inMesh->mat1].rho;
	double rho2 = inMat[inMesh->mat2].rho;

	// mass sensitivity factor
	// rho1 = (rho1 - rho2)/rho1;
    rho1 = (rho2 - rho1)/rho1;

    double Ea; // current elastic modulus

	// for each design variable
	for(i=0;i<numVar;i++)
	{
        if(alpha[inMesh->mat_elems[i]] > aMin)
        {
            // compute element indices
            m = floor((double)inMesh->mat_elems[i] / (double)elemX);
            n = inMesh->mat_elems[i] - m*elemX;

            // Read in grid node numbers of current element
            tnodes[0] = Number[n][m].a;
            tnodes[1] = Number[n][m].b;
            tnodes[2] = Number[n][m].c;
            tnodes[3] = Number[n][m].d;

            // elastic modulus factor
            Ea = dE_dalpha(inMesh->mat_vars[i], 0.5, &inMat[inMesh->mat1], &inMat[inMesh->mat2]);
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
                mat_sens[ind] = -alpha[inMesh->mat_elems[i]]*(uKu -eig_vals[j]*uMu); // reverse sign for maximize
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
{
    double dKmax, dKmin, dK, dGmax, dGmin, dG;
    double kmax, kmin, gmax, gmin, khs, ghs;
    double fact, ftemp, ftemp2, ftemp3;
    double amin = 1.0-alpha;

    // bulk modulus
    fact = amin*mat1->k + alpha*mat2->k;
    ftemp = mat2->k - mat1->k;
    ftemp *= ftemp;
    ftemp *= amin*alpha;

    kmax = fact - (ftemp / (amin*mat2->k + alpha*mat1->k + mat2->g));
    kmin = fact - (ftemp / (amin*mat2->k + alpha*mat1->k + mat1->g));

    khs = hs_int*kmax + (1.0-hs_int)*kmin;

    // shear modulus
    fact = amin*mat1->g + alpha*mat2->g;
    ftemp = mat2->g - mat1->g;
    ftemp *= ftemp;
    ftemp *= amin*alpha;

    ftemp2 = (mat2->k * mat2->g) / (mat2->k + 2.0*mat2->g);
    gmax = fact - (ftemp / (amin*mat2->g + alpha*mat1->g + ftemp2));

    ftemp2 = (mat1->k * mat1->g) / (mat1->k + 2.0*mat1->g);
    gmin = fact - (ftemp / (amin*mat2->g + alpha*mat1->g + ftemp2));

    ghs = hs_int*gmax + (1.0-hs_int)*gmin;

    // bulk modulus derivatives
    fact = mat2->k - mat1->k;
    ftemp = amin*mat2->k + alpha*mat1->k + mat2->g; // denom
    ftemp2 = (1.0-2.0*alpha)*fact;
    ftemp3 = (alpha-alpha*alpha)*fact*fact;

    dKmax = 1.0 - (ftemp2/ftemp) - (ftemp3/(ftemp*ftemp));
    dKmax *= fact;

    ftemp = amin*mat2->k + alpha*mat1->k + mat1->g; // denom
    dKmin = 1.0 - (ftemp2/ftemp) - (ftemp3/(ftemp*ftemp));
    dKmin *= fact;

    // shear modulus derivatives
    fact = mat2->g - mat1->g;
    ftemp = amin*mat2->g + alpha*mat1->g;
    ftemp += (mat2->k * mat2->g) / (mat2->k + 2.0*mat2->g);

    ftemp2 = (1.0-2.0*alpha)*fact;
    ftemp3 = (alpha-alpha*alpha)*fact*fact;

    dGmax = 1.0 - (ftemp2/ftemp) - (ftemp3/(ftemp*ftemp));
    dGmax *= fact;

    ftemp = amin*mat2->g + alpha*mat1->g;
    ftemp += (mat1->k * mat1->g) / (mat1->k + 2.0*mat1->g);

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
