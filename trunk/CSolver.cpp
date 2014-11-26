/*
	CSolver.cpp

	Created on: Nov 24, 2014
	Author: jeehang

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

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "CMesh.h"
#include "CMathUtility.h"
#include "CSolver.h"

//
// Constructor / Destructor
//

CSolver::CSolver() {
	// TODO Auto-generated constructor stub

}

CSolver::~CSolver() {
	// TODO Auto-generated destructor stub
}

//
// Implementation - Interfaces
//

// Function that sets up and calls the MA57 to solve Ku = f
void CSolver::FE_Solve(sp_mat *Kg, int *dofMap, double *disp, int freeDof, int order, int numRhs)
{
	int i,j,temp,temp2;

	// copy disp array (input as the rhs load)
	temp2 = freeDof*numRhs;
	double *disp_temp = malloc(temp2*sizeof(double));
	cblas_dcopy(temp2, disp, 1, disp_temp, 1);

	// Solve the equation
	solve(freeDof, Kg->ne, Kg->irn, Kg->jcn, Kg->A, numRhs, disp_temp, 0);

	// Re-insert zero displacements into the displacement array
	for(i=0;i<order;i++)
	{
		if(dofMap[i] != -1)	// if dof not fixed copy accross the value
		{
			for(j=0;j<numRhs;j++)
			{
				temp = j * order; // to get to correct place in disp
				temp2 = j * freeDof; // to get to correct place in disp_temp
				disp[i+temp] = disp_temp[dofMap[i]+temp2];
			}
		}
		else	// Else dof is fixed! i.e zero displacement
		{
			for(j=0;j<numRhs;j++)
			{
				temp = j * order; // to get to correct place in disp
				disp[i+temp] = 0.0;
			}
		}
	}

	// clear memory
	free(disp_temp);
}

// Function that sets up and runs the MA57 Multifrontal Solver
void CSolver::solve(int n, int ne, int *irn, int *jcn, double *A, int n_rhs, double *rhs, int pinfo)
{
	// Lots of variables for MA57 fortran multi-frontal solver
	int lkeep,lfact,lifact,lrhs,lwork,liwork;
	int *keep,*ifact,*iwork;
	double *fact,*work;

	// set arrays for control data
	int icntl[20];
	double cntl[5];

	// initalise  control data arrays
	ma57id_(cntl,icntl);

	icntl[4] = pinfo; // info printout on screen, 0 = none, 4 = Full
	icntl[7] = 1; // don't restart if space for factorization runs out

	// setup infomation arrays
	int info[40];
	double rinfo[20];

	int job = 0; // solve linear system AX=B

	// setup workspace arrays
	lkeep = 7*n+2*ne+42;
	keep = malloc(lkeep * sizeof(int));
	liwork = 5 * n;
	iwork = malloc(liwork * sizeof(int));

	// Do Analysis
	ma57ad_(&n,&ne,irn,jcn,&lkeep,keep,iwork,icntl,info,rinfo);

	// re-size integer work array for next phase
	free(iwork);
	iwork = malloc(n * sizeof(int));

	// Setup Factorization arrays
	lfact = 2 * info[8];
	fact = malloc(lfact * sizeof(double));

	lifact = 2 * info[9];
	ifact = malloc(lifact * sizeof(int));

	// Do factorization
	ma57bd_(&n,&ne,A,fact,&lfact,ifact,&lifact,&lkeep,keep,iwork,icntl,cntl,info,rinfo);

	// Setup doubleing point workspace
	lwork = n * n_rhs;
	work  = malloc(lwork * sizeof(double));

	 // Do solve. NB: lrhs = n
	lrhs = n;
	ma57cd_(&job,&n,fact,&lfact,ifact,&lifact,&n_rhs,rhs,&lrhs,work,&lwork,iwork,icntl,info);

	// free workspace
	free(keep);
	free(iwork);
	free(work);
	free(fact);
	free(ifact);
}

// function to solve a linear least squares problem (by factorization)
void CSolver::d_lsLPK(int N, int M, int NRHS, double *A, double *b)
{
    int INFO;
    char trans = 'N';
    //get optimum block size
    /*int ispec = 1; // return optimal block size
    int n2,n3,n4; // dummy - not used
    char name [] = "dgel"; // name of routine
    char opts [] = "N";*/
	// optimal block size for N = 6
    int nBlock = 48; //ilaenv_(&ispec, name, opts, &N, &n2, &n3, &n4,4,1);

    // set work array
    int LWORK = (N < M) ? N : M;
    LWORK += (LWORK > NRHS) ? (LWORK*nBlock) : (NRHS*nBlock);
    double *WORK = malloc(LWORK*sizeof(double));

    // call the LAPACK routine
    dgels_(&trans, &M, &N, &NRHS, A, &M, b, &M, WORK, &LWORK, &INFO);

    free(WORK);

    if(INFO !=0) {
        if(INFO < 0) {
            printf("\nError in DGEL: %i argument had an illegal value",(int)-INFO);
        }
        else {
            printf("\nError in DGEL: the %i-th diagonal element of the triangular factor of A is zero, so that A does not have full rank; the least squares solution could not be computed.",(int)INFO);
        }
    }
}

// fucntion to solve a non-symmetric system (LU decomposition)
int CSolver::dun_slvLPK(char trans, int N, int NRHS, double *A, double *B)
{
	int i,j;
	if(N==1)
	{
		if(fabs(A[0]) < 1.0e-12){ return -1; }
		for(i=0;i<NRHS;i++)
		{
			B[i] /= A[0];
		}
	}
	else if(N==2)
	{
		double ftemp = A[0];
		double det = A[0]*A[3] - A[1]*A[2];
		if(fabs(det) < 1.0e-12){ return -1; }
		A[0] = A[3] / det;
		A[1] /= -det;
		A[2] /= -det;
		A[3] = ftemp / det;

		for(i=0;i<NRHS;i++)
		{
			j=2*i;
			ftemp = B[j];
			B[j] = A[0]*B[j] + A[1]*B[j+1];
			B[j+1] = A[2]*ftemp + A[3]*B[j+1];
		}
	}
	else if(N==3)
	{
		double Minv[9];
		double Btemp[3];

		Minv[0] = A[4]*A[8] - A[5]*A[7];
		Minv[3] = A[5]*A[6] - A[3]*A[8];
		Minv[6] = A[3]*A[7] - A[4]*A[6];
		Minv[1] = A[2]*A[7] - A[1]*A[8];
		Minv[4] = A[0]*A[8] - A[2]*A[6];
		Minv[7] = A[6]*A[1] - A[0]*A[7];
		Minv[2] = A[1]*A[5] - A[2]*A[4];
		Minv[5] = A[2]*A[3] - A[0]*A[5];
		Minv[8] = A[0]*A[4] - A[1]*A[3];

		double det = A[0]*Minv[0] + A[1]*Minv[3] + A[2]*Minv[6];
		if(fabs(det)<1.0e-20){det = (det < 0.0) ? -1.0e20 : 1.0e20;} // 1/det
		else {det = 1.0/det;}

		for(i=0;i<NRHS;i++)
		{
			for(j=0;j<3;j++)
			{
				Btemp[j] = B[i*3+j];
			}
			j = i*3;
			B[j++] = det*(Btemp[0]*Minv[0] + Btemp[1]*Minv[1] + Btemp[2]*Minv[2]);
			B[j++] = det*(Btemp[0]*Minv[3] + Btemp[1]*Minv[4] + Btemp[2]*Minv[5]);
			B[j]   = det*(Btemp[0]*Minv[6] + Btemp[1]*Minv[7] + Btemp[2]*Minv[8]);
		}
	}

	else
	{
		int *IPIV = malloc(N*sizeof(int));
		int INFO;

		dgetrf_(&N,&N,A,&N,IPIV,&INFO); // factroize

		if(INFO !=0) {
			if(INFO < 0) {
				printf("\nError in DGETRF: %i argument had an illegal value",(int)-INFO);
			}
			else {
				printf("\nError in DGETRF: the matrix is singular and its factorization could not be computed");
				free(IPIV);
				return -1;
			}
		}

		dgetrs_(&trans,&N,&NRHS,A,&N,IPIV,B,&N,&INFO); // solve

		if(INFO !=0) {
			if(INFO < 0) {
				printf("\nError in DGETRS: %i argument had an illegal value",(int)-INFO);
			}
			else {
				printf("\nError in DGETRS!");
			}
		}

		free(IPIV);
	}
	return 0;
}

// function to solve a symmetric system (Cholesky factorization)
int CSolver::dsy_slvLPK(int N, int NRHS, double *A, double *B)
{
	int i,j;

	if(N==1)
	{
		if(fabs(A[0]) < 1.0e-12){ return -1; }
		for(i=0;i<NRHS;i++)
		{
			B[i] /= A[0];
		}
	}
	else if(N==2)
	{
		double ftemp = A[0];
		double det = A[0]*A[3] - A[1]*A[2];
		if(fabs(det) < 1.0e-12)
        {
            return -1;
        }
        det = 1.0/det;
		A[0] = A[3] * det;
		A[1] *= -det;
		A[2] *= -det;
		A[3] = ftemp * det;

		for(i=0;i<NRHS;i++)
		{
			j=2*i;
			ftemp = B[j];
			B[j] = A[0]*B[j] + A[1]*B[j+1];
			B[j+1] = A[2]*ftemp + A[3]*B[j+1];
		}
	}
	else if(N==3)
	{
		double Minv[9];
		double Btemp[3];

		Minv[0] = A[4]*A[8] - A[5]*A[7];
		Minv[3] = A[5]*A[6] - A[3]*A[8];
		Minv[6] = A[3]*A[7] - A[4]*A[6];
		Minv[1] = A[2]*A[7] - A[1]*A[8];
		Minv[4] = A[0]*A[8] - A[2]*A[6];
		Minv[7] = A[6]*A[1] - A[0]*A[7];
		Minv[2] = A[1]*A[5] - A[2]*A[4];
		Minv[5] = A[2]*A[3] - A[0]*A[5];
		Minv[8] = A[0]*A[4] - A[1]*A[3];

		double det = A[0]*Minv[0] + A[1]*Minv[3] + A[2]*Minv[6];
		if(fabs(det)<1.0e-20){det = (det < 0.0) ? -1.0e20 : 1.0e20;} // 1/det
		else {det = 1.0/det;}

		for(i=0;i<NRHS;i++)
		{
			for(j=0;j<3;j++)
			{
				Btemp[j] = B[i*3+j];
			}
			j = i*3;
			B[j++] = det*(Btemp[0]*Minv[0] + Btemp[1]*Minv[1] + Btemp[2]*Minv[2]);
			B[j++] = det*(Btemp[0]*Minv[3] + Btemp[1]*Minv[4] + Btemp[2]*Minv[5]);
			B[j]   = det*(Btemp[0]*Minv[6] + Btemp[1]*Minv[7] + Btemp[2]*Minv[8]);
		}
	}

	else
	{
		int INFO;
		char uplo = 'U';

		// factorization
		dpotf2_(&uplo, &N, A, &N, &INFO);

		if(INFO !=0) {
			if(INFO < 0) {
				printf("\nError in DPOTRF: %i argument had an illegal value",(int)-INFO);
			}
			else {
				printf("\nError in DPOTRF: the matrix is not pos def and its factorization could not be computed");
				return -1;
				for(i=0;i<N;i++)
				{
					j=i*(N+1); // diagonal entry
					if(fabs(A[j]) < 1.0e-12){ A[j] = 1.0e-12; }
				}

			}
		}

		// solve
		dpotrs_(&uplo, &N, &NRHS, A, &N, B, &N, &INFO);

		if(INFO < 0) {
			printf("\nError in DPOTRS: %i argument had an illegal value",(int)-INFO);
		}
	}
	return 0;
}

// function to invert a matrix
void CSolver::din_LPK(double* A, int N)
{
    int *IPIV = malloc((N+1)*sizeof(int));
    int LWORK = N*N;
    double *WORK = malloc(LWORK*sizeof(double));
    int INFO;

    dgetrf_(&N,&N,A,&N,IPIV,&INFO);
    dgetri_(&N,A,&N,IPIV,WORK,&LWORK,&INFO);

    free(IPIV);
    free(WORK);
}

// function to solve a linear program using HSL revised Simplex method
int CSolver::lp_simplex(int n, int m, int nu, double *x, double *c, double *A, double *b, double *u, int pinfo)
{
	int i,j,k,err;
	double f;
	// set up and call HSL LA01
	int totCon = m + nu; // total number of constraints (equal + inequal)

	// combine b & u
	double *b_in = malloc(totCon*sizeof(double));
    //double b_in[3];
	cblas_dcopy(nu, u, 1, b_in, 1);
	cblas_dcopy(m, b, 1, &b_in[nu], 1);

	// add I to A (for inequal constraints)
	double *A_in = calloc(totCon*n,sizeof(double));
    //double A_in[6];
    //for(i=0;i<6;i++){A_in[i]=0.0;}
	for(i=0;i<nu;i++)
	{
		j=i*totCon + i; // identity matrix
		A_in[j] = 1.0;
	}
	// transpose for FORTRAN column major order
	for(i=0;i<m;i++) // equal constraint num
	{
		k = nu + i;
		for(j=0;j<n;j++) // variable num
		{
			A_in[(j*totCon) + k] = A[(i*n) + j];
		}
	}

	// other arrays
	j=(totCon+1)*(totCon+3);
	double *work = malloc(j*sizeof(double));
	int *ind = malloc((totCon+1)*sizeof(int));

	// call the routine
    int db = pinfo;
	la01bd_(&n,&totCon,&totCon,A_in,b_in,c,x,&f,&totCon,&db,ind,work,&err);

	// clear memory
	free(b_in);
	free(A_in);
	free(work);
	free(ind);

	if(err != 0){return -1;}

	return err; // return error flag
}

// function to solve linear prog using interior point method
int CSolver::LPsolve(int n, int m, int nu, double *x, double *c, double *A, double *b, double *u, int pinfo)
{
	double *delx = malloc(n*sizeof(double));
	// array of slack varibles (s, length nu)
	double *s, *dels;
	if(nu > 0){s = malloc(nu*sizeof(double));
		dels = malloc(nu*sizeof(double));}
	else {s=0; dels=0;}
	// array of adjoint variables (y, length m)
	double *y = malloc(m*sizeof(double));
	double *dely = malloc(m*sizeof(double));
	// more adjoint variables (z, length n)
	double *z = malloc(n*sizeof(double));
	double *delz = malloc(n*sizeof(double));
	// array of adjoint slack variables (w, length nu)
	double *w, *delw;
	if(nu > 0){w = malloc(nu*sizeof(double));
		delw = malloc(nu*sizeof(double));}
	else {w=0; delw=0;}
	// assemble all variables into one array using pointers
	double *v[5] = {x,z,s,w,y};
	double *delv[5] = {delx,delz,dels,delw,dely};

	// compute lengths of b, c & u arrays (used in stopping criterion)
	double lenb = cblas_dnrm2(m, b, 1); // L2 norm of b
	lenb = (1.0 > lenb) ? 1.0 : lenb;

	double lenc = cblas_dnrm2(n, c, 1); // L2 norm of c
	lenc = (1.0 > lenc) ? 1.0 : lenc;

	double lenu = 0.0;
	if(nu>0) {
		cblas_dnrm2(nu, u, 1); // L2 norm of u
	}
	lenu = (1.0 > lenu) ? 1.0 : lenu;

    //  --- scale input gradients --- //
    /*int i,j;
    double *maxG = calloc(m+1, sizeof(double));
    double ftemp;

    // find max gradients
    for(i=0;i<(n-m);i++) // each variable
    {
        ftemp = fabs(c[i]);
        maxG[0] = (ftemp > maxG[0]) ? ftemp : maxG[0]; //objective
        for(j=0;j<m;j++) // each constraint
        {
            ftemp = fabs(A[j*n + i]);
            maxG[j+1] = (ftemp > maxG[j+1]) ? ftemp : maxG[j+1];
        }
    }

    if(pinfo==3)
    {
        for(j=0;j<=m;j++){printf("\nmaxG[%i]=%12.4e",j,maxG[j]);}
    }

    // now scale to one
    for(i=0;i<n;i++) // each variable
    {
        c[i] /= maxG[0]; // objective
        for(j=0;j<m;j++) // each constraint
        {
            A[j*n + i] /= maxG[j+1];
        }
    }
    for(j=0;j<m;j++) // each constraint value
    {
        b[j] /= maxG[j+1];
    }*/

	// Initialize variables: v = (x,z,s,w,y)
	// ---------------------
	int err = GetMathUtility().initialize(n, m, nu, v, A, b, u, c);
	if(err != 0) {
		// Clean up and exit

        // re-scale gradients
        /*for(i=0;i<(n-m);i++) // each variable
        {
            c[i] *= maxG[0]; // objective
            for(j=0;j<m;j++) // each constraint
            {
                A[j*n + i] *= maxG[j+1];
            }
        }
        for(j=0;j<m;j++) // each constraint value
        {
            b[j] *= maxG[j+1];
        }*/

		if(nu>0)
		{
			free(s);
			free(w);
			free(dels);
			free(delw);
		}
		free(y);
		free(z);
		free(dely);
		free(delz);
		free(delx);
        //free(maxG);
		return -1; // exit on error
	}

	//printf("\nReport initial values");
	//report(n, m, nu, v);

	// Main Algorithm
	// ---------------
	// Do until ||F(v)|| is small
	double mu, lim, stop;
    double stop_min = 1.0e20;
	lim = 0.2;
	int count = 0;
	do {
		//printf("\nIteration %i -------\n",count+1);
		// Compute delP(v) = -F'(v)F(v) : special function
		err = GetMathUtility().predictor(n, m, nu, v, delv, A, b, u, c);
		if (err==-1) {
			// Clean up and exit

            // re-scale gradients
            /*for(i=0;i<n;i++) // each variable
            {
                c[i] *= maxG[0]; // objective
                for(j=0;j<m;j++) // each constraint
                {
                    A[j*n + i] *= maxG[j+1];
                }
            }
            for(j=0;j<m;j++) // each constraint value
            {
                b[j] *= maxG[j+1];
            }*/

			if(nu>0)
			{
				free(s);
				free(w);
				free(dels);
				free(delw);
			}
			free(y);
			free(z);
			free(dely);
			free(delz);
			free(delx);
            //free(maxG);
			return -1; // exit on error
		}

		// determine mu (special fucntion)
		mu = GetMathUtility().centering(n, m, nu, lim, v, delv);

		// Compute delC(v) = -F'(v)F(v+delP(v)) - mu*e (e = 0 or 1)
		GetMathUtility().corrector(n, nu, mu, v, delv);

		// choose alpha > 0 : special function linked to update step
		// update v(k+1) = v(k) + alpha( delP(v) + delC(v) )
		// step length to ensure v(k+1) > 0
		GetMathUtility().update(n, nu, m, v, delv);

		// Assess stopping crieria : special function
		stop = GetMathUtility().stopping(n, m, nu, v, A, lenb, b, lenu, u, lenc, c);
        stop_min = (stop < stop_min) ? stop : stop_min;
		count++;

		if(stop > 2.0*stop_min && stop > 0.1)
		{
			//printf("\nReporting v - update");
			if(pinfo==3)
			{
				COutput cOutput;
				cOutput.report(n, m, nu, v);
			}
			count = 200;
		}

		lim = (stop < 0.2) ? stop : 0.2;
		//if(dual_gap(n,nu,v) < TINY){lim=TINY;}
		if(pinfo==3){printf("\n%i itt, function = %12.4e, constraint = %12.4e, stop = %12.4e",count,cblas_ddot(n, v[0], 1, c, 1),cblas_ddot(n, v[0], 1, A, 1),stop);}
	} while (stop > TINY && count < 200);

	// report solution
	if(pinfo==3 && stop < TINY){
		printf("\nSolution found after %i itterations",count);
		//report(n, m, nu, v);
	}

    // re-scale gradients
    /*for(i=0;i<(n-m);i++) // each variable
    {
        c[i] *= maxG[0]; // objective
        for(j=0;j<m;j++) // each constraint
        {
            A[j*n + i] *= maxG[j+1];
        }
    }
    for(j=0;j<m;j++) // each constraint value
    {
        b[j] *= maxG[j+1];
    }*/

	// Clean up and exit
	if(nu>0)
	{
		free(s);
		free(w);
		free(dels);
		free(delw);
	}
	free(y);
	free(z);
	free(dely);
	free(delz);
	free(delx);
    //free(maxG);

	if( fabs(b[0] - cblas_ddot(n, v[0], 1, A, 1)) > TINY){return -1;}
	//if(stop > TINY){ return -1; }
	return 0;
}

// funciton to use ARPACK reverse communication to solve generalized eigenvalue problem (mode 3)
int CSolver::eig_solve(int nev_in, sp_mat *Kg, sp_mat *Mg, int n, int order, int *dofMap, double *vals, double *vecs, int pinfo)
{
	int nev = 2 * nev_in; // always compute twice the required amount (for repeated values)
	int i,j;
	double sigma = 0.0; // zero shift (want eig vals closest to zero)

	// NB: [Kg - sigma*Mg] = Kg (sigma = 0)
	// then factorize Kg using MA57 (use for solve ... )

	// Lots of variables for MA57 fortran multi-frontal solver
	int num_ent = Kg->ne;
	int lkeep,lfact,lifact,lwork,liwork;
	int *keep,*ifact,*iwork;
	double *fact,*work;

	// set arrays for control data
	int icntl[20];
	double cntl[5];

	// initalise  control data arrays
	ma57id_(cntl,icntl);

	icntl[4] = 0; // info printout on screen, 0 = none, 4 = Full
	icntl[7] = 1; // don't restart if space for factorization runs out

	// setup infomation arrays
	int info[40];
	double rinfo[20];

	int n_rhs = 1; // only solve for 1 rhs at a time
	int job = 0; // solve linear system AX=B

	// setup workspace arrays
	lkeep = 7*n+2*num_ent+42;
	keep = malloc(lkeep * sizeof(int));
	liwork = 5 * n;
	iwork = malloc(liwork * sizeof(int));

	// Do Analysis
	ma57ad_(&n,&num_ent,Kg->irn,Kg->jcn,&lkeep,keep,iwork,icntl,info,rinfo);

	// re-size integer work array for next phase
	free(iwork);
	iwork = malloc(n * sizeof(int));

	// Setup Factorization arrays
	lfact = 2 * info[8];
	fact = malloc(lfact * sizeof(double));

	lifact = 2 * info[9];
	ifact = malloc(lifact * sizeof(int));

	// Do factorization
	ma57bd_(&n,&num_ent,Kg->A,fact,&lfact,ifact,&lifact,&lkeep,keep,iwork,icntl,cntl,info,rinfo);

	// Setup doubleing point workspace
	lwork = n;
	work  = malloc(lwork * sizeof(double));

	int ido = 0; // initial value
	int ncv = 2*nev +2; // default
	double tol = 0.0; // default (use machine precision)
	char which [] = "LM";
	char bmat = 'G'; // generalized system
	double *resid = malloc(n*sizeof(double)); // residual
	double *v = malloc(n*ncv*sizeof(double)); // Lanczos basis vectors
	int ldv = n; // leading dimension of v (n * nev)
	int ipntr[11]; // locations of various arrays used "under the hood"
	double *workd = malloc(3*n*sizeof(double));
	int lworkl = ncv*(8+ncv); // length of lwork
	double *workl = malloc(lworkl*sizeof(double));
	//int ar_info = 0; // start with random residual vector
    int ar_info = 1; // start with set residual vector
    for(i=0;i<n;i++){ resid[i] = 1.0; }

	int iparam[11]; // parameters
	iparam[0] = 1; // ISHIFT = 1 - exact shifts
	iparam[2] = 100*nev; // max number of iterations (default)
	iparam[3] = 1; // block size
	iparam[6] = 3; // MODE = 3 (shift invert, as M assumed semi-pos-def)

	// iterate to find eigenvalues
	int stop=0;
	do
	{
		dsaupd_(&ido, &bmat, &n, which, &nev, &tol, resid, &ncv, v, &ldv, iparam, ipntr, workd, workl, &lworkl, &ar_info);
		if(ido == -1)
		{
			Msp_Vec(n, Mg, &workd[ipntr[0]-1], &workd[ipntr[1]-1]);
			ma57cd_(&job,&n,fact,&lfact,ifact,&lifact,&n_rhs,&workd[ipntr[1]-1],&n,work,&lwork,iwork,icntl,info);
		}
		else if(ido == 1)
		{
			// copy workd[ipntr[2]] -> workd[ipntr[1]]
			cblas_dcopy(n, &workd[ipntr[2]-1], 1, &workd[ipntr[1]-1], 1);
			ma57cd_(&job,&n,fact,&lfact,ifact,&lifact,&n_rhs,&workd[ipntr[1]-1],&n,work,&lwork,iwork,icntl,info);
		}
		else if (ido == 2)
		{
			Msp_Vec(n, Mg, &workd[ipntr[0]-1], &workd[ipntr[1]-1]);
		}
		else{stop=1;}
	} while(stop==0);

	// now post-process to obtain eigenvalues & vectors
	int rvec = 1; // compute vectors and values
	int *select = malloc(ncv*sizeof(int));
	for(i=0;i<ncv;i++){select[i]=1;}
	double *vec_temp = malloc(n*nev*sizeof(double)); // temp array
	dseupd_(&rvec,"All",select,vals,vec_temp,&n,&sigma,&bmat,&n,which,&nev,&tol,resid,&ncv,v,&ldv,iparam,ipntr,workd,workl,&lworkl,&ar_info);

	if(pinfo>1) {
        for(i=0;i<nev;i++) { printf("\nEig %i = %12.12e",i+1,vals[i]); }
		for(i=0;i<nev;i++) { printf("\nFreq %i = %12.4e",i+1,sqrt(vals[i])); }
	}

	// Re-insert zero displacements into the eigenvector arrays
	int temp,temp2;

	for(i=0;i<order;i++)
	{
		if(dofMap[i] != -1)	// if dof not fixed copy accross the value
		{
			for(j=0;j<nev;j++)
			{
				temp = j * order; // to get to correct place in disp
				temp2 = j * n; // to get to correct place in disp_temp
				vecs[i+temp] = vec_temp[dofMap[i]+temp2]; // / scale;
			}
		}
		else	// Else dof is fixed! i.e zero displacement
		{
			for(j=0;j<nev;j++)
			{
				temp = j * order; // to get to correct place in vecs
				vecs[i+temp] = 0.0;
			}
		}
	}

	// clear memory
	free(select);
	free(vec_temp);
	free(v);
	free(workd);
	free(workl);
	free(resid);
	free(keep);
	free(iwork);
	free(work);
	free(fact);
	free(ifact);

	return ar_info;
}

// multiply a sparse matrix by a vector
void CSolver::Msp_Vec(int n, sp_mat *mat, double *vin, double *vout)
{
	int i,rn,cn;
	int ne = mat->ne;
	// reset vout to zero
	for(i=0;i<n;i++){vout[i]=0.0;}

	// multiply
	for(i=0;i<ne;i++)
	{
		rn = mat->irn[i]-1;
		cn = mat->jcn[i]-1;

		vout[rn] += mat->A[i] * vin[cn];

		// symmetry
		if(rn != cn)
		{
			vout[cn] += mat->A[i] * vin[rn];
		}
	}
}

