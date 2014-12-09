/*
	CSolver.h

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

#ifndef CSOLVER_H_
#define CSOLVER_H_

//
// Includes
//
#include "CommonTypes.h"
#include "CMathUtility.h"
//
// Pre-Definitions
//

//
// Statics
//

// declarations of MA57 Fortran functions (double precision)
void ma57id_(double*,int*);
void ma57ad_(int*,int*,int*,int*,int*,int*,int*,int*,int*,double*);
void ma57bd_(int*,int*,double*,double*,int*,int*,int*,int*,int*,int*,int*,double*,int*,double*);
void ma57cd_(int*,int*,double*,int*,int*,int*,int*,double*,int*,double*,int*,int*,int*,int*);

// declaration of ARPACK Fortran functions (double precision, general symmetric)
void dsaupd_(int*,char*,int*,char*,int*,double*,double*,int*,double*,int*,int*,int*,double*,double*,int*,int*);
void dseupd_(int*,char*,int*,double*,double*,int*,double*,char*,int*,char*,int*,double*,double*,int*,double*,int*,int*,int*,double*,double*,int*,int*);

class CSolver {
public:
	CSolver();
	virtual ~CSolver();

public:
	// Function that sets up and calls the MA57 to solve Ku = f
	void FE_Solve(sp_mat *Kg, int *dofMap, double *disp, int freeDof, int order, int numRhs);

	// Function that sets up and runs the MA57 Multifrontal Solver
	void solve(int n, int ne, int *irn, int *jcn, double *A, int n_rhs, double *rhs, int pinfo);

	// function to solve a linear least squares problem (by factorization)
	void d_lsLPK(int N, int M, int NRHS, double *A, double *b);

	// fucntion to solve a non-symmetric system (LU decomposition)
	int dun_slvLPK(char trans, int N, int NRHS, double *A, double *B);

	// function to solve a symmetric system (Cholesky factorization)
	int dsy_slvLPK(int N, int NRHS, double *A, double *B);

	// function to invert a matrix
	void din_LPK(double* A, int N);

	// function to solve a linear program using HSL revised Simplex method
	void la01bd_(int*,int*,int*,double*,double*,double*,double*,double*,int*,int*,int*,double*,int*);
	int lp_simplex(int n, int m, int nu, double *x, double *c, double *A, double *b, double *u, int pinfo);

	// function to solve linear prog using interior point method
	int LPsolve(int n, int m, int nu, double *x, double *c, double *A, double *b, double *u, int pinfo);

	// funciton to use ARPACK reverse communication to solve generalized eigenvalue problem (mode 3)
	int eig_solve(int nev_in, sp_mat *Kg, sp_mat *Mg, int n, int order, int *dofMap, double *vals, double *vecs, int pinfo);

	// multiply a sparse matrix by a vector
	void Msp_Vec(int n, sp_mat *mat, double *vin, double *vout);

public:
	CMathUtility GetMathUtility() {return m_cMathUtil;}
	void SetMathUtility(CMathUtility cMathUtil) { m_cMathUtil = cMathUtil; }

private:
	CMathUtility m_cMathUtil;
};

#endif /* CSOLVER_H_ */
