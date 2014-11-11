/*
 *  Solve.h
 *
 *  Created by Peter Dunning
 *	Version 5.4 @ 23/04/2014
 */

// declarations of MA57 Fortran functions (double precision)
void ma57id_(double*,int*);
void ma57ad_(int*,int*,int*,int*,int*,int*,int*,int*,int*,double*);
void ma57bd_(int*,int*,double*,double*,int*,int*,int*,int*,int*,int*,int*,double*,int*,double*);
void ma57cd_(int*,int*,double*,int*,int*,int*,int*,double*,int*,double*,int*,int*,int*,int*);

// declaration of ARPACK Fortran functions (double precision, general symmetric)
void dsaupd_(int*,char*,int*,char*,int*,double*,double*,int*,double*,int*,int*,int*,double*,double*,int*,int*);
void dseupd_(int*,char*,int*,double*,double*,int*,double*,char*,int*,char*,int*,double*,double*,int*,double*,int*,int*,int*,double*,double*,int*,int*);

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
