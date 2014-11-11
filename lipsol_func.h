/*
 *  lipsol_func.h
 *
 *  Created by Peter Dunning on 22/06/2013.
 *	Version 5.4 @ 23/04/2014
 *
 */

#define TINY 1.0E-5
#define SMALL 1.0E-8
#define BIG 1.0E6

// function to handle divide by zero
double divZero(double num, double denom);

// function to find an initial point
int initialize(int n,  int m, int nu, double **v, double *A, double *b, double *u, double *c);

int initialize2(int n,  int m, int nu, double **v, double *A, double *b, double *u, double *c);

// function to compute del(v) - predictor
int predictor(int n,  int m, int nu, double **v, double **delP, double *A, double *b, double *u, double *c);

// function to compute delC(v) - corrector
void corrector(int n, int nu, double mu, double **v, double **delP);

// function to compute the centering parameter
double centering(int n, int m, int nu, double lim, double **v, double **delP);

// fucntion to find the duality gap
double dual_gap(int n, int nu, double **v);

// update the variables
void update(int n, int nu, int m, double **v, double **delP);

// compute the stopping criterion
double stopping(int n,  int m, int nu, double **v, double *A, double lenb, double *b,
			 double lenu, double *u, double lenc, double *c);

// function to print solution to screen
void report(int n,  int m, int nu, double **v);
