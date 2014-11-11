/*
 *  lipsol_func.c
 *
 *  Created by Peter Dunning on 22/06/2013.
 *	Version 5.4 @ 23/04/2014
 *
 */
#include "ls_types.h"
#include <stdlib.h>
#include <stdio.h>
#include "lipsol_func.h"
#include "Solve.h"

// function to handle divide by zero
double divZero(double num, double denom)
{
	if(fabs(denom) < SMALL) {
		return 1.0E+8;
	}
	else {
		return (num / denom);
	}
}

// function to find an initial point
int initialize(int n,  int m, int nu, double **v, double *A, double *b, double *u, double *c)
{
	int i;
	double dtemp;
	
	// first compute AA^T matrix and initial y and x values
	double *AAT = malloc(m*m*sizeof(double));
	
	if(m==1)
	{
		// AAT is simply the dot product
		dtemp = cblas_ddot(n, A, 1, A, 1);
		
		// initial y values	
		v[4][0] = cblas_ddot(n, A, 1, c, 1) / dtemp;
		
		// initial x values (primary variables)
		dtemp = *b / dtemp;
		for(i=0;i<n;i++)
		{
			v[0][i] = A[i] * dtemp;
		}
		
		// initial z values (dual variables)
		dtemp = v[4][0]; // y value
		for(i=0;i<n;i++)
		{
			v[1][i] = c[i] - A[i] * dtemp; // c - AT * y 
		}
		
	}
	else
	{
		char uplo = 'U';
		int one = 1;
		int info;
		// otherwise perform matrix matrix multiplication
		cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, m, m, n, 1.0, A, n, A, n, 0.0, AAT, m);
		
		// factor the matrix ready for inversion
		dpotrf_(&uplo, &m, AAT, &m, &info);
		
		if(info < 0) {
			printf("\nError in dpotrf - argument %i has illegal value - aborting",-info);
			free(AAT);
			return -1;
		}
		else if (info > 0) {
			printf("\nError in dpotrf - the leading minor of order %i is not positive definite,\nand the factorization could not be completed. - aborting",info);
			free(AAT);
			return -1;
		}
		
		// initial y values
		
		// create rhs for the solve
		cblas_dgemv(CblasRowMajor, CblasNoTrans, m, n, 1.0, A, n, c, 1, 0.0, v[4], 1); // A * c
		
		dpotrs_(&uplo, &m, &one, AAT, &m, v[4], &m, &info); //  complete sovlve: AAT^-1 * A * c = y
		
		if(info < 0) {
			printf("\nError in dpotrs - argument %i has illegal value - aborting",-info);
			free(AAT);
			return -1;
		}
		
		// initial x values (primary variables)
		
		// copy b into a temp array for the solve
		double *btemp = malloc(m * sizeof(double));
		cblas_dcopy(m, b, 1, btemp, 1);
		
		dpotrs_(&uplo, &m, &one, AAT, &m, btemp, &m, &info); //  complete sovlve: AAT^-1 * b = btemp
		
		cblas_dgemv(CblasRowMajor, CblasTrans, m, n, 1.0, A, n, btemp, 1, 0.0, v[0], 1); // AT * btemp
		
		free(btemp);
		
		// initial z values (dual variables)
		cblas_dgemv(CblasRowMajor, CblasTrans, m, n, -1.0, A, n, v[4], 1, 0.0, v[1], 1); // AT * y
		
		for(i=0;i<n;i++)
		{
			v[1][i] += c[i]; // c - AT * y 
		}
	}
	
	free(AAT);
	
	// find min in x and z
	double minx = BIG;
	double minz = BIG;
	for(i=0;i<n;i++)
	{
		minx = (v[0][i] < minx) ? v[0][i] : minx;
		minz = (v[1][i] < minz) ? v[1][i] : minz;
	}
	minx *= -1.5;
	minz *= -1.5;
	
	double delx = (minx > 0.0) ? minx : 0.0;
	double delz = (minz > 0.0) ? minz : 0.0;
	 
	double *xtemp = malloc(n * sizeof(double));
	double *ztemp = malloc(n * sizeof(double));
	
	// xtemp = x + delx
	// ztemp = z + delz

	for(i=0;i<n;i++)
	{
		xtemp[i] = v[0][i] + delx;
		ztemp[i] = v[1][i] + delz;
	}
	
	// compute dot product xtemp * ztemp
	dtemp = cblas_ddot(n, xtemp, 1, ztemp, 1);
	
	// compute sum of xtemp and sum ztemp
	double xsum = 0.0;
	double zsum = 0.0;
	
	for(i=0;i<n;i++)
	{
		xsum += xtemp[i];
		zsum += ztemp[i];
	}
	
	// update delx and delz
	delx += 0.5*(dtemp/xsum);
	delz += 0.5*(dtemp/zsum);
	
	// update initial x and z values
	for(i=0;i<n;i++)
	{
		v[0][i] += delx;
		v[1][i] += delz;
	}
	
	// repeat above for s & w ??
	if(nu > 0)
	{
		// first set initial s & w values from starting x & z values
		cblas_dgemv(CblasRowMajor, CblasTrans, m, n, 1.0, A, n, v[4], 1, 0.0, xtemp, 1); // AT * y = xtemp
		for(i=0;i<nu;i++)
		{
			v[2][i] = (u[i] - v[0][i]); // s = u - x
			v[3][i] = (xtemp[i] + v[1][i] - c[i]); // w = AT*y + z - c 
		}
		
		// reuse varible for x (=s) & z (=w)

		// find min in s and w
		minx = BIG;
		minz = BIG;
		for(i=0;i<nu;i++)
		{
			minx = (v[2][i] < minx) ? v[2][i] : minx;
			minz = (v[3][i] < minz) ? v[3][i] : minz;
		}
		minx *= -1.5;
		minz *= -1.5;
		
		delx = (minx > 0.0) ? minx : 0.0;
		delz = (minz > 0.0) ? minz : 0.0;
		
		// xtemp = s + delx
		// ztemp = z + delz
		for(i=0;i<nu;i++)
		{
			xtemp[i] = v[2][i] + delx;
			ztemp[i] = v[3][i] + delz;
		}
		
		// compute dot product xtemp * ztemp
		dtemp = cblas_ddot(nu, xtemp, 1, ztemp, 1);
		
		// compute sum of xtemp and sum ztemp
		xsum = 0.0;
		zsum = 0.0;
		
		for(i=0;i<nu;i++)
		{
			xsum += xtemp[i];
			zsum += ztemp[i];
		}
		
		// update delx and delz
		delx += 0.5*(dtemp/xsum);
		delz += 0.5*(dtemp/zsum);
		
		// update initial s and w values
		for(i=0;i<nu;i++)
		{
			v[2][i] += delx;
			v[3][i] += delz;
		}
	}
	
	free(xtemp);
	free(ztemp);
	
	return 0;
}

// function to find an initial point
int initialize2(int n,  int m, int nu, double **v, double *A, double *b, double *u, double *c)
{
	int i;
	double dtemp,fact1,fact2;
	double sfact = 1.0;
	
	double sumb = 0.0; // l1 norm of b
	for(i=0;i<m;i++)
	{
		sumb += fabs(b[i]);
	}
	sumb /= sfact; // norm / sfact
	
	fact2=0.0;
	for(i=0;i<n;i++)
	{
		fact2 += fabs(c[i]);
	}
	fact2 += 1.0;
	
	// first compute AA^T matrix and initial y and x values
	double *AAT = malloc(m*m*sizeof(double));
	
	if(m==1)
	{
		// AAT is simply the dot product
		dtemp = cblas_ddot(n, A, 1, A, 1);
		
		// initial x values (primary variables)
		dtemp = *b / dtemp;
		for(i=0;i<n;i++)
		{
			v[0][i] = A[i] * dtemp;
		}
		
		fact1 = -v[0][i];
	}
	else
	{
		char uplo = 'U';
		int one = 1;
		int info;
		// otherwise perform matrix matrix multiplication
		cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, m, m, n, 1.0, A, n, A, n, 0.0, AAT, m);
		
		// factor the matrix ready for inversion
		dpotrf_(&uplo, &m, AAT, &m, &info);
		
		if(info < 0) {
			printf("\nError in dpotrf - argument %i has illegal value - aborting",-info);
			free(AAT);
			return -1;
		}
		else if (info > 0) {
			printf("\nError in dpotrf - the leading minor of order %i is not positive definite,\nand the factorization could not be completed. - aborting",info);
			free(AAT);
			return -1;
		}
		
		// initial x values (primary variables)
		
		// copy b into a temp array for the solve
		double *btemp = malloc(m * sizeof(double));
		cblas_dcopy(m, b, 1, btemp, 1);
		
		dpotrs_(&uplo, &m, &one, AAT, &m, btemp, &m, &info); //  complete sovlve: AAT^-1 * b = btemp
		
		cblas_dgemv(CblasRowMajor, CblasTrans, m, n, 1.0, A, n, btemp, 1, 0.0, v[0], 1); // AT * btemp
		
		free(btemp);

		fact1 = BIG;
		for(i=0;i<n;i++)
		{
			fact1 = (v[0][i] < fact1) ? v[0][i] : fact1; // find min
		}
		fact1 *= -1.0; // change sign
	}
	
	fact1 = (fact1 > sfact) ? fact1 : sfact;
	fact1 = (fact1 > sumb) ? fact1 : sumb;
	
	free(AAT);
	
	// update initial x
	for(i=0;i<n;i++)
	{
		v[0][i] = (v[0][i] > fact1) ? v[0][i] : fact1;
	}
	
	// set initial s
	for (i=0;i<nu;i++) {
		dtemp = u[i] - v[0][i];
		v[2][i] = (dtemp > fact1) ? dtemp : fact1;
	}
	
	// initial y values (zero)
	for(i=0;i<m;i++)
	{
		v[4][i] = 0.0;
	}
	
	// set intial z & w
	for(i=0;i<n;i++)
	{
		if(c[i] >= 0.0)
		{
			v[1][i] = c[i] + fact2;
			if(i<nu) {
				v[3][i] = fact2;
			}
		}
		
		else if(c[i] < -fact2)
		{
			v[1][i] = -c[i];
			if(i<nu) {
				v[3][i] = -2.0*c[i];
			}
		}
		
		else
		{
			v[1][i] = fact2;
			if(i<nu) {
				v[3][i] = fact2 - c[i];
			}
		}
	}
	
	return 0;
}

// function to compute delP(v) - predictor
int predictor(int n,  int m, int nu, double **v, double **delP, double *A, double *b, double *u, double *c)
{
	int i,j,ind;
	double dtemp;
	
	// compute D matrix
	// as D is a diagonal matrix, store as an array
	double *D = malloc(n*sizeof(double));
	
	/*for(i=0;i<n;i++)
	{
		if(v[0][i] == 0.0) { printf("\nx=0"); }
		if(i<nu && v[2][i] == 0.0) { printf("\ns=0"); }
	}*/
	
	for(i=0;i<n;i++)
	{
		dtemp = v[1][i] / v[0][i]; // z/x
		
		if(i<nu) { dtemp += v[3][i] / v[2][i]; } // w/s
		D[i] = 1.0 / dtemp;
	}
	
	// -------------- compute rhs arrays (rb, ru, rc, rxz, rsw) -------------- //
	
	// compute rb
	double *rb = malloc(m*sizeof(double));
	// multiply: rb = A * x
	cblas_dgemv(CblasRowMajor, CblasNoTrans, m, n, 1.0, A, n, v[0], 1, 0.0, rb, 1);
	// compute rb = rb - b
	for(i=0;i<m;i++){ rb[i]-=b[i]; }
	
	// compute ru
	double *ru;
	if(nu>0)
	{
		ru = malloc(nu*sizeof(double));
		for(i=0;i<nu;i++){ ru[i] = v[0][i] + v[2][i] -  u[i]; } // x + s - u
	}
	
	// compute rc
	double *rc = malloc(n*sizeof(double));
	// multiply: rc = A^T * y
	cblas_dgemv(CblasRowMajor, CblasTrans, m, n, 1.0, A, n, v[4], 1, 0.0, rc, 1);
	// compute rc = rc + z - app(w) - c
	for(i=0;i<n;i++)
	{
		rc[i] += v[1][i] - c[i];
		if(i<nu){ rc[i] -= v[3][i]; }
	}
	
	// compute rxz
	double *rxz = malloc(n*sizeof(double));
	for(i=0;i<n;i++)
	{
		rxz[i] = v[0][i] * v[1][i]; // x * z
	}
	
	// compute rsw
	double *rsw;
	if(nu>0)
	{
		rsw = malloc(nu*sizeof(double));
		for(i=0;i<nu;i++)
		{
			rsw[i] = v[2][i] * v[3][i]; // s * w
		}
	}
	
	// update rc
	for(i=0;i<n;i++)
	{
		//dtemp = divZero(rxz[i], v[0][i]); // rxz / x
		dtemp = rxz[i] / v[0][i];
		if(i<nu){ dtemp -= (rsw[i] - v[3][i]*ru[i]) / v[2][i]; }
		rc[i] -= dtemp;
	}
	
	// -------------- compute the update steps (delP) -------------- //
	
	// compute dely (delP[4])
	// compute the ADA^T matrix (assume dense)
	double *ADA = malloc(m*m*sizeof(double));
	double *row_temp = malloc(n*sizeof(double));
	for(i=0;i<m;i++) // row
	{
		ind = n*i; // start of row in A
		for(j=0;j<n;j++)
		{
			row_temp[j] = A[ind+j] * D[j]; // this will be row i in AD
		}
			
		// only compute upper triangle (symmetric matrix)
		ind = m*i; // start of row in ADA
		for(j=i;j<m;j++) // column
		{
			ADA[ind+j] = cblas_ddot(n, row_temp, 1, &A[n*j], 1); // &A[n*j] is the start of row j
		}
	}
	// compute the rhs for the dely solve
	double *y_rhs = malloc(m*sizeof(double));
	// first compute row_temp = D*rc
	for(i=0;i<n;i++)
	{
		row_temp[i] = rc[i] * D[i]; // this will be row i in AD
	}
	// then multiply -A * row_temp = y_rhs
	cblas_dgemv(CblasRowMajor, CblasNoTrans, m, n, -1.0, A, n, row_temp, 1, 0.0, y_rhs, 1);
	// finally minus rb from y_rhs
	for(i=0;i<m;i++)
	{
		y_rhs[i] -= rb[i];
	}
	
	// solve for dely = ADA^-1 *y_rhs (special considerations!)
	if(m==1)
	{
		double fact;
        // error!
		if(ADA[0] < 0.0){
            free(D);
            free(ADA);
            free(rb);
            free(rc);
            free(rxz);
            free(row_temp);
            free(y_rhs);
            if(nu>0) {
                free(ru);
                free(rsw);
            }
            return -1;
        }
		if(ADA[0] < SMALL){fact=BIG;}	// analogous to the Cholesky-infinity factorization
		else {fact=1.0/ADA[0];}			// pseudo factorization of a 1x1 matrix
		delP[4][0] = y_rhs[0] * fact;
	}
	
	else if(m==2)
	{
		double det;
		det = ADA[0]*ADA[3] - ADA[1]*ADA[1]; // ADA is symmetric: ADA[2] = ADA[1]
		det = divZero(1.0, det); // 1/det
		delP[4][0] = det*(y_rhs[0]*ADA[3] - y_rhs[1]*ADA[1]);
		delP[4][1] = det*(y_rhs[1]*ADA[0] - y_rhs[0]*ADA[1]);
	}
	
	else if(m==3)
	{
		double Minv[9];
		
		Minv[0] = ADA[4]*ADA[8] - ADA[5]*ADA[7];
		//Minv[3] = ADA[5]*ADA[6] - ADA[3]*ADA[8];
		//Minv[6] = ADA[3]*ADA[7] - ADA[4]*ADA[6];
		Minv[1] = ADA[2]*ADA[7] - ADA[1]*ADA[8];
		Minv[4] = ADA[0]*ADA[8] - ADA[2]*ADA[6];
		//Minv[7] = ADA[6]*ADA[1] - ADA[0]*ADA[7]; // symmetry
		Minv[2] = ADA[1]*ADA[5] - ADA[2]*ADA[4];
		Minv[5] = ADA[2]*ADA[3] - ADA[0]*ADA[5];
		Minv[8] = ADA[0]*ADA[4] - ADA[1]*ADA[3];
		
		double det = ADA[0]*Minv[0] + ADA[1]*Minv[1] + ADA[2]*Minv[2];
		det = divZero(1.0, det); // 1/det
		
		delP[4][0] = det*(y_rhs[0]*Minv[0] + y_rhs[1]*Minv[1] + y_rhs[2]*Minv[2]);
		delP[4][1] = det*(y_rhs[0]*Minv[1] + y_rhs[1]*Minv[4] + y_rhs[2]*Minv[5]);
		delP[4][2] = det*(y_rhs[0]*Minv[2] + y_rhs[1]*Minv[5] + y_rhs[2]*Minv[8]);
	}

	// more than 3 constraints - use MA57
	else
	{
		int count;
		int len = ((m+1)*m)/2; // storage space
		
		// convert ADA to triplet form
		int *irn = malloc(len*sizeof(int));
		int *jcn = malloc(len*sizeof(int));
		double *At = malloc(len*sizeof(double));
		ind=0;
		for(i=0;i<m;i++) // row
		{
			count = i*m; // start of row
			for(j=i;j<m;j++) // col
			{
				irn[ind] = i+1;
				jcn[ind] = j+1;
				At[ind++] = ADA[count+j];
			}
		}
		
		// solve using MA57
		solve(m, len, irn, jcn, At, 1, y_rhs, 0);
		
		// copy solution
		cblas_dcopy(m, y_rhs, 1, delP[4], 1);
		
		// free memory
		free(irn);
		free(jcn);
		free(At);
	}
	
	// compute delx (delP[0])
	// frist  compute A^T dely = delx
	cblas_dgemv(CblasRowMajor, CblasTrans, m, n, 1.0, A, n, delP[4], 1, 0.0, delP[0], 1);
	// then add rc and times by D
	for(i=0;i<n;i++)
	{
		delP[0][i] += rc[i];
		delP[0][i] *= D[i];
	}
	
	// compute delz (delP[1])
	// frist  compute Z delx + rxz = delz
	// then divide by -X
	for(i=0;i<n;i++)
	{
		dtemp = v[1][i]*delP[0][i] + rxz[i];
		delP[1][i] = -dtemp / v[0][i];
	}
	
	// compute dels (delP[2])
	// dels = -(delx + ru)
	for(i=0;i<nu;i++)
	{
		delP[2][i] = -delP[0][i] - ru[i];
	}
	
	// compute delw (delP[3])
	// delw = -(W*dels + rsw) / S
	for(i=0;i<nu;i++)
	{
		delP[3][i] = v[3][i]*delP[2][i] + rsw[i];
		delP[3][i] = -delP[3][i] / v[2][i];
	}
	
	// clean up and return
	free(D);
	free(ADA);
	free(rb);
	free(rc);
	free(rxz);
	free(row_temp);
	free(y_rhs);
	if(nu>0) {
		free(ru);
		free(rsw);
	}
	
	return 0;
}

// function to compute delC(v) - corrector
void corrector(int n, int nu, double mu, double **v, double **delP)
{
	int i;
	// compute rxz
	double *rxz = malloc(n*sizeof(double));
	for(i=0;i<n;i++)
	{
		rxz[i] = (delP[0][i] * delP[1][i]) - mu; // delx * delz - mu
		rxz[i] = (rxz[i] > 0.0) ? rxz[i] : 0.0; // ensure positivity
	}
	
	// compute rsw
	double *rsw;
	if(nu>0)
	{
		rsw = malloc(nu*sizeof(double));
		for(i=0;i<nu;i++)
		{
			rsw[i] = (delP[2][i] * delP[3][i]) - mu; // dels * delw - mu
			rsw[i] = (rsw[i] > 0.0) ? rsw[i] : 0.0; // ensure positivity
		}
	}
	
	// update delz (delP[1])
	// frist  compute Z delx + rxz = delz
	for(i=0;i<n;i++){ delP[1][i] -= rxz[i] / v[0][i]; }
	
	// update delw (delP[3])
	// delw = -(W*dels + rsw) / S
	for(i=0;i<nu;i++){ delP[3][i] -= rsw[i] / v[2][i]; }
	
	free(rxz);
	if(nu>0) {
		free(rsw);
	}
}

// function to compute the centering parameter
double centering(int n, int m, int nu, double lim, double **v, double **delP)
{	
	double ap = 1.0;
	double ad = 1.0;
	double *z = malloc(n*sizeof(double));
	double *x = malloc(n*sizeof(double));
	double *s, *w;
	if(nu > 0)
	{
		s = malloc(nu*sizeof(double));
		w = malloc(nu*sizeof(double));
	}
	
	int i;
	double dtemp;
	for(i=0;i<n;i++)
	{
		// x update
		dtemp = v[0][i] + delP[0][i];
		if(dtemp < 0.0)
		{
			dtemp = divZero(-v[0][i], delP[0][i]);
			ap = (dtemp < ap) ? dtemp : ap;
		}
		
		// z update
		dtemp = v[1][i] + delP[1][i];
		if(dtemp < 0.0)
		{
			dtemp = divZero(-v[1][i], delP[1][i]);
			ad = (dtemp < ad) ? dtemp : ad;
		}
	}
	for(i=0;i<nu;i++)
	{
		// s update
		dtemp = v[2][i] + delP[2][i];
		if(dtemp < 0.0)
		{
			dtemp = divZero(-v[2][i], delP[2][i]);
			ap = (dtemp < ap) ? dtemp : ap;
		}
		
		// w update
		dtemp = v[3][i] + delP[3][i];
		if(dtemp < 0.0)
		{
			dtemp = divZero(-v[3][i], delP[3][i]);
			ad = (dtemp < ad) ? dtemp : ad;
		}
	}
	// forcast update
	for(i=0;i<n;i++)
	{
		x[i] = v[0][i] + ap * delP[0][i];
		z[i] = v[1][i] + ad * delP[1][i];
	}
	for(i=0;i<nu;i++)
	{
		s[i] = v[2][i] + ap * delP[2][i];
		w[i] = v[3][i] + ad * delP[3][i];
	}
	
	// compute sigma (then round to upper limit)
	double g = dual_gap(n, nu, v);
	
	double sigma;
	sigma = cblas_ddot(n, x, 1, z, 1); // xTz
	if(nu > 0) {
		sigma += cblas_ddot(nu, s, 1, w, 1); // sTw
	}
	sigma /= g;
	sigma *= sigma;
	sigma = (sigma < lim) ? sigma : lim;
	
	g /= (double)(n+nu); // normalize the duality gap
	g *= sigma;
	
	free(x);
	free(z);
	if(nu>0) {
		free(s);
		free(w);
	}

	return g;
}

// function to find the duality gap
double dual_gap(int n, int nu, double **v)
{
	double g =  cblas_ddot(n, v[0], 1, v[1], 1); // xTz
	if(nu > 0) {
		g += cblas_ddot(nu, v[2], 1, v[3], 1); // sTw
	}
	
	return g;
}

// update the variables
void update(int n, int nu, int m, double **v, double **delP)
{
	double g,ming;
	
	// set ap & ad based on non-negativity condition
	double ap = 1.0;
	double ad = 1.0;
	
	int i;
	double dtemp;
	for(i=0;i<n;i++)
	{
		// x update
		dtemp = v[0][i] + delP[0][i];
		if(dtemp < 0.0)
		{
			dtemp = -v[0][i] / delP[0][i];
			ap = (dtemp < ap) ? dtemp : ap;
		}
		
		// z update
		dtemp = v[1][i] + delP[1][i];
		if(dtemp < 0.0)
		{
			dtemp = -v[1][i] / delP[1][i];
			ad = (dtemp < ad) ? dtemp : ad;
		}
	}
	for(i=0;i<nu;i++)
	{
		// s update
		dtemp = v[2][i] + delP[2][i];
		if(dtemp < 0.0)
		{
			dtemp = -v[2][i] / delP[2][i];
			ap = (dtemp < ap) ? dtemp : ap;
		}
		
		// w update
		dtemp = v[3][i] + delP[3][i];
		if(dtemp < 0.0)
		{
			dtemp = -v[3][i] / delP[3][i];
			ad = (dtemp < ad) ? dtemp : ad;
		}
	}
	
	// update v
	cblas_daxpy(n, ap, delP[0], 1, v[0], 1); // x update
	cblas_daxpy(n, ad, delP[1], 1, v[1], 1); // z update
	cblas_daxpy(nu, ap, delP[2], 1, v[2], 1); // s update
	cblas_daxpy(nu, ad, delP[3], 1, v[3], 1); // w update
	cblas_daxpy(m, ad, delP[4], 1, v[4], 1); // y update
	
	// then loop and check duality gap condition
	// reduce ap & ad as required - can change v on the fly
	int stop;
	int count = 0;
	double pstep, dstep;
	double reduce = 0.9995;
	double red2 = 1.0-reduce;
	do
	{
		stop = 1;
		// check on the duality gap condition
		ming = BIG;
		for(i=0;i<n;i++)
		{
			dtemp = v[0][i] * v[1][i]; // xz
			ming = (dtemp < ming) ? dtemp : ming;
		}
		for(i=0;i<nu;i++)
		{
			dtemp = v[2][i] * v[3][i]; // sw
			ming = (dtemp < ming) ? dtemp : ming;
		}
		
		g = dual_gap(n,nu,v) / (double)(n+nu);
		g *= 1.0E-5; // multiply by small number

		if(ming < g) // if condition not met
		{
			pstep = -red2 * ap; dstep = -red2 * ad;
			
			// back track the update
			cblas_daxpy(n, pstep, delP[0], 1, v[0], 1); // x update
			cblas_daxpy(n, dstep, delP[1], 1, v[1], 1); // z update
			cblas_daxpy(nu, pstep, delP[2], 1, v[2], 1); // s update
			cblas_daxpy(nu, dstep, delP[3], 1, v[3], 1); // w update
			cblas_daxpy(m, dstep, delP[4], 1, v[4], 1); // y update
			
			ap *= reduce;
			ad *= reduce; // reduce update step lengths
			stop = 0;
		}
		count++;
		/*if(count > 10) {
			reduce = 0.9;
		}
		else if(count > 5) {
			reduce = 0.95;
		}*/
		red2 = 1.0 - reduce;
		
	} while(stop == 0 && count < 10);
	//printf("\nupdate in %i itt : %12.4e , %12.4e",count,ming,g);
}

// compute the stopping criterion
double stopping(int n,  int m, int nu, double **v, double *A, double lenb, double *b,
			 double lenu, double *u, double lenc, double *c)
{
	int i;
	// compute size of rb, rc, ru
	
	// compute rb
	double *rb = malloc(m*sizeof(double));
	// multiply: rb = A * x
	cblas_dgemv(CblasRowMajor, CblasNoTrans, m, n, 1.0, A, n, v[0], 1, 0.0, rb, 1);
	// compute rb = rb - b
	for(i=0;i<m;i++)
	{
		rb[i]-=b[i];
	}
	
	double stop = cblas_dnrm2(m, rb, 1) / lenb;
	
	// compute rc
	double *rc = malloc(n*sizeof(double));
	// multiply: rc = A^T * y
	cblas_dgemv(CblasRowMajor, CblasTrans, m, n, 1.0, A, n, v[4], 1, 0.0, rc, 1);
	// compute rc = rc + z - app(w) - c
	for(i=0;i<n;i++)
	{
		rc[i] += v[1][i] - c[i];
		if(i<nu)
		{
			rc[i] -= v[3][i];
		}
	}
	
	stop += cblas_dnrm2(n, rc, 1) / lenc;
	
	// compute ru
	double *ru;
	if(nu>0)
	{
		ru = malloc(nu*sizeof(double));
		for(i=0;i<nu;i++)
		{
			ru[i] = v[0][i] + v[2][i] -  u[i]; // x + s - u
		}
		
		stop += cblas_dnrm2(nu, ru, 1) / lenu;
	}
	
	// compute final part of stopping criterion
	double cx = cblas_ddot(n, c, 1, v[0], 1); // c^Tx
	double by = cblas_ddot(m, b, 1, v[4], 1); // b^Ty
	double uw = 0.0;
	if(nu>0) {
		uw = cblas_ddot(nu, u, 1, v[3], 1); // b^Ty
	}
	double dtemp = fabs(by - uw);
	double denom = (fabs(cx) > dtemp) ? fabs(cx) : dtemp;
	denom = (1.0 > denom) ? 1.0 : denom;
	
	stop += fabs(cx - by + uw) / denom;
	
	// free memory
	free(rb);
	free(rc);
	if(nu>0) {
		free(ru);
	}
	
	return stop;
}

// function to print solution to screen
void report(int n,  int m, int nu, double **v)
{
	int i;
	int lim = (nu < 10) ? nu : 10;
	int lim2 = (n < 10) ? n : 10;
	
	printf("\nPrinting solution ...");
	printf("\n\nX values\n--------");
	for(i=0;i<lim2;i++)
	{
		printf("\nx[%i] = %12.6e",i,v[0][i]);
	}
	
	printf("\n\nZ values\n--------");
	for(i=0;i<lim2;i++)
	{
		printf("\nz[%i] = %12.6e",i,v[1][i]);
	}
	
	if(nu > 0)
	{
		printf("\n\nS values\n--------");
		for(i=0;i<lim;i++)
		{
			printf("\ns[%i] = %12.6e",i,v[2][i]);
		}
		
		printf("\n\nW values\n--------");
		for(i=0;i<lim;i++)
		{
			printf("\nw[%i] = %12.6e",i,v[3][i]);
		}
	}
	
	printf("\n\nY values\n--------");
	for(i=0;i<m;i++)
	{
		printf("\ny[%i] = %12.6e",i,v[4][i]);
	}
	
	printf("\n\n-------------------");
}
