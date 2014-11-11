/*
 *  Sens.h
 *
 *  Created by Peter Dunning
 *	Version 5.4 @ 23/04/2014
 */
 
#include "ls_types.h"

#define gp 1.577350269189626 // 1 + 1/sqrt(3)
#define gn 0.422649730810374 // 1 - 1/sqrt(3)

// non-dimensional stran-disp matrices at each gauss point
// only valid for Q4 elements
const static double BMat[4][24] = { {-gp, 0.0, gp, 0.0, gn, 0.0, -gn, 0.0,
	0.0, -gp, 0.0, -gn, 0.0, gn, 0.0, gp,
	-gp, -gp, -gn, gp,  gn, gn, gp, -gn},
	
	{-gp, 0.0, gp, 0.0, gn, 0.0, -gn, 0.0,
		0.0, -gn, 0.0, -gp, 0.0, gp, 0.0, gn,
		-gn, -gp, -gp, gp,  gp, gn, gn, -gn},
	
	{-gn, 0.0, gn, 0.0, gp, 0.0, -gp, 0.0,
		0.0, -gn, 0.0, -gp, 0.0, gp, 0.0, gn,
		-gn, -gn, -gp, gn,  gp, gp, gn, -gp},
	
	{-gn, 0.0, gn, 0.0, gp, 0.0, -gp, 0.0,
		0.0, -gp, 0.0, -gn, 0.0, gn, 0.0, gp,
		-gp, -gn, -gn, gn,  gn, gp, gp, -gp} };

// interpolation matrix (each row is a guass point)
const static double Q4_inter[4][4] = {
	{gp*gp*0.25, gp*gn*0.25, gn*gn*0.25, gp*gn*0.25},
	{gp*gn*0.25, gp*gp*0.25, gp*gn*0.25, gn*gn*0.25},
	{gn*gn*0.25, gp*gn*0.25, gp*gp*0.25, gp*gn*0.25},
	{gp*gn*0.25, gn*gn*0.25, gp*gn*0.25, gp*gp*0.25} };

// Function that calculates the sensitivity of a node by a least squares (2nd order) filter of near-by gauss points
int Lsens(Coord *pt, int xMax, int xMin, int yMax, int yMin, double aMin, double *alpha, double r2,
				  Coord *gCoord, double *gSens, Elem **Number, int wFlag, int numDual, double *out);
					
// Function to calculate sensitivity values for an element at 4 gauss points
// for for a plane 4-node element (Q4)
void GaSens_Q4(int *tnodes, double **prim, double **dual, double alpha, double h, double t,
			isoMat *mat, int Gcount, double *gSens, int numCase, int numDual, double *wgt, bool sw, Coord *acc);

// Function to calculate additional sensitivity part for eigenvalues
// for for a plane 4-node element (Q4)
void GaEigSens_Q4(int *tnodes, double **prim, double **dual, double alpha, double h, double t,
				  isoMat *inMat, int Gcount, double *gSens, int num_eig, double *eig);

// function that computes bar senstivites for compliance (possible multi-load case)
void barSens(mesh *inMesh, double *bar_sens, double **prim, int numCase, double *wgt);

// function to compute designable bc sensitvities for compliance (possible multi-load case)
void bcSens(mesh *inMesh, double *bc_sens, double **prim, int numCase, double *wgt);

// function to compute designable material design varibles for compliance
void matSens_comp(mesh *inMesh, isoMat *inMat, double *KE, double *mat_sens, int numCase, double *wgt,
                  double **prim, double **dual, double *alpha, double aMin, bool sw, Coord *acc);

// function to compute designable material design varibles for eigenvalues
void matSens_eig(mesh *inMesh, isoMat *inMat, double *KE, double *ME, double *mat_sens, 
				 int numEig, double *eig_vals, double **eig_vecs, double *alpha, double aMin);

// function to compute designable material design H-S varibles for eigenvalues
void HS_Sens_eig(mesh *inMesh, isoMat *inMat, double *KE, double *ME, double *mat_sens,
                 int numEig, double *eig_vals, double **eig_vecs, double *alpha, double aMin);

// derivative of E for H-S bound material model
double dE_dalpha(double alpha, double hs_int, isoMat *mat1, isoMat *mat2);
