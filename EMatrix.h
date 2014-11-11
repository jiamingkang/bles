/*
 *  EMatrix.h
 *
 *  Created by Peter Dunning
 *	Version 5.4 @ 23/04/2014
 */
 
#include "ls_types.h"

// For odd (i + j = odd) Matrix enrites	
double Aodd(int i,int j,double e,double g);
// For even (i + j = even) Matrix enrites	
double Aeven(int i,int j,double e,double g);

// Complete square IN element stiffness matix computation
// only for isotropic material
void KEMatrix(double *KE, isoMat *mat, double t);

// consistent mass matrix for a 2D plane element
void MEMatrix(double *ME, double rho, double h, double t);

// Function that assemebles the element matricies into the global matrix (in triplet form)
// Symmetric matrix - upper triangle
void Assemble2(int *K_begin, int *M_begin, int nNod, int nDof, int *tnodes, double *KE, double *ME, 
			   sp_mat *Kg, sp_mat *Mg, int mass);

// function to compute elastic modulus from two materials
double HS_mat(double alpha, double hs_int, isoMat *mat1, isoMat *mat2);

// function to compute self-weight load vector
// function to compute self-weight load vector
void self_weight(mesh *inMesh, isoMat *inMat, double aMin, double mMin, double *alpha,
                 int freeDof, int *dofMap, int numCase, double *load_in, double *load_out, Coord *acc);
