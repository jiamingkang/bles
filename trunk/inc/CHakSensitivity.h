/*
	CHakSensitivity.h

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

#ifndef CHakSensitivity_H_
#define CHakSensitivity_H_

// include
#include "CommonTypes.h"
#include "CHakMaterial.h"

// pre-processors
#define GP 1.577350269189626 // 1 + 1/sqrt(3)
#define GN 0.422649730810374 // 1 - 1/sqrt(3)

// non-dimensional stran-disp matrices at each gauss point
// only valid for Q4 elements
const static double BMat[4][24] = { {-GP, 0.0, GP, 0.0, GN, 0.0, -GN, 0.0,
	0.0, -GP, 0.0, -GN, 0.0, GN, 0.0, GP,
	-GP, -GP, -GN, GP,  GN, GN, GP, -GN},

	{-GP, 0.0, GP, 0.0, GN, 0.0, -GN, 0.0,
		0.0, -GN, 0.0, -GP, 0.0, GP, 0.0, GN,
		-GN, -GP, -GP, GP,  GP, GN, GN, -GN},

	{-GN, 0.0, GN, 0.0, GP, 0.0, -GP, 0.0,
		0.0, -GN, 0.0, -GP, 0.0, GP, 0.0, GN,
		-GN, -GN, -GP, GN,  GP, GP, GN, -GP},

	{-GN, 0.0, GN, 0.0, GP, 0.0, -GP, 0.0,
		0.0, -GP, 0.0, -GN, 0.0, GN, 0.0, GP,
		-GP, -GN, -GN, GN,  GN, GP, GP, -GP} };

// interpolation matrix (each row is a guass point)
const static double Q4_inter[4][4] = {
	{GP*GP*0.25, GP*GN*0.25, GN*GN*0.25, GP*GN*0.25},
	{GP*GN*0.25, GP*GP*0.25, GP*GN*0.25, GN*GN*0.25},
	{GN*GN*0.25, GP*GN*0.25, GP*GP*0.25, GP*GN*0.25},
	{GP*GN*0.25, GN*GN*0.25, GP*GN*0.25, GP*GP*0.25} };

//
// class CHakSensitivity
// 		Description: This class is originated from Sens.h and ABFG.h
//
class CHakSensitivity {

public:
	CHakSensitivity();
	virtual ~CHakSensitivity();

public:
	// calculate sensitivies using least squares of integration points for AFG method
	void AFG_Sens(mesh *inMesh, boundary *bound_in, double *alpha, isoMat *inMat,  double *Nsens,
				  double **prim, double **dual, int numDual, int numCase, double *wgt, Coord *gCoord,
					double aMin, int mode, double *fact, bool sw, Coord *acc);

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

public:
	CMaterial GetMaterial() { return m_cMaterial; }
	void SetMaterial(CMaterial cMat) { m_cMaterial = cMat; }

private:
	CMaterial m_cMaterial;

};

#endif /* CHakSensitivity_H_ */
