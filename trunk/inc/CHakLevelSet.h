/*
   CHakLevelSet.h

    Created on: 24 Nov 2014
    Author: Peter Dunning, Khalid Ismail

	-- This file is part of Topology Optimisation Opensource Project,
 	owned by MSO (Multidisciplinary and Structural Optimisation) Research Group
 	(http://people.bath.ac.uk/ens/MSORG/index.html) at University of Bath.

 	The project is led by Dr. Alicia Kim
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

#ifndef CHakLevelSet_H_
#define CHakLevelSet_H_

#include "CommonTypes.h"
#include "CHakSolver.h"

class CHakLevelSet {
public:
	CHakLevelSet();
	virtual ~CHakLevelSet();

public:
	// function to create the initial signed distance function - inc holes
	void initialLsf(mesh *inMesh, levSet *levelset, int NumHole, CirH *holes, int NumRect, Coord *Rect, double lBand);

	// Function to calcualte lsf around a rectangular hole
	void RectHole(int NumNodes, Coord *NodeCoord, int NumRect, Coord *Rect, double *lsf);

	// Function to calculate active and mine nodes for the narrow band method
	void NarBand(mesh *inMesh, levSet *levelset, double lBand);

	// function to weight boundary segments so that lengths are more even
	void BsegWgt(boundary *bound_in, mesh *inMesh);

	// Function to perfrom boundary intergration of objective and constraint shape sens
	void BoundInt(mesh *inMesh, levSet *levelset, boundary *bound_in, int numFunc, double **Nsens,
					int *Lbound_nums, int *numLbound,  double *Lbound);

	// Function to calculate extension velocities using the fast marching method
	void Vext(mesh *inMesh, levSet *levelset, boundary *bound_in, double *Vnorm);

	// function that works out Velocity for nodes close to the boundary
	void LocalVext2(int NodeX, int NodeY, int **Nodes, double *Vnorm, double *Vtemp, double *lsf, double *lsf_temp, double tol,
					bool *known, bool *trial, boundary *bound_in, double h, int NumNodes, Coord *NodeCoord, int sign);

	// Function to calcualte gradient using WENO scheme
	double GradWENO(int Xi, int Yj, int num, double *lsf, mesh *inMesh, int sign, double Vn, double dt);

	// sub-function for GradWENO
	double GWsub(double v1,double v2,double v3,double v4,double v5);

	// Function to re-initalise the lsf as a signed distance function - similar to Vext above
	void ReInt(mesh *inMesh, levSet *levelset);

	// Function to reinitalise lsf for inital set of trial nodes - similar to LocalVext above
	void LocalInt(int NodeX, int NodeY, int **Nodes, double *lsf, double *lsf_temp,
				  bool *known, bool *trial, double h, double tol, int sign);

	// ---- set of functions for SLP Level-set method ---- //

	// obtain a boundary move vector from input of lam, s and move limits
	void get_delD(int n, int m, double *x, double *lam, double *s, double *up_lim, double *low_lim);

	// function to obtain gradients by finite difference
	void get_slpGrad(int n, int m, int numVar, double *lam, double *s, double *c, double *up_lim,
					 double *low_lim, double *max_lam, double *min_lam, double *grad, int pinfo);

	// sub-function for sorting in ascending order
	static int dcmpfunc (const void * p1, const void * p2);

	// function to get limits for lambda
	void getLamLim(int n, int m, double *inmax, double *inmin, double *s, double *c,
						double *up_lim, double *low_lim, int pinfo);

	// function to set velocity (or move dist) using SLP - filter method
	int SLPsubSol4(mesh *inMesh, levSet *levelset, boundary *bound_in, double alt, double *delCon, int numCon,
	                double **sens, double *cA, int n, int *bound, double *lam_in, int *active, double *Vnorm, double *pred,
					int numAdd, double *add_sens, double *add_min, double *add_max, double *add_change, int pinfo);

	// sub-solve function for the trust region method
	int trust_sub(int n, int m, int nu, double *x, double *c, double *A, double *b, double *u, int pinfo);

	// function to minimise constraint violation using Newtons method (more efficient)
	int con_min3(int n, int m, int numCon, double *lam,  double *s, double *A, double *b,
				 double *lam_min, double *lam_max, double *up_lim, double *low_lim, int pinfo);

	// function to find an initial set of lambda values
	void get_lam0(int n, int m, int numCon, double *lam,  double *s, double *cA, double *b,
				  double *lam_min, double *lam_max);

// jeehanglee@gmail.com: temp code. Refactoring required.
public:
	CHakSolver getSolver() { return m_solver; }
	void setSolver(CHakSolver solver) { m_solver = solver; }

private:
	CHakSolver m_solver;
};

#endif /* CHakLevelSet_H_ */
