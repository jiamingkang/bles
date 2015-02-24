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
	void initialLsf(CHakMesh *inMesh, int NumHole, CirH *holes, int NumRect, Coord *Rect, double lBand);

	// Function to calcualte lsf around a rectangular hole
	void RectHole(int NumNodes, Coord *NodeCoord, int NumRect, Coord *Rect, double *lsf);

	// Function to calculate active and mine nodes for the narrow band method
	void NarBand(CHakMesh *inMesh, double lBand);

	// function to weight boundary segments so that lengths are more even
	void BsegWgt(CHakBoundary *bound_in, CHakMesh *inMesh);

	// Function to perfrom boundary intergration of objective and constraint shape sens
	void BoundInt(CHakMesh *inMesh, CHakBoundary *bound_in, int numFunc, double **Nsens,
					int *Lbound_nums, int *numLbound,  double *Lbound);

	// Function to calculate extension velocities using the fast marching method
	void Vext(CHakMesh *inMesh, CHakBoundary *bound_in, double *Vnorm);

	// function that works out Velocity for nodes close to the boundary
	void LocalVext2(int NodeX, int NodeY, int **Nodes, double *Vnorm, double *Vtemp, double *lsf, double *lsf_temp, double tol,
					bool *known, bool *trial, CHakBoundary *bound_in, double h, int NumNodes, Coord *NodeCoord, int sign);

	// Function to calcualte gradient using WENO scheme
	double GradWENO(int Xi, int Yj, int num, double *lsf, CHakMesh *inMesh, int sign, double Vn, double dt);

	// sub-function for GradWENO
	double GWsub(double v1,double v2,double v3,double v4,double v5);

	// Function to re-initalise the lsf as a signed distance function - similar to Vext above
	void ReInt(CHakMesh *inMesh);

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
	int SLPsubSol4(CHakMesh *inMesh, CHakBoundary *bound_in, double alt, double *delCon, int numCon,
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

	// ---- set of functions for hole insertion algorithm ---- //
public:
	//Get the index and map of elements and nodes involved in hole insertion
	int HoleMap(mesh *inMesh, int *h_index, int *h_EmapX, int *h_EmapY, levSet *levelset, double *alpha);

	//Intialise hole level set function
	void Intialise_hole_lsf(mesh *inMesh, int h_count, double holeCFL, int *h_index, int *h_EmapX, int *h_EmapY,
			double *h_Nsens, double *h_lsf, double *h_area, double *h_gMin, double *h_gMax, int num_gam);

	//Get hole lsf from current gammer and sensitvities
	void Get_h_lsf(int NumNodes, int num_gam, int *h_index, double *Gam, double *h_Nsens, double *h_lsf);

	//Get area of region after hole removal
	double Get_h_area(int h_count, int *h_EmapX, int *h_EmapY, Elem **Number, double h, double *h_area, double *h_lsf);

	//Get area of cut hole element
	double Get_Hole_Area_Elem(double lsf1, double lsf2, double lsf3, double lsf4, int eNum);

	//Finite difference for intial gammer set up
	double Hole_FD_SG(int NumNodes, int num_gam, int *h_index, double *Gam, double *h_Nsens,
			double *h_lsf, int h_count, int *h_EmapX, int *h_EmapY, Elem **Number, double h,
			double *h_area, double DA, double TargetV, double perb);

	//Function to perform finite difference anaylsis to get the sensivity of each objective function to a change in gammer.
	void Hole_FD(int NumNodes, int num_gam, int num_sens, int *h_index, double *Gam,
			double *h_Nsens, double *h_Esens, double *h_lsf, int h_count, int *h_EmapX,
			int *h_EmapY, Elem **Number, int NumElem, double h, double *h_area,
			double *h_gMin, double *GamS);


	// function to set velocity (or move dist) using SLP - filter method with jole insertion
	int SLPsubSol4_hole(mesh *inMesh, levSet *levelset, boundary *bound_in, double alt, double *delCon, int numCon,
	                    double **sens, double *cA, int n, int *bound, double *lam_in, int *active, double *Vnorm,
						double *pred, int numAdd, double *add_sens, double *add_min, double *add_max,
						double *add_change, int pinfo, int *h_index, double *h_Nsens, double *h_Esens,
						double *h_lsf, int h_count, int *h_EmapX, int *h_EmapY, double *h_area, double *h_gMin,
						double *h_gMax, int *Reint);

//  temp code. Refactoring required.
public:
	CHakSolver getSolver() { return m_solver; }
	void setSolver(CHakSolver solver) { m_solver = solver; }

private:
	CHakSolver m_solver;

public:
//private:
	// number of nodes for level set discretization
	int m_numNode;	 //num

	// pointer to array containing nodal lsf values
	double *m_pNodalLsf;  //*lsf
	
	//bool *bound; // nodes that are on fixed boundaries (not currently used)
	
	// array for fixed lsf values (NULL if none fixed)
	bool *m_pFixedLsf; //*fixed

	// array for nodes within narrow band (and not fixed)
	bool *m_pActive; //*active

	// number of mines on edge of narrow band
	int m_numMine; //numMine

	// mine node numbers
	int *m_pMine;    //*mine
};

#endif /* CHakLevelSet_H_ */
