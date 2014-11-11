/*
 *  ABFG.h
 *
 *  Created by Peter Dunning
 *	Version 5.4 @ 23/04/2014
 */
 
 #include "ls_types.h"

// Function to find the nearest grid node number to a set of co-ordinates
int closeNode(mesh *inMesh, double xp, double yp);

// Function to extract the structure from the lsf
// 1. determines the node and element status
// 2. discretizes the boundary
// 3. computes area ratios for AFG method
void Find_struct(mesh *inMesh, levSet *levelset, boundary *bound_in, double *alpha, double aMin);

// compute area ratio for all elements
void AFG_area(mesh *inMesh, double *alpha, short *NodeStat, short *ElemStat, Coord *AuxNodes,
			  int NumBound, Bseg *Boundary, double aMin);

// compute a cut element area
double InArea(int eStat, int eNum, int *Lnodes, short *NodeStat, int NumNodes, Coord *NodeCoord,
			  Coord *AuxNodes, int NumBound, Bseg *Boundary);

// Function that caculates the area of any Polygon
double PolyArea(int N, Coord *point);

// Function to determine if two lines (p1 -> p2 & p3 -> p4) cross
short LineCross(Coord *pts);

// compute determinant of a 2x2 matrix
double det2(double a, double b, double c, double d);

// Assembles global stiffness (& maybe mass) matrix for AFG method in triplet format (for MA57 solver)
void AFG_Matrix(int mass, double **KE, double **ME, sp_mat *Kg, sp_mat *Mg, sp_mat *lump_mass, double *alpha, 
				mesh *inMesh, isoMat *inMat, double aMin, double mMin);

// function to compute gauss point coords
void Gauss_Coord(mesh *inMesh, Coord *gCoord);

// calculate sensitivies using least squares of integration points for AFG method
void AFG_Sens(mesh *inMesh, boundary *bound_in, double *alpha, isoMat *inMat,  double *Nsens,
			  double **prim, double **dual, int numDual, int numCase, double *wgt, Coord *gCoord,
				double aMin, int mode, double *fact, bool sw, Coord *acc);
