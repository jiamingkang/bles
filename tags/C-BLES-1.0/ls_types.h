/*
 *  ls_types.h
 *
 *  Created by Peter Dunning
 *	Version 5.4 @ 23/04/2014
 *  Data structures used in BLES V5
 */

#ifndef __LS_TYPES_H
#define __LS_TYPES_H

#include <stdbool.h>
#include <Accelerate/Accelerate.h>
//#include "cblas.h"

#define NUM_DOF 2	// current code for 2d planar elements
#define KE_SIZE (NUM_DOF * NUM_DOF * 16) // number of entries in elem matrix

// Element and node Numbering information structure
typedef struct 
{
	int n;	// element number
	int a,b,c,d; // node numbers
} Elem;

// 2D Co-ordinate inforamtion structure
typedef struct
{
	double x,y; // 2D Co-ordinates
} Coord;
				
// Structure for inital circular hole data
typedef struct
{
	double x,y; // center co-ordinates
	double r;	// hole radius
} CirH;

// Structure to hold data for a boundary segment
typedef struct
{
	int n1, n2; // node numbers for end points
	int e;  // associated element
}	Bseg;

// isotropic material structure
typedef struct
{
	double e,v,rho; // modulus, Poisson's ratio, density
    double k,g;     // bulk and shear modulii
	double mat[9];  // material property matrix
} isoMat;

// mesh information structure
typedef struct
{
	int elemX, elemY;	// num elements in x & y
	int NumElem, NumNodes; // total num elements & nodes
	double h,t;			// element edge length & thickness
	double tol;			// rounding tollerance
	double maxX, maxY;	// max domain dimensions
	Coord *NodeCoord;	// pointer to node coordinat array
	Elem **Number;		// pointer to element node numbering (2d array)
	int NodeX, NodeY;	// number of nodes in X & Y
	int **Nodes2;		// structured node numbering (2d array)
	int *mat_type;		// indicates which material each element is made from
	// additional variables for bar reinforcements
	bool bars;			// flag to indicate if bar elements are in the mesh
	int NumBars;		// number of bars
	Bseg *bar_nums;		// array to hold bar numbers (n1 = dof 1, n2 = dof 2, e = elem no.)
	double *bar_areas;	// array of bar areas
	double bar_max, bar_min; // max and min bar areas
	isoMat *bar_mat;	// bar material type
	// additional variables for designable boundry conditions
	bool des_bc;		// flag to indicate designable bcs are present
	int NumBC;			// number of elements with designable bcs
	int *BC_nums;		// array to store elements with designable bcs
	double *K_bc;		// array to store designable bc variables
	double K0_bc;		// maximum stiffness of bc springs
	// additional variables for elemental material design variables
	bool des_mat;		// flag to indicate designable material is present
	int mat1, mat2;		// pointers the two materials
	int NumDesMat;		// number of elements wit hdesignable material
	int *mat_elems;		// array of element numbers with designable material
	double *mat_vars;	// array to store designable material variables
    bool dm_sim;        // flag for simultaneous (true) or sequential (false) optimization of material
    bool mat_lin;       // linear (true) or H-S bound (false) material model
} mesh;

// boundary discretization info
typedef struct
{
	int NumAux;		// number of auxillary nodes
	Coord *AuxNodes; // aux node coords
	int NumBound; // number of boundary segments
	Bseg *Bound; // boundary segment data
	int *na_conn; // grid node -> aux node connectivity
	int *na_conn_ind; // indicator array for na_conn (compressed storage)
	double *BsegLen, *Bwgt; // boundary segment lengths & weights
} boundary;

// sparse matrix in triplet form
typedef struct
{
	int ne; // number of entries (in arrays)
	int *irn, *jcn; // row and column indicators
	double *A; // matrix entry value
} sp_mat;

// structure containing optimization controls
typedef struct
{
	int maxItt; // max number of iterations (>0)
	int pinfo;  // amount of output (1->3, less->more)
	double gm;  // convergence criterion (gamma)
	double lband; // narrow band width (2h -> large)
	double aMin;  // minimum area ratio (i.e. small or zero) - stiffness
	double mMin;  // minimum area ratio for mass
} ctrl;

// structure to hold constraint info
typedef struct
{
	int type; // constraint type identifier (vol, freq, etc..)
	int sign; // < , > or = (-ve, +ve, 0) constraint
	double data[4]; // data for constraint (first value is value)
} cnst;

// structure to hold optimization problem defintion
typedef struct
{
	int obj; // objective type identifier
	int num; // number of constraints
	cnst *con; // pointer to array of constraint structs
} prob;

// structure for the level set function
typedef struct
{
	int num;	 // number of nodes for level set discretization
	double *lsf; // pointer to array containing nodal lsf values
	//bool *bound; // nodes that are on fixed boundaries (not currently used)
	bool *fixed; // array for fixed lsf values (NULL if none fixed)
	bool *active; // array for nodes within narrow band (and not fixed)
	int numMine;  // number of mines on edge of narrow band
	int *mine;    // mine node numbers
} levSet;
	
#endif
