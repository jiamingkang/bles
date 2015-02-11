/*
   CHakMesh.h

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

#ifndef CHakMesh_H_
#define CHakMesh_H_

#include "CommonTypes.h"
#include "CHakFiniteElement.h"

class CHakMesh {
public:
	CHakMesh();
	virtual ~CHakMesh();

public:
	// Function that numbers all elements and nodes in the FG domain
	//void Numbering(mesh *inMesh);

	// Function to find the nearest grid node number to a set of co-ordinates
	//int closeNode(mesh *inMesh, double xp, double yp);

	// Node Co-ordinate calculation function
	//void Coordinates(mesh *inMesh);

	// function that orders node numbers into a 2D based on their relative positions
	//void NodeNums2(mesh *inMesh);

	// function to number bars elements
	//void Bar_numbering(mesh *inMesh);

//
// OOD version
//
public:

	// Function to extract the structure from the lsf
	// 1. determines the node and element status
	// 2. discretizes the boundary
	// 3. computes area ratios for AFG method
	void FindStruct(CHakLevelSet *m_levelset, CHakBoundary *m_boundary, double *alpha, double aMin);

	// function to compute area ratio for all elements
	void ComputeAreaRatio(double *alpha, short *NodeStat, short *ElemStat, Coord *AuxNodes, int NumBound, Bseg *Boundary, double aMin);

	// Function that numbers all elements and nodes in the FG domain
	void Numbering();

	// Node Co-ordinate calculation function
	void Coordinates();

	// function that orders node numbers into a 2D based on their relative positions
	void NodeNums2D();

	// function to number bars elements
	void NumberingBarElements();

	// Function to find the nearest grid node number to a set of co-ordinates
	int CloseNode(double xp, double yp);

public:
//private:
	// num elements in x & y
	int m_elemX, m_elemY; //elemX, elemY;

	// total number of elements
	int m_numElem;  //NumElem

	// total number of nodes
	int m_numNodes; //NumNodes

	// element edge length & thickness
	double m_lenEdge, m_thinkness; //originally h,t

	// rounding tolerance
	double m_tolerance;; //originlly tol

	// max domain dimensions
	double m_maxX, m_maxY;	 //original maxX and maxY

	// pointer to node coordinate array
	Coord *m_pNodeCoord; //originally NodeCoord

	// pointer to element node numbering (2d array)
	Elem **m_pNumber; // originally Number

	// number of nodes in X & Y
	int m_nodeX, m_nodeY;  //Original NodeX, NodeY

	// structured node numbering (2d array)
	int **m_pNodes2D;  //Original Nodes2

	// indicates which material each element is made from
	int *m_pMaterialType; // Orignal mat_type

//private:
	//
	// additional variables for bar reinforcements
	// --> may be a separated class (e.g. design variable for bar-reinforcement) or structure
	//
	//

	// flag to indicate if bar elements are in the mesh
	bool m_bars;  //Original bars

	// number of bars
	int m_numBars;  //Original NumBars

	// array to hold bar numbers (n1 = dof 1, n2 = dof 2, e = elem no.)
	Bseg *m_pBarNums; //original bar_nums

	// array of bar areas
	double *m_pBarAreas; //original bar_areas

	// max and min bar areas
	double m_barMax, m_barMin; //original bar_max, bar_min

	// bar material type
	CHakMaterial *m_pBarMaterialType; //original bar_mat

//private:
	//
	// additional variables for designable boundary conditions
	// --> may be a separated class (e.g. design variable for boundary conditions)
	//

	// flag to indicate designable bcs are present
	bool m_bDesignableBc; //  des_bc

	// number of elements with designable bcs
	int m_numBc; //NumBC

	// array to store elements with designable bcs
	int *m_pBcNums;  //BC_nums

	// array to store designable bc variables
	double *m_K_bc; //K_bc

	// maximum stiffness of bc springs
	double m_K0_bc;  //K0_bc

//private:
	//
	// additional variables for elemental material design variables
	// 	--> may be a separated class (e.g. design variable),
	//

	// flag to indicate designable material is present
	bool m_bDesignableMaterial; //des_mat

	// pointers the two materials
	int m_materialOne, m_materialTwo; //mat1, mat2

	// number of elements wit hdesignable material
	int m_numDesignableMaterialt; //NumDesMat

	// array of element numbers with designable material
	int *m_pMaterialElems; //mat_elems

	// array to store designable material variables
	double *m_pMaterialVars; //mat_vars

	// flag for simultaneous (true) or sequential (false) optimization of material
	bool m_bDmSim; //dm_sim

	// linear (true) or H-S bound (false) material model
	bool m_bMaterialLin; //mat_lin
};

#endif /* CHakMesh_H_ */


