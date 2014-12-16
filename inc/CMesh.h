/*
   CMesh.h

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

#ifndef CMESH_H_
#define CMESH_H_

#include "CommonTypes.h"
#include "CFiniteElement.h"

class CMesh {
public:
	CMesh();
	virtual ~CMesh();

public:
	// Function that numbers all elements and nodes in the FG domain
	void Numbering(mesh *inMesh);

	// Function to find the nearest grid node number to a set of co-ordinates
	int closeNode(mesh *inMesh, double xp, double yp);

	// Function to extract the structure from the lsf
	// 1. determines the node and element status
	// 2. discretizes the boundary
	// 3. computes area ratios for AFG method
	void Find_struct(mesh *inMesh, levSet *levelset, boundary *bound_in, double *alpha, double aMin);

	// function to compute area ratio for all elements
	void AFG_area(mesh *inMesh, double *alpha, short *NodeStat, short *ElemStat, Coord *AuxNodes,
					int NumBound, Bseg *Boundary, double aMin);

	// Node Co-ordinate calculation function
	void Coordinates(mesh *inMesh);

	// function that orders node numbers into a 2D based on their relative positions
	void NodeNums2(mesh *inMesh);

	// function to number bars elements
	void Bar_numbering(mesh *inMesh);

private:
	// num elements in x & y
	int m_elemX, m_elemY;

	// total number of elements
	int m_numElem;

	// total number of nodes
	int m_numNodes;

	// element edge length & thickness
	double m_lenEdge, m_thinkness;

	// rounding tolerance
	double m_tolerance;;
	double m_maxX, m_maxY;	// max domain dimensions

	// pointer to node coordinate array
	Coord *m_pNodeCoord;

	// pointer to element node numbering (2d array)
	Elem **m_pNumber;

	// number of nodes in X & Y
	int NodeX, NodeY;

	// structured node numbering (2d array)
	int **Nodes2;

	// indicates which material each element is made from
	int *mat_type;

private:
	//
	// additional variables for bar reinforcements
	// --> may be a separated class (e.g. design variable for bar-reinforcement)
	//

	// flag to indicate if bar elements are in the mesh
	bool bars;

	// number of bars
	int NumBars;

	// array to hold bar numbers (n1 = dof 1, n2 = dof 2, e = elem no.)
	Bseg *bar_nums;

	// array of bar areas
	double *bar_areas;

	// max and min bar areas
	double bar_max, bar_min;

	// bar material type
	isoMat *bar_mat;

private:
	//
	// additional variables for designable boundary conditions
	// --> may be a separated class (e.g. design variable for boundary conditions), jeehanglee@gmail.com
	//

	// flag to indicate designable bcs are present
	bool des_bc;

	// number of elements with designable bcs
	int NumBC;

	// array to store elements with designable bcs
	int *BC_nums;

	// array to store designable bc variables
	double *K_bc;

	// maximum stiffness of bc springs
	double K0_bc;

private:
	//
	// additional variables for elemental material design variables
	// 	--> may be a separated class (e.g. design variable), jeehanglee@gmail.com
	//

	// flag to indicate designable material is present
	bool des_mat;

	// pointers the two materials
	int mat1, mat2;

	// number of elements wit hdesignable material
	int NumDesMat;

	// array of element numbers with designable material
	int *mat_elems;

	// array to store designable material variables
	double *mat_vars;

	// flag for simultaneous (true) or sequential (false) optimization of material
	bool dm_sim;

	// linear (true) or H-S bound (false) material model
	bool mat_lin;
};

#endif /* CMESH_H_ */


