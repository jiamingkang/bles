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
};

#endif /* CMESH_H_ */


