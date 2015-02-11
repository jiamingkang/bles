/*
	CHakBoundary.h

	Created on: Nov 24, 2014
	Author: Peter Dunning, JeeHang Lee

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

#ifndef CHakBoundary_H_
#define CHakBoundary_H_

//
// includes
//

#include "CommonTypes.h"

class CHakBoundary {
public:
	CHakBoundary();
	virtual ~CHakBoundary();

public:
	// function to weight boundary segments so that lengths are more even
	void BsegWgt(CHakBoundary *bound_in, CHakMesh *inMesh);

	// Function to perfrom boundary intergration of objective and constraint shape sens
	void BoundInt(CHakMesh *inMesh, CHakLevelSet *levelset, CHakBoundary *bound_in, int numFunc, double **Nsens,
				  int *Lbound_nums, int *numLbound,  double *Lbound);

private:
	// number of auxillary nodes
	int m_numAux;	//	NumAux

	// aux node coords
	Coord *m_pAuxNodes;  //*AuxNodes

	// number of boundary segments
	int m_numBound; //NumBound
	
	// boundary segment data
	Bseg *Bound; //*Bound

	// grid node -> aux node connectivity
	int *m_pConnAuxNode; //*na_conn
	
	// indicator array for na_conn (compressed storage)
	int *m_pIndConnAuxNode; //*na_conn_ind

	// boundary segment lengths 
	double *m_pLenBseg; //*BsegLen
	
	// boundary segment weights
	double *m_pWeightBseg; //*Bwgt
};

#endif /* CHakBoundary_H_ */
