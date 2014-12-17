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
	void BsegWgt(boundary *bound_in, mesh *inMesh);

	// Function to perfrom boundary intergration of objective and constraint shape sens
	void BoundInt(mesh *inMesh, levSet *levelset, boundary *bound_in, int numFunc, double **Nsens,
				  int *Lbound_nums, int *numLbound,  double *Lbound);

};

#endif /* CHakBoundary_H_ */
