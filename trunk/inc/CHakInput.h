/*
	CInput.h

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

#ifndef CINPUT_H_
#define CINPUT_H_

#include "CommonTypes.h"

#define MAX_CHARS_PER_LINE 512
#define MAX_TOKENS_PER_LINE 16
#define MAX_HOLES 500

class CHakInput {
public:
	CHakInput();
	virtual ~CHakInput();

public:
	// function to read new-style input file
	int read_input(char *datafile, mesh *inMesh, int *numMat, isoMat *inMat, levSet *levelset, prob *lsprob,
				   ctrl *control, int **fixDof, int *numCase, double **load, int *freeDof, sp_mat *lump_mass, bool *sw, Coord **acc);

protected:
	static int icmpfunc(const void *p1, const void *p2);

private:

};

#endif /* CINPUT_H_ */
