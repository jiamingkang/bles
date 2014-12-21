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
#include "CHakMesh.h"
//#include "CHakIsoMaterial.h"
#include <string>

#define MAX_CHARS_PER_LINE 512
#define MAX_TOKENS_PER_LINE 16
#define MAX_HOLES 500

//
// Class CHakInput 
//	- Mainly read the input file and retrieve initial data
//	- All initial data are preserved in this scope. No manipulation is allowed.
//	- To use the data, external entities should get reference of the member.
//
class CHakInput {
public:
	CHakInput();
	CHakInput(char *filename) { m_filename = filename; }
	virtual ~CHakInput();

public:
	// function to read new-style input file & assign to classes
	int ReadFile(char *datafile);

	// function to read new-style input file
	int read_input(char *datafile, mesh *inMesh, int *numMat, isoMat *inMat, levSet *levelset, prob *lsprob,
				   ctrl *control, int **fixDof, int *numCase, double **load, int *freeDof, sp_mat *lump_mass, bool *sw, Coord **acc);

protected:
	static int icmpfunc(const void *p1, const void *p2);

private:
	// Mesh 
	CHakMesh m_mesh;

	// isotropic material including material data and its count
	//CHakIsoMaterial m_isoMat;

private:
	// input file name
	std::string m_filename;
};

#endif /* CINPUT_H_ */
