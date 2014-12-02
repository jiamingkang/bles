/*
	main.cpp

 	Created on: Dec 01, 2014
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

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "CExMinimiseCompliance.h"

// argv[1] = input file

int main(int argc, char *argv[])
{
	printf("\n\n------*---Start of BLES Version 5.4 Program---*------\n");

	int i,j,i2,j2,k,temp,temp2;	// incrementors
	double ftemp;

	/*----------------------------------------------------------------------------------/
	/																					/
	/		Read input file.															/
	/		Define: mesh, loads & BCs, level set, problem, controls, element matrices	/
	/																					/
	/----------------------------------------------------------------------------------*/

	char filename[100];
	for(i=0;i<100;i++)
	{
		if(argv[1][i] == '.') {break;} // chop off file extension (if found)
		else if(argv[1][i] == '\0'){break;}
		else {filename[i] = argv[1][i];}
	}
	filename[i] = '\0';

	CExMinimiseCompliance cMC;
	cMC.Solve(argv[1], filename);

}
