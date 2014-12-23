/*
	main.cpp

 	Created on: Dec 01, 2014
	Author: Peter Dunning, JeeHang Lee, Khalid Ismail

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
	printf("\n\n------*---Start of BLES Version 5.4 Program---*------\n\n");

	CExMinimiseCompliance cMC;
	int maxIter = 0;
	int iter = 0;
	int convergence = 0;
	    
	// jeehanglee@gmail.com: with Exception handling...
	if (argc > 1)
	{
		// return control.maxIter from initialise
		maxIter = cMC.Initialise(argv[1]);
	}
	else
	{
		printf("Failure in Reading Input File... Application Terminated... \n\n");
		return 0;
	}

	// automatically stop after max iterations
	// otherwise continue the interation(s)
    while (iter < maxIter) 
	{
		// do analysis including FEA...
        cMC.Analyse(iter);
        
		// adjusting sensitivities
        convergence = cMC.Sensitivity();
        
		// terminate the process if converged...
		if (convergence == 1) 
			break;
        
		// do optimisation at this iteration step
        cMC.Optimise();
        
        iter++;
    } 
    
    cMC.Output();
    
    return 0;
}
