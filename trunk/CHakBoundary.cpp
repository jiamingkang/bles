/*
	CBoundary.cpp

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

#include <stdio.h>
#include <stdlib.h>
#include <math.h>


#include "CBoundary.h"
#include "CSolver.h"

//
// Constructor / Destructor
//

CBoundary::CBoundary() {
	// TODO Auto-generated constructor stub

}

CBoundary::~CBoundary() {
	// TODO Auto-generated destructor stub
}

// function to weight boundary segments so that lengths are more even
void CBoundary::BsegWgt(boundary *bound_in, mesh *inMesh)
{
	// read in data
	int NumNodes = inMesh->NumNodes;
	Coord *NodeCoord = inMesh->NodeCoord;
	int NumBound = bound_in->NumBound;
	Bseg *Lbound = bound_in->Bound;
	Coord *AuxNodes = bound_in->AuxNodes;

	// memory allocation
    free(bound_in->BsegLen);
    bound_in->BsegLen = (double *) malloc(NumBound*sizeof(double)); // boundary segment lengths
    free(bound_in->Bwgt);
	bound_in->Bwgt = (double *) malloc(NumBound*sizeof(double));	// boundary segment weights

	int i,j,temp,n1,n2;
	int s1,s2;
	double la,lb;
	double x1,x2,y11,y2,x,y;

	// first work out lengths of the boundary segments
	for(i=0;i<NumBound;i++)
	{
		// read in end node numbers
		n1 = Lbound[i].n1;
		n2 = Lbound[i].n2;

		// read in node co-ordinates
		if(n1<NumNodes) {
			x1 = NodeCoord[n1].x;
			y11 = NodeCoord[n1].y;
		}
		else {
			temp = n1 - NumNodes;
			x1 = AuxNodes[temp].x;
			y11 = AuxNodes[temp].y;
		}

		if(n2<NumNodes) {
			x2 = NodeCoord[n2].x;
			y2 = NodeCoord[n2].y;
		}
		else {
			temp = n2 - NumNodes;
			x2 = AuxNodes[temp].x;
			y2 = AuxNodes[temp].y;
		}

		x = x1-x2; y = y11-y2;
		x*=x; y*=y;
		bound_in->BsegLen[i] = sqrt(x+y); // boundary segment length
	}

	for(i=0;i<NumBound;i++)
	{
		// find the length of the two neighbouring boundary segments
		s1=0; s2=0;
		for(j=0;j<NumBound;j++)
		{
			if(i != j)
			{
				if(Lbound[i].n1 == Lbound[j].n1 || Lbound[i].n1 == Lbound[j].n2)
				{
					la = bound_in->BsegLen[j];
					s1=1;
				}
				if(Lbound[i].n2 == Lbound[j].n1 || Lbound[i].n2 == Lbound[j].n2)
				{
					lb = bound_in->BsegLen[j];
					s2=1;
				}
			}
			if(s1+s2 == 2)
			{
				j=NumBound; // break out of for loop
			}
		}

		bound_in->Bwgt[i] = lb / (la+lb); // weight for length associated with node 1
		//bound_in->Bwgt[i] = 0.5; // for even weighting
	}
}

// Function to perfrom boundary intergration of objective and constraint shape sens
void CBoundary::BoundInt(mesh *inMesh, levSet *levelset, boundary *bound_in, int numFunc, double **Nsens,
				int *Lbound_nums, int *numLbound,  double *Lbound)
{
	// read in data
	bool *fixed = levelset->fixed;
	double *lsf = levelset->lsf;
	int NumNodes = inMesh->NumNodes;
	double tol = inMesh->tol;

	int NumBound = bound_in->NumBound;
	Bseg *bptr = bound_in->Bound;
	double *Bseglen = bound_in->BsegLen;
	double *wgt = bound_in->Bwgt;
	int Ntot = NumNodes + bound_in->NumAux;

	int i,j,count; // incrementors
	int n1,n2,p1,p2; // boundary point node numbers
	double ftemp,len,len2,sens1,sens2; // temporary variables
	double *Ltemp, *b_len; // temporary array to store boundary data
	Ltemp = (double *) calloc(Ntot*numFunc,sizeof(double));
	b_len = (double *) calloc(Ntot, sizeof(double));

	for(i=0;i<NumBound;i++)
	{
		// read in end node numbers
		n1 = bptr[i].n1;
		n2 = bptr[i].n2;

		len = Bseglen[i] * wgt[i];
		len2 = Bseglen[i] * (1.0 - wgt[i]);

		b_len[n1] += len;
		b_len[n2] += len2; // boundary length integration

		for(j=0;j<numFunc;j++)
		{
			p1 = Ntot*j + n1;
			p2 = Ntot*j + n2; // point to correct locations in Ltemp

			sens1 = Nsens[j][n1]; sens2 = Nsens[j][n2];
			ftemp = (sens2 - sens1)* wgt[i] + sens1; // interpolate

			Ltemp[p1] += 0.5 * (ftemp + sens1) * len;
			Ltemp[p2] += 0.5 * (ftemp + sens2) * len2; // Update shape sens intergral
		}
	}

	// define bound array (boundary point numbers)
	count = 0;
	for(i=0;i<Ntot;i++)
	{
		// if associated boundary length is > than small value then include
		if(b_len[i] > tol)
		{
			if(i >= NumNodes || (!fixed || !fixed[i]) ) // only check fixed for grid nodes
			{ Lbound_nums[count++] = i; }
		}
		// if node designated a boundary node, then include anyway
		else if( (i < NumNodes) && (!fixed || !fixed[i]) && (fabs(lsf[i])<tol) )
		{
			Lbound_nums[count++] = i;
		}
	}

	*numLbound = count; // store length of Lbound array as first entry of array

	for(i=0;i<count;i++)
	{
		for(j=0;j<numFunc;j++)
		{
			p1 = j*Ntot + Lbound_nums[i]; // point to correct place in Ltemp
			p2 = j*count + i; // point to correct place in Lbound
			Lbound[p2] = Ltemp[p1];
		}
	}

	free(Ltemp);
	free(b_len);
}
