/*
   CHakLevelSet.cpp

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

#include "CHakSolver.h"
#include "CHakLevelSet.h"

CHakLevelSet::CHakLevelSet() {
	// TODO Auto-generated constructor stub

}

CHakLevelSet::~CHakLevelSet() {
	// TODO Auto-generated destructor stub
}

//
// Implementation - Interfaces
//

// function to create the initial signed distance function - inc holes
void CHakLevelSet::initialLsf(mesh *inMesh, levSet *levelset, int NumHole, CirH *holes, int NumRect, Coord *Rect, double lBand)
{
	int i,j;
	double ftemp,dist,minX,minY;

	// read in data
	int NumNodes = inMesh->NumNodes;
	Coord *NodeCoord = inMesh->NodeCoord;
	double maxX = inMesh->maxX;
	double maxY = inMesh->maxY;
	double tol = inMesh->tol;
	double *lsf = levelset->lsf;

	//For all nodes find the closest edge of the outer domain (initial lsf)
	for(i=0;i<NumNodes;i++)
	{
		// Is node closer to the right or left edge
		ftemp = maxX - NodeCoord[i].x;
		minX = ( (NodeCoord[i].x - ftemp) < tol ) ? NodeCoord[i].x : ftemp;

		// Is node closer to the top or bottom edge
		ftemp = maxY - NodeCoord[i].y;
		minY = ( (NodeCoord[i].y - ftemp) < tol )  ? NodeCoord[i].y : ftemp;

		// Signed distance function is then the lower of minX & minY
		lsf[i] = (minX < minY) ? minX : minY;
	}

	// Insert initial circular holes
	// For all circular holes change signed distance fucntion if less that it current value
	for(i=0;i<NumHole;i++)
	{
		for(j=0;j<NumNodes;j++)
		{
			// calculate x & y distance from node to hole center
			minX = holes[i].x - NodeCoord[j].x;
			minY = holes[i].y - NodeCoord[j].y;

			// calcualte distance from hole center to node
			ftemp = minX * minX;
			ftemp += minY * minY;
			dist = sqrt(ftemp);
			dist -= holes[i].r; // signed distance from the hole edge to the node

			// If signed distance to hole edge is less than current value, then update
			if((dist - lsf[j]) < tol) { lsf[j] = dist; }
		}
	}

	//call function to update lsf if a rectangular holes exists
	if(NumRect > 0) {
		RectHole(NumNodes, NodeCoord, NumRect, Rect, lsf);
	}

	//re-initialize signed distance funciton after initial hole insetion
	ReInt(inMesh, levelset);

	// create the inital narrow band region
	NarBand(inMesh, levelset, lBand);
}

// Function to calcualte lsf around a rectangular hole
void CHakLevelSet::RectHole(int NumNodes, Coord *NodeCoord, int NumRect, Coord *Rect, double *lsf)
{
	int i,j,k;
	double xtemp,ytemp,ftemp,dist;
	double minX, minY, maxX, maxY;

	// For all rectangular holes change signed distance function if less that it current value
	for(i=0;i<NumRect;i++)
	{
		// read in rectangular hole dimensions
		k = 2 * i;
		minX = Rect[k].x;
		maxX = Rect[k+1].x;
		minY = Rect[k].y;
		maxY = Rect[k+1].y;

		for(j=0;j<NumNodes;j++)
		{
			// read in node coords
			xtemp = NodeCoord[j].x;
			ytemp = NodeCoord[j].y;
			// if node to left or right of the rectangle
			if(((xtemp > maxX) || (xtemp < minX)) && (ytemp <= maxY) && (ytemp >= minY))
			{
				dist = minX - xtemp;
				ftemp = xtemp - maxX;
				dist = (fabs(dist) < fabs(ftemp)) ? dist : ftemp;
			}
			// else if node above or below the rectangle
			else if(((ytemp > maxY) || (ytemp < minY)) && (xtemp <= maxX) && (xtemp >= minX))
			{
				dist = minY - ytemp;
				ftemp = ytemp - maxY;
				dist = (fabs(dist) < fabs(ftemp)) ? dist : ftemp;
			}
			// If node inside the rectagle, find closest edge
			else if((xtemp <= maxX) && (xtemp >= minX) && (ytemp <= maxY) && (ytemp >= minY))
			{
				dist = minX - xtemp;
				ftemp = minY - ytemp;
				dist = (fabs(dist) < fabs(ftemp)) ? dist : ftemp;
				ftemp = minY - ytemp;
				dist = (fabs(dist) < fabs(ftemp)) ? dist : ftemp;
				ftemp = xtemp - maxX;
				dist = (fabs(dist) < fabs(ftemp)) ? dist : ftemp;
				ftemp = ytemp - maxY;
				dist = (fabs(dist) < fabs(ftemp)) ? dist : ftemp;
			}

			// if node to bottom left of rectangle
			else if ((xtemp < minX) && (ytemp < minY))
			{
				ftemp = (minX - xtemp) * (minX - xtemp);
				ftemp += (minY - ytemp) * (minY - ytemp);
				dist = sqrt(ftemp);
			}
			// if node to bottom right of rectangle
			else if ((xtemp > maxX) && (ytemp < minY))
			{
				ftemp = (maxX - xtemp) * (maxX - xtemp);
				ftemp += (minY - ytemp) * (minY - ytemp);
				dist = sqrt(ftemp);
			}
			// if node to top left of rectangle
			else if ((xtemp < minX) && (ytemp > maxY))
			{
				ftemp = (minX - xtemp) * (minX - xtemp);
				ftemp += (maxY - ytemp) * (maxY - ytemp);
				dist = sqrt(ftemp);
			}
			// if node to top right of rectangle
			else if ((xtemp > maxX) && (ytemp > maxY))
			{
				ftemp = (maxX - xtemp) * (maxX - xtemp);
				ftemp += (maxY - ytemp) * (maxY - ytemp);
				dist = sqrt(ftemp);
			}
			else	{
				dist = 1.0e+20;
				printf("\nERROR! cant locate node %i in RectHole function!",j+1);
			}
			// If signed distance to hole edge is less than current value, then update
			if(dist < lsf[j]) {
				lsf[j] = dist;
			}
		}
	}
}

// Function to calculate active and mine nodes for the narrow band method
void CHakLevelSet::NarBand(mesh *inMesh, levSet *levelset, double lBand)
{
	int NumNodes = inMesh->NumNodes;
	double h = inMesh->h;
	double tol = inMesh->tol;
	double *lsf = levelset->lsf;
	bool *active = levelset->active;
	bool *fixed = levelset->fixed;

	// reset mine array memory
	levelset->mine = (int *) realloc(levelset->mine,NumNodes * sizeof(int));
	int mineCnt = 0; // count number of mine nodes

	int i;
	double ftemp;
	double lMine = lBand - h;	// Mine "band-width"

	for(i=0;i<NumNodes;i++)
	{
		ftemp = fabs(lsf[i]); // absolute value of the signed distance function
		// If node within narrow band then it is active - unless it is fixed!
		if( ((ftemp - lBand) < tol) && (!fixed || !fixed[i]))
		{
			active[i] = true;
			// If node is also outdie lMine then define it as a mine
			if((ftemp - lMine) > -tol)
			{
				levelset->mine[mineCnt++] = i; // add node num to mine array
			}
		}
		else { active[i] = false; }
	}

	// reallocate memory for mine array
    if(mineCnt>0)
    {
        levelset->mine = (int *) realloc(levelset->mine,mineCnt * sizeof(int));
        levelset->numMine = mineCnt;
    }
    else
    {
        levelset->mine = (int *) realloc(levelset->mine,sizeof(int));
        levelset->numMine = 0;
    }
}

// Function to perfrom boundary intergration of objective and constraint shape sens
void CHakLevelSet::BoundInt(mesh *inMesh, levSet *levelset, boundary *bound_in, int numFunc, double **Nsens,
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

// Function to calculate extension velocities using the fast marching method
void CHakLevelSet::Vext(mesh *inMesh, levSet *levelset, boundary *bound_in, double *Vnorm)
{
	// read in data
	int NodeX = inMesh->NodeX - 1;
	int NodeY = inMesh->NodeY - 1;
	int **Nodes2 = inMesh->Nodes2;
	int NumNodes = inMesh->NumNodes;
	double tol = inMesh->tol;
	double *lsf = levelset->lsf;
	double h = inMesh->h;
	bool *active = levelset->active;

	int i,j,sign; // incrementors
	int Nnum,Nnum2,Nnum3; // node numbers
	int Ci,Cj; // variables to store x,y location of node with minimum lsf_temp values in trial set
	double side;
	double lsf_min; // variable to track minimum lsf_temp value in trial set
	double Af, Bf, Cf; // variables for lsf_temp calculation by solution of a quadratic equn.
	double vt1,vt2;	// variables for tempary velocities
	double ftemp, dtemp; // temp variables
	double Lsum, Fsum; //varibales to calcualte Vext from neighbour nodes

	bool *known = (bool *) calloc(NumNodes,sizeof(bool)); // array to for computed velocity values
	bool *trial = (bool *) calloc(NumNodes,sizeof(bool)); // array for trail nodes in fast marching method
	double *lsf_temp = (double *) calloc(NumNodes,sizeof(double)); // temp lsf during fast marching method
	double *Vtemp = (double *) calloc(NumNodes,sizeof(double)); // array to store velocity values

	short flag; // used to detect when all extension velocities for active node set have been calcualted
	double *lsf_trial = (double *) calloc(NumNodes,sizeof(double)); // array to store trial values of the lsf_temp during the fast marching method

	// set two arrays to store nodes to be updated during each iteration of the fast marching method
	int *xind, *yind, ncount;
	xind = (int *) malloc(NumNodes * sizeof(int));
	yind = (int *) malloc(NumNodes * sizeof(int));

	// Initalise Known set
	for(i=0;i<NumNodes;i++)
	{
		if(active[i]) // if active
		{
			// If node is on the boundary, then add to known set
			if( fabs(lsf[i]) < tol )
			{
				known[i] = true;
				lsf_temp[i] = 0.0; //lsf[i]; // this should be boundary values i,e lsf ~= 0.0
				Vtemp[i] = Vnorm[i];
			}
		}
	}

	// run the extension velocity in both directions
	for(sign=-1;sign<2;sign+=2)
	{
		side = (sign > 0) ? 1.0 : -1.0;

		// Initalise Trial set
		for(j=1;j<NodeY;j++)
		{
			for(i=1;i<NodeX;i++)
			{
				Nnum = Nodes2[i][j]; // read in node number
				// If current node not in known set,
				// but is in the narrow band active set and is the correct side of the boundary
				if( !known[Nnum] && (lsf[Nnum]*side > 0.0) && active[Nnum] )
				{
					if(lsf[Nodes2[i-1][j]]*side < 0.0) {
						trial[Nnum] = true;
						}
					else if(lsf[Nodes2[i+1][j]]*side < 0.0) {
						trial[Nnum] = true;
						}
					else if(lsf[Nodes2[i][j-1]]*side < 0.0) {
						trial[Nnum] = true;
						}
					else if(lsf[Nodes2[i][j+1]]*side < 0.0) {
						trial[Nnum] = true;
						}
				}
			}
		}

		// calculate lsf_temp and Vnorm for all current trial functions
		LocalVext2(NodeX, NodeY, Nodes2, Vnorm, Vtemp, lsf, lsf_temp, tol, known, trial, bound_in, h, NumNodes, inMesh->NodeCoord, sign);

		// update trial set by considering neighbours of nodes in known set
		for(j=1;j<NodeY;j++)
		{
			for(i=1;i<NodeX;i++)
			{
				Nnum = Nodes2[i][j]; // read in node number
				// If current node not in known set
				// but is in the narrow band active set and is the correct side of the boundary
				if( !known[Nnum] && active[Nnum] && (lsf[Nnum]*side > 0.0) )
				{
					// If any neighbouring node is in known set then curret node is a trial node
					if(known[Nodes2[i-1][j]]) {
						trial[Nnum] = true;
						}
					else if(known[Nodes2[i+1][j]]) {
						trial[Nnum] = true;
						}
					else if(known[Nodes2[i][j-1]]) {
						trial[Nnum] = true;
						}
					else if(known[Nodes2[i][j+1]]) {
						trial[Nnum] = true;
						}
				}
			}
		}

		// Do until trial set is empty

		// set a maximum lsf value to larger than maximum domain length
		double lsf_max = (NodeX > NodeY) ? (100.0 * h * (double)NodeX) : (100.0 * h * (double)NodeY);

		do // do unitl all requied velocities are calculated
		{
			flag = 0;
			lsf_min = lsf_max; // initalize to a value much bigger than the element edge length

			// For all trial nodes calculate lsf temp
			for(j=1;j<NodeY;j++)
			{
				for(i=1;i<NodeX;i++)
				{
					Nnum = Nodes2[i][j]; // read in node number
					// If current node not in known set,
					// but is in the narrow band active set and is the correct side of the boundary
					if(trial[Nnum])
					{
						flag = 1; // indicates that more velocites required calculating

						if(lsf_trial[Nnum] < tol) // if trial values needs to be calculated
						{
							Af = 0.0; Bf = 0.0; Cf = 0.0;
							ftemp = 0.0;  // Initalise values

							// look at all neighbouring nodes
							Nnum2 = Nodes2[i][j-1];
							Nnum3 = Nodes2[i][j+1]; // read in node number (above and below)

							// If neigbouring node is in known set update co-efficents for lsf_temp calc
							// ensure only closest node to boundary is used for upwind scheme
							if( known[Nnum2] && known[Nnum3] )	{
								Af += 1.0;
								ftemp =  (lsf_temp[Nnum2] < lsf_temp[Nnum3]) ? lsf_temp[Nnum2] : lsf_temp[Nnum3]; // choose closest known node
							}
							else if(known[Nnum2])	{
								Af += 1.0;
								ftemp = lsf_temp[Nnum2];
							}
							else if(known[Nnum3]) {
								Af += 1.0;
								ftemp = lsf_temp[Nnum3];
							}

							Bf += ftemp;
							Cf += ftemp * ftemp; // Update quadratic co-efficients

							Nnum2 = Nodes2[i-1][j];
							Nnum3 = Nodes2[i+1][j]; // read in node number (left and right)
							ftemp = 0.0; // re-initalise for other direction

							// If neigbouring node is in known set update co-efficents for lsf_temp calc
							// ensure only closest node to boundary is used for upwind scheme
							if( known[Nnum2] && known[Nnum3] )	{
								Af += 1.0;
								ftemp =  (lsf_temp[Nnum2] < lsf_temp[Nnum3]) ? lsf_temp[Nnum2] : lsf_temp[Nnum3]; // choose closest known node
							}
							else if(known[Nnum2])	{
								Af += 1.0;
								ftemp = lsf_temp[Nnum2];
							}
							else if(known[Nnum3])	{
								Af += 1.0;
								ftemp = lsf_temp[Nnum3];
							}

							Bf += ftemp;
							Cf += ftemp * ftemp; // Update quadratic co-efficients

							// Need to do final update of the Bf and Cf co-efficients
							Bf *= 2.0;
							Cf -= h * h;

							if(Bf < 0.0) {
								printf("\nERROR! In Vext for node %i, Bf = %12.4e",Nnum+1,Bf);
							}

							// Calculate lsf_temp by solving the quadratic equation
							ftemp = Bf * Bf;
							ftemp -= 4 * Af * Cf;
							dtemp = sqrt(ftemp);
							ftemp = Bf + dtemp;
							lsf_trial[Nnum] = ftemp / (2.0 * Af);
							if(Af < tol){ printf("\nError in Vext Af = 0!"); }

							if( (lsf_trial[Nnum] < tol) ) {
								printf("\nERROR! In Vext for node %i, lsf_trial = %lf",Nnum+1,lsf_trial[Nnum]);
								printf("\nAf=%lf,Bf=%lf,Cf=%lf",Af,Bf,Cf);
							}
						}

						// update minimum lsf_value
						lsf_min = (lsf_trial[Nnum] < lsf_min) ? lsf_trial[Nnum] : lsf_min;
					}
				}
			}

			if(flag == 0) {	break; } // if there are no trial nodes, then end the loop

			ncount = 0; // re-initialize
			lsf_min *= 1.001; // increase slightly
			// check to see which nodes need updating this iteration
			for(j=1;j<NodeY;j++)
			{
				for(i=1;i<NodeX;i++)
				{
					Nnum = Nodes2[i][j]; // read in node number

					if(lsf_trial[Nnum] < 0.0)	{
						printf("\nERROR IN ReInt: lsf_trial = %lf, for node %i",lsf_trial[Nnum],Nnum);
					}

					// if a trial node has a trail lsf value <= minimum then update
					if( trial[Nnum] && lsf_trial[Nnum] < lsf_min )
					{
						// For trial nodes with lowest lsf temp value
						Nnum = Nodes2[i][j]; // read in node number
						lsf_temp[Nnum] = lsf_trial[Nnum]; // update lsf_temp for node with trial value
						// look at all neighbouring nodes
						Lsum = 0.0;
						Fsum = 0.0; // initaize sum variables to zero

						// look at all neighbouring nodes
						Nnum2 = Nodes2[i][j-1];
						Nnum3 = Nodes2[i][j+1]; // read in node number (above and below)

						// If neigbouring node is in known set - potentially use for Vext calc
						// ensure only closest node to boundary is used for upwind scheme
						if( known[Nnum2] && known[Nnum3] )
						{
							Af = fabs(lsf_temp[Nnum] - lsf_temp[Nnum2]); // absolute diff to node below
							Bf = fabs(lsf_temp[Nnum] - lsf_temp[Nnum3]); // absolute diff to node above
							Cf = lsf_temp[Nnum3] - lsf_temp[Nnum2];		 // diff between above and below nodes

							// choose values in shortest direction
							if(fabs(Cf) < tol)  {	// if distance is same both sides
								Lsum += Af; // could be Bf, doesn't matter
								vt1 = Vtemp[Nnum2];
								vt2 = Vtemp[Nnum3]; // find absolute values of neighbouring velocities
								// choose velocity that will move lsf towards zero, or least distance away from zero
								Fsum += ((sign * vt1) > (sign * vt2)) ? (Af * vt1) : (Bf * vt2);
							}
							else if(Cf > -tol) {	// else if Af < Bf, i.e. node below is closer
								Lsum += Af;
								Fsum += Af * Vtemp[Nnum2];
							}
							else if(Cf < tol)  {	// else if Bf < Af, i.e. node above is closer
								Lsum += Bf;
								Fsum += Bf * Vtemp[Nnum3];
							}
							else	{
								printf("\nERROR! Can't calculate Vext for Node %i",Nnum);
							}
						}
						else if(known[Nnum2])
						{
							Af = fabs(lsf_temp[Nnum] - lsf_temp[Nnum2]);
							Lsum += Af;
							Fsum += Af * Vtemp[Nnum2];
						}
						else if(known[Nnum3])
						{
							Af = fabs(lsf_temp[Nnum] - lsf_temp[Nnum3]);
							Lsum += Af;
							Fsum += Af * Vtemp[Nnum3];
						}

						Nnum2 = Nodes2[i-1][j];
						Nnum3 = Nodes2[i+1][j]; // read in node number (left and right)

						// If neigbouring node is in known set- potentially use for Fext calc
						// ensure only closest node to boundary is used for upwind scheme
						if( known[Nnum2] && known[Nnum3] )
						{
							Af = fabs(lsf_temp[Nnum] - lsf_temp[Nnum2]); // absolute diff to left node
							Bf = fabs(lsf_temp[Nnum] - lsf_temp[Nnum3]); // absolute diff to right node
							Cf = lsf_temp[Nnum3] - lsf_temp[Nnum2];		 // diff between left and right nodes

							// choose values in shortest direction
							if(fabs(Cf) < tol) {	// if distance is same both sides
								Lsum += Af; // could be Bf, doesn't matter
								vt1 = Vtemp[Nnum2];
								vt2 = Vtemp[Nnum3]; // find absolute values of neighbouring velocities
								// choose velocity that will move lsf towards zero, or least distance away from zero
								Fsum += ((sign * vt1) > (sign * vt2)) ? (Af * vt1) : (Bf * vt2);
							}
							else if(Cf > -tol) {	// else if Af < Bf, i.e. left node is closer
								Lsum += Af;
								Fsum += Af * Vtemp[Nnum2];
							}
							else if(Cf < tol) {	// else if Bf < Af, i.e. right node is closer
								Lsum += Bf;
								Fsum += Bf * Vtemp[Nnum3];
							}
							else	{
								printf("\nERROR! Can't calculate Vext for Node %i",Nnum);
							}
						}
						else if(known[Nnum2])
						{
							Af = fabs(lsf_temp[Nnum] - lsf_temp[Nnum2]);
							Lsum += Af;
							Fsum += Af * Vtemp[Nnum2];
						}
						else if(known[Nnum3])
						{
							Af = fabs(lsf_temp[Nnum] - lsf_temp[Nnum3]);
							Lsum += Af;
							Fsum += Af * Vtemp[Nnum3];
						}

						// Calculate Extension velocity and move node to known set
						Vtemp[Nnum] = Fsum / Lsum;
						xind[ncount] = i;
						yind[ncount++] = j; // store update node for next bit
					}
				}
			}

			// update known set for next iteration
			for(i=0;i<ncount;i++)
			{
				Ci = xind[i]; Cj = yind[i]; // read in updated node indicators
				Nnum = Nodes2[Ci][Cj];	// read in node number
				trial[Nnum] = false;
				known[Nnum] = true; // update known & trial sets
			}

			// update trial set for next iteration
			for(i=0;i<ncount;i++)
			{
				Ci = xind[i]; Cj = yind[i]; // read in updated node indicators

				// look at all neighbouring nodes
				Nnum2 = Nodes2[Ci][Cj-1];
				Nnum3 = Nodes2[Ci][Cj+1]; // read in node number (above and below)

				// Check to see if node below needs to be added to trial set
				if(!known[Nnum2] && active[Nnum2] && lsf[Nnum2]*side > 0.0)
				{
					trial[Nnum2] = true;
					lsf_trial[Nnum2] = 0.0; // reset trial value
				}

				// Check to see if node above needs to be added to trial set
				if(!known[Nnum3] && active[Nnum3] && lsf[Nnum3]*side > 0.0)
				{
					trial[Nnum3] = true;
					lsf_trial[Nnum3] = 0.0; // reset trial value
				}

				Nnum2 = Nodes2[Ci-1][Cj];
				Nnum3 = Nodes2[Ci+1][Cj]; // read in node number (left and right)

				// Check to see if node to left needs to be added to trial set
				if(!known[Nnum2] && active[Nnum2] && lsf[Nnum2]*side > 0.0)
				{
					trial[Nnum2] = true;
					lsf_trial[Nnum2] = 0.0; // reset trial value
				}

				// Check to see if node to right needs to be added to trial set
				if(!known[Nnum3] && active[Nnum3] && lsf[Nnum3]*side > 0.0)
				{
					trial[Nnum3] = true;
					lsf_trial[Nnum3] = 0.0; // reset trial value
				}
			}

			if(ncount == 0)	{
				printf("\nERROR! ncount=0 in Vext! Aborting");
			}
		}
		while(ncount > 0);
	}

	// finally copy calcualted Vtemp values to Vnorm array
	cblas_dcopy(NumNodes, Vtemp, 1, Vnorm, 1);

	// Free memory used
	free(known);
	free(trial);
	free(lsf_temp);
	free(Vtemp);
	free(lsf_trial);
	free(xind);
	free(yind);
}

// function that works out Velocity for nodes close to the boundary
void CHakLevelSet::LocalVext2(int NodeX, int NodeY, int **Nodes, double *Vnorm, double *Vtemp, double *lsf, double *lsf_temp, double tol,
				bool *known, bool *trial, boundary *bound_in, double h, int NumNodes, Coord *NodeCoord, int sign)
{
	// read data
	Coord *AuxNodes = bound_in->AuxNodes;
	int *na_conn = bound_in->na_conn;
	int *na_conn_ind = bound_in->na_conn_ind;

	double s,s1,s2,t,t1,t2; // distance variables
	double vs,vs1,vs2,vt,vt1,vt2; // velocity variables
	int i,j,k,k2,stop; // incrementors
	double ftemp; // temp varibale
	int Nnum,Nnum2; // node number variables
	int check; // variable to check aux node has been correctly identified

	// for all nodes
	for(j=1;j<NodeY;j++)
	{
		for(i=1;i<NodeX;i++)
		{
			Nnum = Nodes[i][j];
			if(trial[Nnum])
			{
				stop = na_conn_ind[Nnum+1];
				// look at each neighbouring node in turn

				// -----------look at node below-----------
				Nnum2 = Nodes[i][j-1];
				if(fabs(lsf[Nnum2]) < tol) // If neigbouring node is on the boundary
				{
					t1 = h; // distance to boundary = element edge length
					vt1 = Vnorm[Nnum2]; // store velocity of the node
				}
				else if((lsf[Nnum] * lsf[Nnum2]) < 0.0) // otherwise check to see is grid line is intersected
				{
					check = -1;
					for(k2=na_conn_ind[Nnum];k2<stop;k2++) // For all auxillary nodes (connected to Nnum node)
					{
						k = na_conn[k2]; // aux node num
						if(fabs(NodeCoord[Nnum].x - AuxNodes[k].x) < tol)	// If aux node has same x-coord
						{
							// If aux node lies between the two grid nodes
							if( (NodeCoord[Nnum].y > AuxNodes[k].y) && (AuxNodes[k].y > NodeCoord[Nnum2].y) )
							{
								t1 = fabs(NodeCoord[Nnum].y - AuxNodes[k].y);
								vt1 = Vnorm[NumNodes+k]; // store velocity of the auxillary node
								check = 1;
								break;
							}
						}
					}

					if(check == -1)	{
						printf("\nERROR! couldn't find aux node below in LocalVext near node %i",Nnum+1);
						t1 = h + 1.0; // set distance to h+1.0 to make algorithm work later
						vt1 = 0.0; // set velocity to 0.0, also for later
					}
				}
				else	{
					t1 = h + 1.0; // set distance to h+1.0 to make algorithm work later
					vt1 = 0.0; // set velocity to 0.0, also for later
				}

				if(t1 < tol)	{
					printf("\nERROR! for node %i in LocalVext t1 = %lf",Nnum+1,t1);
				}

				// ----------look at node above----------
				Nnum2 = Nodes[i][j+1];
				if(fabs(lsf[Nnum2]) < tol) // If neigbouring node is on the boundary
				{
					t2 = h; // distance to boundary = element edge length
					vt2 = Vnorm[Nnum2]; // store velocity of the node
				}
				else if((lsf[Nnum] * lsf[Nnum2]) < 0.0) // otherwise check to see is grid line is intersected
				{
					check = -1;

					for(k2=na_conn_ind[Nnum];k2<stop;k2++) // For all auxillary nodes (connected to Nnum node)
					{
						k = na_conn[k2]; // aux node num
						if( fabs(NodeCoord[Nnum].x - AuxNodes[k].x) < tol )	// If aux node has same x-coord
						{
							// If aux node lies between the two grid nodes
							if( (AuxNodes[k].y < NodeCoord[Nnum2].y) && (AuxNodes[k].y > NodeCoord[Nnum].y) )
							{
								t2 = fabs(AuxNodes[k].y - NodeCoord[Nnum].y);
								vt2 = Vnorm[NumNodes+k]; // store velocity of the auxillary node
								check = 1;
								break;
							}
						}
					}

					if(check == -1)
					{
						printf("\nERROR! couldn't find aux node above in LocalVext near node %i",Nnum+1);
						t2 = h + 1.0; // set distance to h+1.0 to make algorithm work later
						vt2 = 0.0; // set velocity to 0.0, also for later
					}

				}
				else	{
					t2 = h + 1.0; // set distance to h+1.0 to make algorithm work later
					vt2 = 0.0; // set velocity to 0.0, also for later
				}

				if(t2 < tol)	{
					printf("\nERROR! for node %i in LocalVext t2 = %lf",Nnum+1,t2);
				}

				// ----------look at node to left----------
				Nnum2 = Nodes[i-1][j];
				if(fabs(lsf[Nnum2]) < tol) // If neigbouring node is on the boundary
				{
					s1 = h; // distance to boundary = element edge length
					vs1 = Vnorm[Nnum2]; // store velocity of the node
				}
				else if((lsf[Nnum] * lsf[Nnum2]) < 0.0) // otherwise check to see is grid line is intersected
				{
					check = -1;

					for(k2=na_conn_ind[Nnum];k2<stop;k2++) // For all auxillary nodes (connected to Nnum node)
					{
						k = na_conn[k2]; // aux node num
						if( fabs(NodeCoord[Nnum].y - AuxNodes[k].y) < tol )// If aux node has same y-coord
						{
							// If aux node lies between the two grid nodes
							if( (AuxNodes[k].x < NodeCoord[Nnum].x) && (AuxNodes[k].x > NodeCoord[Nnum2].x) )
							{
								s1 = fabs(NodeCoord[Nnum].x - AuxNodes[k].x);
								vs1 = Vnorm[NumNodes+k]; // store velocity of the auxillary node
								check = 1;
								break;
							}
						}
					}

					if(check == -1)
					{
						printf("\nERROR! couldn't find aux node to left in LocalVext near node %i",Nnum+1);
						s1 = h + 1.0; // set distance to h+1.0 to make algorithm work later
						vs1 = 0.0; // set velocity to 0.0, also for later
					}

				}
				else	{
					s1 = h + 1.0; // set distance to h+1.0 to make algorithm work later
					vs1 = 0.0; // set velocity to 0.0, also for later
				}
				if(s1 < tol)	{
					printf("\nERROR! for node %i in LocalVext s1 = %lf",Nnum+1,s1);
				}

				// ----------look at node to right----------
				Nnum2 = Nodes[i+1][j];
				if(fabs(lsf[Nnum2]) < tol) // If neigbouring node is on the boundary
				{
					s2 = h; // distance to boundary = element edge length
					vs2 = Vnorm[Nnum2]; // store velocity of the node
				}
				else if((lsf[Nnum] * lsf[Nnum2]) < 0.0) // otherwise check to see is grid line is intersected
				{
					check = -1;

					for(k2=na_conn_ind[Nnum];k2<stop;k2++) // For all auxillary nodes (connected to Nnum node)
					{
						k = na_conn[k2]; // aux node num
						if( fabs(NodeCoord[Nnum].y - AuxNodes[k].y) < tol ) // If aux node has same y-coord
						{
							// If aux node lies between the two grid nodes
							if( (AuxNodes[k].x < NodeCoord[Nnum2].x)  && (AuxNodes[k].x > NodeCoord[Nnum].x) )
							{
								s2 = fabs(AuxNodes[k].x - NodeCoord[Nnum].x);
								vs2 = Vnorm[NumNodes+k]; // store velocity of the auxillary node
								check = 1;
								break;
							}
						}
					}

					if(check == -1)
					{
						printf("\nERROR! couldn't find aux node to right in LocalVext near node %i",Nnum+1);
						s2 = h + 1.0; // set distance to h+1.0 to make algorithm work later
						vs2 = 0.0; // set velocity to 0.0, also for later
					}

				}
				else	{
					s2 = h + 1.0; // set distance to h+1.0 to make algorithm work later
					vs2 = 0.0; // set velocity to 0.0, also for later
				}
				if(s2 < tol)	{
					printf("\nERROR! for node %i in LocalVext s2 = %lf",Nnum+1,s2);
				}

				// ---Now calculate lsf temp and Vnorm for node---
				//choose lowest s and t values, i.e. nearest to boundary
				if(fabs(s1 - s2) < tol) // if distances are the same then choose greatest velocity
				{
					s = s1; // could be s2, doesn't make much difference
					// choose velocity that will move lsf towards zero, or least distance away from zero
					vs = ((sign * vs1) > (sign * vs2)) ? vs1 : vs2;
				}
				else if(s1 < s2)	{
					s = s1;
					vs = vs1;
				}
				else if(s1 > s2) {
					s = s2;
					vs = vs2;
				}

				// t values
				if(fabs(t1 - t2) < tol) // if distances are the same then choose greatest velocity
				{
					t = t1;  // could be t2, doesn't make much difference
					// choose velocity that will move lsf towards zero, or least distance away from zero
					vt = ((sign * vt1) > (sign * vt2)) ? vt1 : vt2;
				}
				else if(t1 < t2)	{
					t = t1;
					vt = vt1;
				}
				else if(t1 > t2) {
					t = t2;
					vt = vt2;
				}

				// Calculate Extension velocity
				s2 = ((s - h) > 0.1) ? 0.0 : (1.0 / (s * s));
				t2 = ((t - h) > 0.1) ? 0.0 : (1.0 / (t * t)); // weighting for vs and vt velocities
				ftemp = s2 * vs;
				ftemp += t2 * vt;
				ftemp /= (s2 + t2); // Vnorm = Vext, by using inverse square weighting of distances
				Vtemp[Nnum] = ftemp;

				// Now calculate lsf temp
				if((s - h) > 0.1) {	// If horizontal distance > edge length only use vertical distance
					lsf_temp[Nnum] = t;
				}
				else if((t - h) > 0.1)	{	// If vertical distance > edge length only use horizontal distance
					lsf_temp[Nnum] = s;
				}
				else
				{
					s2 = s * s;
					t2 = t * t; // squared distances

					ftemp = s * t;
					ftemp /= sqrt(s2 + t2);
					lsf_temp[Nnum] = ftemp; // calculate perpendicular distance to the boundary (by pythag)
				}

				if((s > h) && (t > h)) {	// print error if both horizontal and vertical distances greater than edge length
					printf("\nERROR! For node %i, both s & t > h in LocalVext, s=%lf, t=%lf", Nnum+1,s,t);
				}
				// Move node into known set
				trial[Nnum] = false;
				known[Nnum] = true;
			}
		}
	}
}

// Function to calcualte gradient using WENO scheme
double CHakLevelSet::GradWENO(int Xi, int Yj, int num, double *lsf, mesh *inMesh, int sign, double Vn, double dt)
{
	// read data
	int NodeX = inMesh->NodeX;
	int NodeY = inMesh->NodeY;
	int **Nodes2 = inMesh->Nodes2;
	double h = inMesh->h;
	double tol = inMesh->tol;

	double v1,v2,v3,v4,v5; // varibales for finite differences used
	double grad, ftemp;
	int flag = 0;
	double phi = lsf[num]; // read in lsf value of current node

	if(Xi == 1) // If the node is on the left edge
	{
		if(Yj == 1) // node is bottom left corner
		{
			// if lsf of node above and to right have the same value as the corner node then use node on diagonal
			if( (fabs(lsf[Nodes2[Xi][Yj+1]] - phi) < tol) && (fabs(lsf[Nodes2[Xi+1][Yj]] - phi) < tol) )
			{
				grad = fabs(phi - lsf[Nodes2[Xi+1][Yj+1]]); // distance to diagonal node
				grad *= 1.414213562; // multiply by sqrt(2)
				grad /= h;	// divide by edge length
				flag = 1;	// flag gradient already calcualted
			}
		}

		else if(Yj == NodeY -2) // node is top left corner
		{
			// if lsf of node below and to right have the same value as the corner node then use node on diagonal
			if( (fabs(lsf[Nodes2[Xi][Yj-1]] - phi) < tol) && (fabs(lsf[Nodes2[Xi+1][Yj]] - phi) < tol) )
			{
				grad = fabs(phi - lsf[Nodes2[Xi+1][Yj-1]]); // distance to diagonal node
				grad *= 1.414213562; // multiply by sqrt(2)
				grad /= h;	// divide by edge length
				flag = 1;	// flag gradient already calcualted
			}
		}
	}

	else if(Xi == NodeX -2) // If the node is on the right edge
	{
		if(Yj == 1) // node is bottom right corner
		{
			// if lsf of node above and to left have the same value as the corner node then use node on diagonal
			if( (fabs(lsf[Nodes2[Xi][Yj+1]] - phi) < tol) && (fabs(lsf[Nodes2[Xi-1][Yj]] - phi) < tol) )
			{
				grad = fabs(phi - lsf[Nodes2[Xi-1][Yj+1]]); // distance to diagonal node
				grad *= 1.414213562; // multiply by sqrt(2)
				grad /= h;	// divide by edge length
				flag = 1;	// flag gradient already calcualted
			}
		}

		else if(Yj == NodeY -2) // node is top right corner
		{
			// if lsf of node below and to left have the same value as the corner node then use node on diagonal
			if( (fabs(lsf[Nodes2[Xi][Yj-1]] - phi) < tol) && (fabs(lsf[Nodes2[Xi-1][Yj]] - phi) < tol) )
			{
				grad = fabs(phi - lsf[Nodes2[Xi-1][Yj-1]]); // distance to diagonal node
				grad *= 1.414213562; // multiply by sqrt(2)
				grad /= h;	// divide by edge length
				flag = 1;	// flag gradient already calcualted
			}
		}
	}

	if(flag != 1)
	{
		// first assume info travelling to the right
		if(Xi == 1) // if node on left edge
		{
			v1 = (lsf[Nodes2[4][Yj]] - lsf[Nodes2[3][Yj]]) / h;
			v2 = (lsf[Nodes2[3][Yj]] - lsf[Nodes2[2][Yj]]) / h;
			v3 = (lsf[Nodes2[2][Yj]] - lsf[Nodes2[1][Yj]]) / h;
			v4 = v3; // approx gradient outside of domain
			v5 = v3; // approx gradient outside of domain
		}
		else if(Xi == 2)
		{
			v1 = (lsf[Nodes2[5][Yj]] - lsf[Nodes2[4][Yj]]) / h;
			v2 = (lsf[Nodes2[4][Yj]] - lsf[Nodes2[3][Yj]]) / h;
			v3 = (lsf[Nodes2[3][Yj]] - lsf[Nodes2[2][Yj]]) / h;
			v4 = (lsf[Nodes2[2][Yj]] - lsf[Nodes2[1][Yj]]) / h;
			v5 = v4; // approx gradient outside of domain
		}
		else if(Xi == NodeX -2) // if node on right edge
		{
			v5 = (lsf[Nodes2[Xi-1][Yj]] - lsf[Nodes2[Xi-2][Yj]]) / h;
			v4 = (lsf[Nodes2[Xi][Yj]] - lsf[Nodes2[Xi-1][Yj]]) / h;
			v3 = v4;
			v2 = v4;
			v1 = v4;	// approx gradients outside of domain
		}
		else if(Xi == NodeX -3)
		{
			v5 = (lsf[Nodes2[Xi-1][Yj]] - lsf[Nodes2[Xi-2][Yj]]) / h;
			v4 = (lsf[Nodes2[Xi][Yj]] - lsf[Nodes2[Xi-1][Yj]]) / h;
			v3 = (lsf[Nodes2[Xi+1][Yj]] - lsf[Nodes2[Xi][Yj]]) / h;
			v2 = v3;
			v1 = v3;	// approx gradients outside of domain
		}
		else if(Xi == NodeX -4)
		{
			v5 = (lsf[Nodes2[Xi-1][Yj]] - lsf[Nodes2[Xi-2][Yj]]) / h;
			v4 = (lsf[Nodes2[Xi][Yj]] - lsf[Nodes2[Xi-1][Yj]]) / h;
			v3 = (lsf[Nodes2[Xi+1][Yj]] - lsf[Nodes2[Xi][Yj]]) / h;
			v2 = (lsf[Nodes2[Xi+2][Yj]] - lsf[Nodes2[Xi+1][Yj]]) / h;
			v1 = v2;
		}
		else
		{
			v1 = (lsf[Nodes2[Xi+3][Yj]] - lsf[Nodes2[Xi+2][Yj]]) / h;
			v2 = (lsf[Nodes2[Xi+2][Yj]] - lsf[Nodes2[Xi+1][Yj]]) / h;
			v3 = (lsf[Nodes2[Xi+1][Yj]] - lsf[Nodes2[Xi][Yj]]) / h;
			v4 = (lsf[Nodes2[Xi][Yj]] - lsf[Nodes2[Xi-1][Yj]]) / h;
			v5 = (lsf[Nodes2[Xi-1][Yj]] - lsf[Nodes2[Xi-2][Yj]]) / h;
		}

		// calcualte gradient, info going to the right
		double Xr = (double)sign * GWsub(v1, v2, v3, v4, v5);

		// Now assume info travelling to the left
		if(Xi == NodeX -2) // if node on right edge
		{
			v1 = (lsf[Nodes2[Xi-2][Yj]] - lsf[Nodes2[Xi-3][Yj]]) / h;
			v2 = (lsf[Nodes2[Xi-1][Yj]] - lsf[Nodes2[Xi-2][Yj]]) / h;
			v3 = (lsf[Nodes2[Xi][Yj]] - lsf[Nodes2[Xi-1][Yj]]) / h;
			v4 = v3; // approx gradient outside of domain
			v5 = v3; // approx gradient outside of domain
		}
		else if(Xi == NodeX -3)
		{
			v1 = (lsf[Nodes2[Xi-2][Yj]] - lsf[Nodes2[Xi-3][Yj]]) / h;
			v2 = (lsf[Nodes2[Xi-1][Yj]] - lsf[Nodes2[Xi-2][Yj]]) / h;
			v3 = (lsf[Nodes2[Xi][Yj]] - lsf[Nodes2[Xi-1][Yj]]) / h;
			v4 = (lsf[Nodes2[Xi+1][Yj]] - lsf[Nodes2[Xi][Yj]]) / h;
			v5 = v4; // approx gradient outside of domain
		}
		else if(Xi == 1) // if node on left edge
		{
			v5 = (lsf[Nodes2[3][Yj]] - lsf[Nodes2[2][Yj]]) / h;
			v4 = (lsf[Nodes2[2][Yj]] - lsf[Nodes2[1][Yj]]) / h;
			v3 = v4;
			v2 = v4;
			v1= v4;	// approx gradients outside of domain
		}
		else if(Xi == 2)
		{
			v5 = (lsf[Nodes2[4][Yj]] - lsf[Nodes2[3][Yj]]) / h;
			v4 = (lsf[Nodes2[3][Yj]] - lsf[Nodes2[2][Yj]]) / h;
			v3 = (lsf[Nodes2[2][Yj]] - lsf[Nodes2[1][Yj]]) / h;
			v2 = v3;
			v1= v3;	// approx gradients outside of domain
		}
		else if(Xi == 3)
		{
			v5 = (lsf[Nodes2[5][Yj]] - lsf[Nodes2[4][Yj]]) / h;
			v4 = (lsf[Nodes2[4][Yj]] - lsf[Nodes2[3][Yj]]) / h;
			v3 = (lsf[Nodes2[3][Yj]] - lsf[Nodes2[2][Yj]]) / h;
			v2 = (lsf[Nodes2[2][Yj]] - lsf[Nodes2[1][Yj]]) / h;
			v1 = v2;	// approx gradients outside of domain
		}
		else
		{
			v5 = (lsf[Nodes2[Xi+2][Yj]] - lsf[Nodes2[Xi+1][Yj]]) / h;
			v4 = (lsf[Nodes2[Xi+1][Yj]] - lsf[Nodes2[Xi][Yj]]) / h;
			v3 = (lsf[Nodes2[Xi][Yj]] - lsf[Nodes2[Xi-1][Yj]]) / h;
			v2 = (lsf[Nodes2[Xi-1][Yj]] - lsf[Nodes2[Xi-2][Yj]]) / h;
			v1 = (lsf[Nodes2[Xi-2][Yj]] - lsf[Nodes2[Xi-3][Yj]]) / h;
		}

		double Xl = (double)sign * GWsub(v1, v2, v3, v4, v5); // calcualte gradient info going to the left

		// Now assume info travelling upward
		if(Yj == 1) // if node on bottom edge
		{
			v1 = (lsf[Nodes2[Xi][4]] - lsf[Nodes2[Xi][3]]) / h;
			v2 = (lsf[Nodes2[Xi][3]] - lsf[Nodes2[Xi][2]]) / h;
			v3 = (lsf[Nodes2[Xi][2]] - lsf[Nodes2[Xi][1]]) / h;
			v4 = v3; // approx gradient outside of domain
			v5 = v3; // approx gradient outside of domain
		}
		else if(Yj == 2)
		{
			v1 = (lsf[Nodes2[Xi][5]] - lsf[Nodes2[Xi][4]]) / h;
			v2 = (lsf[Nodes2[Xi][4]] - lsf[Nodes2[Xi][3]]) / h;
			v3 = (lsf[Nodes2[Xi][3]] - lsf[Nodes2[Xi][2]]) / h;
			v4 = (lsf[Nodes2[Xi][2]] - lsf[Nodes2[Xi][1]]) / h;
			v5 = v4; // approx gradient outside of domain
		}
		else if(Yj == NodeY -2) // if node on top edge
		{
			v5 = (lsf[Nodes2[Xi][Yj-1]] - lsf[Nodes2[Xi][Yj-2]]) / h;
			v4 = (lsf[Nodes2[Xi][Yj]] - lsf[Nodes2[Xi][Yj-1]]) / h;
			v3 = v4;
			v2 = v4;
			v1 = v4;	// approx gradients outside of domain
		}
		else if(Yj == NodeY -3)
		{
			v5 = (lsf[Nodes2[Xi][Yj-1]] - lsf[Nodes2[Xi][Yj-2]]) / h;
			v4 = (lsf[Nodes2[Xi][Yj]] - lsf[Nodes2[Xi][Yj-1]]) / h;
			v3 = (lsf[Nodes2[Xi][Yj+1]] - lsf[Nodes2[Xi][Yj]]) / h;
			v2 = v3;
			v1 = v3;	// approx gradients outside of domain
		}
		else if(Yj == NodeY -4)
		{
			v5 = (lsf[Nodes2[Xi][Yj-1]] - lsf[Nodes2[Xi][Yj-2]]) / h;
			v4 = (lsf[Nodes2[Xi][Yj]] - lsf[Nodes2[Xi][Yj-1]]) / h;
			v3 = (lsf[Nodes2[Xi][Yj+1]] - lsf[Nodes2[Xi][Yj]]) / h;
			v2 = (lsf[Nodes2[Xi][Yj+2]] - lsf[Nodes2[Xi][Yj+1]]) / h;
			v1 = v2;
		}
		else
		{
			v5 = (lsf[Nodes2[Xi][Yj-1]] - lsf[Nodes2[Xi][Yj-2]]) / h;
			v4 = (lsf[Nodes2[Xi][Yj]] - lsf[Nodes2[Xi][Yj-1]]) / h;
			v3 = (lsf[Nodes2[Xi][Yj+1]] - lsf[Nodes2[Xi][Yj]]) / h;
			v2 = (lsf[Nodes2[Xi][Yj+2]] - lsf[Nodes2[Xi][Yj+1]]) / h;
			v1 = (lsf[Nodes2[Xi][Yj+3]] - lsf[Nodes2[Xi][Yj+2]]) / h;
		}

		// calcualte gradient info going up
		double Yu = (double)sign * GWsub(v1, v2, v3, v4, v5);

		// Now assume info travelling down
		if(Yj == NodeY -2) // if node on right edge
		{
			v1 = (lsf[Nodes2[Xi][Yj-2]] - lsf[Nodes2[Xi][Yj-3]]) / h;
			v2 = (lsf[Nodes2[Xi][Yj-1]] - lsf[Nodes2[Xi][Yj-2]]) / h;
			v3 = (lsf[Nodes2[Xi][Yj]] - lsf[Nodes2[Xi][Yj-1]]) / h;
			v4 = v3; // approx gradient outside of domain
			v5 = v3; // approx gradient outside of domain
		}
		else if(Yj == NodeY -3)
		{
			v1 = (lsf[Nodes2[Xi][Yj-2]] - lsf[Nodes2[Xi][Yj-3]]) / h;
			v2 = (lsf[Nodes2[Xi][Yj-1]] - lsf[Nodes2[Xi][Yj-2]]) / h;
			v3 = (lsf[Nodes2[Xi][Yj]] - lsf[Nodes2[Xi][Yj-1]]) / h;
			v4 = (lsf[Nodes2[Xi][Yj+1]] - lsf[Nodes2[Xi][Yj]]) / h;
			v5 = v4; // approx gradient outside of domain
		}
		else if(Yj == 1) // if node on bottom edge
		{
			v5 = (lsf[Nodes2[Xi][3]] - lsf[Nodes2[Xi][2]]) / h;
			v4 = (lsf[Nodes2[Xi][2]] - lsf[Nodes2[Xi][1]]) / h;
			v3 = v4;
			v2 = v4;
			v1= v4;	// approx gradients outside of domain
		}
		else if(Yj == 2)
		{
			v5 = (lsf[Nodes2[Xi][4]] - lsf[Nodes2[Xi][3]]) / h;
			v4 = (lsf[Nodes2[Xi][3]] - lsf[Nodes2[Xi][2]]) / h;
			v3 = (lsf[Nodes2[Xi][2]] - lsf[Nodes2[Xi][1]]) / h;
			v2 = v3;
			v1= v3;	// approx gradients outside of domain
		}
		else if(Yj == 3)
		{
			v5 = (lsf[Nodes2[Xi][5]] - lsf[Nodes2[Xi][4]]) / h;
			v4 = (lsf[Nodes2[Xi][4]] - lsf[Nodes2[Xi][3]]) / h;
			v3 = (lsf[Nodes2[Xi][3]] - lsf[Nodes2[Xi][2]]) / h;
			v2 = (lsf[Nodes2[Xi][2]] - lsf[Nodes2[Xi][1]]) / h;
			v1 = v2;	// approx gradients outside of domain
		}
		else
		{
			v1 = (lsf[Nodes2[Xi][Yj-2]] - lsf[Nodes2[Xi][Yj-3]]) / h;
			v2 = (lsf[Nodes2[Xi][Yj-1]] - lsf[Nodes2[Xi][Yj-2]]) / h;
			v3 = (lsf[Nodes2[Xi][Yj]] - lsf[Nodes2[Xi][Yj-1]]) / h;
			v4 = (lsf[Nodes2[Xi][Yj+1]] - lsf[Nodes2[Xi][Yj]]) / h;
			v5 = (lsf[Nodes2[Xi][Yj+2]] - lsf[Nodes2[Xi][Yj+1]]) / h;
		}

		 // calcualte gradient info going down
		double Yd = (double)sign * GWsub(v1, v2, v3, v4, v5);

		// now calculate final gradinet using upwind scheme
		ftemp = (Xl > 0.0) ? Xl : 0.0;
		ftemp *= ftemp;
		grad = ftemp;

		ftemp = (Xr < 0.0) ? Xr : 0.0;
		ftemp *= ftemp;
		grad += ftemp;

		ftemp = (Yd > 0.0) ? Yd : 0.0;
		ftemp *= ftemp;
		grad += ftemp;

		ftemp = (Yu < 0.0) ? Yu : 0.0;
		ftemp *= ftemp;
		grad += ftemp;

		ftemp = grad;
		grad = sqrt(ftemp); // need to square root answer
	}

	else	{
		printf("\nGradient approimxated by diagonal node in corner");
	}

	// need to check that gradient won't take the boundary outside the domain
	if( (Yj == 1) || (Yj == NodeY -2) || (Xi == 1) || (Xi == NodeX -2) ) // if node on edge of domain
	{
		if( (phi - (grad * Vn * dt)) > 0.0 ) // if updated hpi value is positive
		{
			grad = phi;
			grad /= (Vn * dt);	// ensure boundary only moves to the edge of the domain
		}
	}
	return(grad);
}

// sub-function for GradWENO
double CHakLevelSet::GWsub(double v1,double v2,double v3,double v4,double v5)
{
	double ftemp;

	// Calculate the smoothness of approx 1
	ftemp = (v1 - (2.0 * v2) + v3);
	ftemp *= ftemp;
	double s1 = 13.0 * ftemp;
	ftemp = (v1 - (4.0 * v2) + (3.0 * v3));
	ftemp *= ftemp;
	s1 += 3.0 * ftemp;

	// Calculate the smoothness of approx 2
	ftemp = (v2 - (2.0 * v3) + v4);
	ftemp *= ftemp;
	double s2 = 13.0 * ftemp;
	ftemp = v2 - v4;
	ftemp *= ftemp;
	s2 += 3.0 * ftemp;

	// Calculate the smoothness of approx 3
	ftemp = (v3 - (2.0 * v4) + v5);
	ftemp *= ftemp;
	double s3 = 13.0 * ftemp;
	ftemp = (v5 - (4.0 * v4) + (3.0 * v3));
	ftemp *= ftemp;
	s3 += 3.0 * ftemp;

	double eps = 0.00001;

	// calcualte alpha weighting values
	ftemp = (s1 < eps) ? (eps * eps) : (s1 * s1);
	double a1 = 1.0 / ftemp;
	ftemp = (s2 < eps) ? (eps * eps) : (s2 * s2);
	double a2 = 6.0 / ftemp;
	ftemp = (s3 < eps) ? (eps * eps) : (s3 * s3);
	double a3 = 3.0 / ftemp;

	double asum = a1 + a2 + a3;

	// calculate normalized weightings
	a1 /= asum;
	a2 /= asum;
	a3 /= asum;

	// Finally calculate the approximate gradient
	double grad = a1 * ( (2.0 * v1) - (7.0 * v2) + (11.0 * v3) );
	grad += a2 * ( (2.0 * v4) - (v2) + (5.0 * v3) );
	grad += a3 * ( (2.0 * v3) - (v5) + (5.0 * v4) );
	grad *= 0.166666666667;

	return(grad);
}

// Function to re-initalise the lsf as a signed distance function - similar to Vext above
void CHakLevelSet::ReInt(mesh *inMesh, levSet *levelset)
{
	// read in data
	int NodeX = inMesh->NodeX - 1;
	int NodeY = inMesh->NodeY - 1;
	int **Nodes2 = inMesh->Nodes2;
	int NumNodes = inMesh->NumNodes;
	double tol = inMesh->tol;
	double *lsf = levelset->lsf;
	double h = inMesh->h;

	//array to temp store lsf function during re-initalisation calcs
	double *lsf_temp = (double *) malloc(NumNodes * sizeof(double));

	int i,j,k,sign; // incrementors
	int Nnum,Nnum2,Nnum3; // node numbers
	int Ci,Cj; // variables to store x,y location of node with minimum lsf_temp values in trial set
	double lsf_min,NS; // variable to track minimum lsf_temp value in trial set
	double Af, Bf, Cf; // variables for lsf_temp calculation by solution of a quadratic equn.
	double ftemp, dtemp; // tempory variables

	bool *known = (bool *) calloc(NumNodes,sizeof(bool));
	bool *trial = (bool *) calloc(NumNodes,sizeof(bool));

	// first ensure all nodes on domain boundary have lsf <= 0.0
	// bottom edge
	for(i=1;i<NodeX;i++)
	{
		Nnum = Nodes2[i][1]; // read in node number
		if(lsf[Nnum] > -tol)	{
			lsf[Nnum] = 0.0; // set lsf value to zero if greater than zero
		}
	}
	// top edge
	k = NodeY-1;
	for(i=1;i<NodeX;i++)
	{
		Nnum = Nodes2[i][k]; // read in node number
		if(lsf[Nnum] > -tol)	{
			lsf[Nnum] = 0.0; // set lsf value to zero if greater than zero
		}
	}
	// left edge
	for(i=1;i<NodeY;i++)
	{
		Nnum = Nodes2[1][i]; // read in node number
		if(lsf[Nnum] > -tol)	{
			lsf[Nnum] = 0.0; // set lsf value to zero if greater than zero
		}
	}
	// right edge
	k = NodeX-1;
	for(i=1;i<NodeY;i++)
	{
		Nnum = Nodes2[k][i]; // read in node number
		if(lsf[Nnum] > -tol)	{
			lsf[Nnum] = 0.0; // set lsf value to zero if greater than zero
		}
	}

	// Initalise Known set
	for(i=0;i<NumNodes;i++)
	{
		// If node is fixed or on boundary, then retain the lsf
		if(fabs(lsf[i]) < tol)	{
			known[i] = true;
			lsf_temp[i] = lsf[i]; // or 0.0 ??
		}
	}

	// put extra loop here for each diraction, changing sign 1, then -1
	for(sign=1;sign>-2;sign-=2)
	{
		ftemp = h*1000.0;
		// Initalise Trial set
		for(j=1;j<NodeY;j++)
		{
			for(i=1;i<NodeX;i++)
			{
				Nnum = Nodes2[i][j]; // read in node number

				// If current node not in known set and is the correct side of the boundary
				if(!known[Nnum] && (sign * lsf[Nnum]) > 0.0)
				{
					NS = ftemp * lsf[Nnum];
					// If any neighbouring node is on opposite side of the boundary, or on boundary, then it is a trial node
					if((lsf[Nodes2[i-1][j]] * NS) < 0.0) {
						trial[Nnum] = true;
						}
					else if((lsf[Nodes2[i+1][j]] * NS) < 0.0) {
						trial[Nnum] = true;
						}
					else if((lsf[Nodes2[i][j-1]] * NS) < 0.0) {
						trial[Nnum] = true;
						}
					else if((lsf[Nodes2[i][j+1]] * NS) < 0.0) {
						trial[Nnum] = true;
						}
				}
			}
		}

		// calculate lsf_temp for all current trial nodes
		LocalInt(NodeX, NodeY, Nodes2, lsf, lsf_temp, known, trial, h, tol, sign);

		// update trial set by considering neighbours of nodes in known set
		for(j=1;j<NodeY;j++)
		{
			for(i=1;i<NodeX;i++)
			{
				Nnum = Nodes2[i][j]; // read in node number
				// If current node not in known set and is the correct side of the boundary
				if(!known[Nnum] && (sign * lsf[Nnum]) > 0.0)
				{
					// If any neighbouring node is in known set then curret node is a trial node
					if(known[Nodes2[i-1][j]]) {
						trial[Nnum] = true;
						}
					else if(known[Nodes2[i+1][j]]) {
						trial[Nnum] = true;
						}
					else if(known[Nodes2[i][j-1]]) {
						trial[Nnum] = true;
						}
					else if(known[Nodes2[i][j+1]]) {
						trial[Nnum] = true;
						}
				}
			}
		}

		// Do until trial set is empty
		short flag;
		double *lsf_trial;
		lsf_trial = (double *) calloc(NumNodes,sizeof(double)); // array to store trial values of the lsf_temp during the fast marching method

		// set two arrays to store nodes to be updatedduring each iteration of the fast marching method
		int *xind, *yind, ncount;
		xind = (int *) malloc(NumNodes * sizeof(int));
		yind = (int *) malloc(NumNodes * sizeof(int));

		double lsf_max = (NodeX > NodeY) ? (100.0 * h * NodeX) : (100.0 * h * NodeY);
		do
		{
			flag = 0;
			lsf_min = lsf_max; // initalize to a value much bigger than the element edge length
			// For all trial nodes calculate lsf temp
			for(j=1;j<NodeY;j++)
			{
				for(i=1;i<NodeX;i++)
				{
					Nnum = Nodes2[i][j]; // read in node number
					// If current node not in known set, but is in the narrow band active set and is the correct side of the boundary
					if(trial[Nnum])
					{
						flag = 1;
						if(fabs(lsf_trial[Nnum]) < tol)
						{
							Af = 0.0; Bf = 0.0; Cf = 0.0;
							ftemp = 0.0; // Initalise values

							// look at all neighbouring nodes
							Nnum2 = Nodes2[i][j-1];
							Nnum3 = Nodes2[i][j+1]; // read in node number (above and below)

							// If neigbouring node is in known set update co-efficents for lsf_temp calc
							// ensure only closest node to boundary is used for upwind scheme
							if(known[Nnum2] && known[Nnum3])	{
								Af += 1.0;
								ftemp =  ( (fabs(lsf_temp[Nnum2])) < (fabs(lsf_temp[Nnum3])) ) ? lsf_temp[Nnum2] : lsf_temp[Nnum3]; // choose closest node
							}
							else if(known[Nnum2])	{
								Af += 1.0;
								ftemp = lsf_temp[Nnum2];
							}
							else if(known[Nnum3])	{
								Af += 1.0;
								ftemp = lsf_temp[Nnum3];
							}

							Bf += sign * ftemp;
							Cf += ftemp * ftemp;	// Update quadratic co-efficients

							Nnum2 = Nodes2[i-1][j];
							Nnum3 = Nodes2[i+1][j]; // read in node number (left and right)
							ftemp = 0.0; // re-initalise for other direction
							// If neigbouring node is in known set update co-efficents for lsf_temp calc
							// ensure only closest node to boundary is used for upwind scheme
							if(known[Nnum2] && known[Nnum3])	{
								Af += 1.0;
								ftemp =  ( (fabs(lsf_temp[Nnum2])) < (fabs(lsf_temp[Nnum3])) ) ? lsf_temp[Nnum2] : lsf_temp[Nnum3];  // choose closest node
							}
							else if(known[Nnum2])	{
								Af += 1.0;
								ftemp = lsf_temp[Nnum2];
							}
							else if(known[Nnum3])	{
								Af += 1.0;
								ftemp = lsf_temp[Nnum3];
							}

							Bf += sign * ftemp;
							Cf += ftemp * ftemp;	// Update quadratic co-efficients

							// Need to do final update of the Bf and Cf co-efficients
							Bf *= 2.0;
							Cf -= h * h;

							// Calculate lsf_temp by solving the quadratic equation
							ftemp = Bf * Bf;
							ftemp -= 4.0 * Af * Cf;
							dtemp = sqrt(ftemp);
							ftemp = Bf + dtemp;
							lsf_trial[Nnum] = ftemp / (2.0 * Af);
							if(Af < tol){ printf("\nError in ReInt Af = 0!"); }

							if( (lsf_trial[Nnum] < tol) || (lsf_trial[Nnum] > lsf_max))
							{
								printf("\nERROR! In ReInt lsf_trial = %lf",lsf_trial[Nnum]);
								printf("\nAf=%lf,Bf=%lf,Cf=%lf",Af,Bf,Cf);
							}
						}

						// If current lsf_temp distance is less than current minimum, update minimum and store node number
						if(lsf_trial[Nnum] < lsf_min)
						{
							lsf_min = lsf_trial[Nnum];
						}
					}
				}
			}

			if(flag == 0) {	break; }// if there are no trial nodes, then end the loop

			ncount = 0; // re-initialize
			// check to see which nodes need updating this iteration
			for(j=1;j<NodeY;j++)
			{
				for(i=1;i<NodeX;i++)
				{
					Nnum = Nodes2[i][j]; // read in node number

					if(lsf_trial[Nnum] < 0.0)	{
						printf("\nERROR IN ReInt: lsf_trial = %lf, for node %i",lsf_trial[Nnum],Nnum+1);
					}

					// if a trial node has a trail lsf value <= minimum then update
					if( trial[Nnum] && (fabs(lsf_trial[Nnum] - lsf_min) < tol) )
					{
						lsf_temp[Nnum] = sign * lsf_trial[Nnum]; // update lsf_temp for node with lsf_min
						// Move node to known set
						xind[ncount] = i;
						yind[ncount++] = j; // store update node for next bit
					}
				}
			}

			// update known set for next iteration
			for(i=0;i<ncount;i++)
			{
				Ci = xind[i];
				Cj = yind[i]; // read in updated node indicators
				Nnum = Nodes2[Ci][Cj];	// read in node number
				trial[Nnum] = false;
				known[Nnum] = true; // update known set
			}

			// update trial set for next iteration
			for(i=0;i<ncount;i++)
			{
				Ci = xind[i]; Cj = yind[i]; // read in updated node indicators

				// look at all neighbouring nodes
				Nnum2 = Nodes2[Ci][Cj-1];
				Nnum3 = Nodes2[Ci][Cj+1]; // read in node number (above and below)

				// Check to see if node below needs to be added to trial set
				if(!known[Nnum2] && sign * lsf[Nnum2] > 0.0)
				{
					trial[Nnum2] = true;
					lsf_trial[Nnum2] = 0.0; // reset trial value
				}
				// Check to see if node above needs to be added to trial set
				if(!known[Nnum3] && (sign * lsf[Nnum3]) > 0.0)
				{
					trial[Nnum3] = true;
					lsf_trial[Nnum3] = 0.0; // reset trial value
				}

				Nnum2 = Nodes2[Ci-1][Cj];
				Nnum3 = Nodes2[Ci+1][Cj]; // read in node number (left and right)

				// Check to see if node to left needs to be added to trial set
				if(!known[Nnum2] && sign * lsf[Nnum2] > 0.0)
				{
					trial[Nnum2] = true;
					lsf_trial[Nnum2] = 0.0; // reset trial value
				}
				// Check to see if node to right needs to be added to trial set
				if(!known[Nnum3] &&  sign * lsf[Nnum3] > 0.0)
				{
					trial[Nnum3] = true;
					lsf_trial[Nnum3] = 0.0; // reset trial value
				}
			}

			if(ncount == 0)	{
				printf("\nERROR! ncount=0 in ReInt! Aborting");
			}
		}
		while(ncount > 0);

		free(lsf_trial);
		free(xind);
		free(yind);
	}

	// update the lsf
	cblas_dcopy(NumNodes, lsf_temp, 1, lsf, 1);

	free(lsf_temp);
	free(known);
	free(trial);
}

// Function to reinitalise lsf for inital set of trial nodes - similar to LocalVext above
void CHakLevelSet::LocalInt(int NodeX, int NodeY, int **Nodes, double *lsf, double *lsf_temp,
			  bool *known, bool *trial, double h, double tol, int sign)
{
	double s,s1,s2,t,t1,t2; // distance variables
	int i,j; // incrementors
	double ftemp; // temp varibale
	int Nnum,Nnum2; // node number variables
	double h2 = h * 0.5;

	for(j=1;j<NodeY;j++)
	{
		for(i=1;i<NodeX;i++)
		{
			Nnum = Nodes[i][j];

			if(trial[Nnum])
			{
				// look at each neighbouring node in turn

				// -----------look at node below-----------
				Nnum2 = Nodes[i][j-1];
				if(fabs(lsf[Nnum2]) < tol) // If neigbouring node is on the boundary
				{
					t1 = h; // distance to boundary = element edge length
				}
				else if((lsf[Nnum] * lsf[Nnum2]) < 0.0) // otherwise check to see if grid line is intersected
				{
					t1 = (lsf[Nnum] + lsf[Nnum2]) / (lsf[Nnum] - lsf[Nnum2]);
					t1 += 1.0;
					t1 *= h2;
					// check to see if t1 > h, i.e. something has gone wrong in update
					// If so then approximate
					if((t1 - h) > -tol)
					{
						if( ((sign * -lsf[Nnum2]) - h) > -tol )	{
							ftemp = lsf[Nnum] - lsf[Nnum2];
							t1 = (lsf[Nnum] * h) / ftemp;
						}
						else	{
							t1 = (sign * lsf[Nnum2]) + h;
						}
						printf("\nERROR in LocalInt for node %i, t1 approx as %lf",Nnum+1,t1);
					}
				}
				else	{
					t1 = h + 1.0; // set distance to h+1.0 to make algorithm work later
				}

				// ----------look at node above----------
				Nnum2 = Nodes[i][j+1];
				if(fabs(lsf[Nnum2]) < tol) // If neigbouring node is on the boundary
				{
					t2 = h; // distance to boundary = element edge length
				}
				else if((lsf[Nnum] * lsf[Nnum2]) < 0.0) // otherwise check to see is grid line is intersected
				{
					t2 = (lsf[Nnum] + lsf[Nnum2]) / (lsf[Nnum] - lsf[Nnum2]);
					t2 += 1.0;
					t2 *= h2;
					// check to see if t2 > h, i.e. something has gone wrong in update
					// If so then approximate
					if((t2 - h) > -tol)
					{
						if( ((sign * -lsf[Nnum2]) - h) > -tol )	{
							ftemp = lsf[Nnum] - lsf[Nnum2];
							t2 = (lsf[Nnum] * h) / ftemp;
						}
						else	{
							t2 = (sign * lsf[Nnum2]) + h;
						}
						printf("\nERROR in LocalInt for node %i, t2 approx as %lf",Nnum+1,t2);
					}
				}
				else	{
					t2 = h + 1.0; // set distance to h+1.0 to make algorithm work later
				}

				// ----------look at node to left----------
				Nnum2 = Nodes[i-1][j];
				if(fabs(lsf[Nnum2]) < tol) // If neigbouring node is on the boundary
				{
					s1 = h; // distance to boundary = element edge length
				}
				else if((lsf[Nnum] * lsf[Nnum2]) < 0.0) // otherwise check to see is grid line is intersected
				{
					s1 = (lsf[Nnum] + lsf[Nnum2]) / (lsf[Nnum] - lsf[Nnum2]);
					s1 += 1.0;
					s1 *= h2;
					// check to see if s1 > h, i.e. something has gone wrong in update
					// If so then approximate
					if((s1 - h) > -tol)
					{
						if( ((sign * -lsf[Nnum2]) - h) > -tol )	{
							ftemp = lsf[Nnum] - lsf[Nnum2];
							s1 = (lsf[Nnum] * h) / ftemp;
						}
						else	{
							s1 = (sign * lsf[Nnum2]) + h;
						}
						printf("\nERROR in LocalInt for node %i, s1 approx as %lf",Nnum+1,s1);
					}
				}
				else	{
					s1 = h + 1.0; // set distance to h+1.0 to make algorithm work later
				}

				// ----------look at node to right----------
				Nnum2 = Nodes[i+1][j];
				if(fabs(lsf[Nnum2]) < tol) // If neigbouring node is on the boundary
				{
					s2 = h; // distance to boundary = element edge length
				}
				else if((lsf[Nnum] * lsf[Nnum2]) < 0.0) // otherwise check to see is grid line is intersected
				{
					s2 = (lsf[Nnum] + lsf[Nnum2]) / (lsf[Nnum] - lsf[Nnum2]);
					s2 += 1.0;
					s2 *= h2;
					// check to see if s2 > h, i.e. something has gone wrong in update
					// If so then approximate
					if((s2 - h) > -tol)
					{
						if( ((sign * -lsf[Nnum2]) - h) > -tol )	{
							ftemp = lsf[Nnum] - lsf[Nnum2];
							s2 = (lsf[Nnum] * h) / ftemp;
						}
						else	{
							s2 = (sign * lsf[Nnum2]) + h;
						}
						printf("\nERROR in LocalInt for node %i, s2 approx as %lf",Nnum+1,s2);
					}
				}
				else	{
					s2 = h + 1.0; // set distance to h+1.0 to make algorithm work later
				}

				// ---------Now calculate lsf temp for node---------
				// choose lowest s and t values, i.e. nearest to boundary
				s = (s1 < s2) ? s1 : s2;
				t = (t1 < t2) ? t1 : t2;

				// Now calculate lsf temp
				if((s - h) > 0.9) {	// If horizontal distance > edge length only use vertical distance
					lsf_temp[Nnum] = (double)sign * t;
				}
				else if((t - h) > 0.9)	{	// If vertical distance > edge length only use horizontal distance
					lsf_temp[Nnum] = (double)sign * s;
				}
				else
				{
					s2 = s * s;
					t2 = t * t; // squared distances

					ftemp = s * t;
					ftemp /= sqrt(s2 + t2);
					lsf_temp[Nnum] = (double)sign * ftemp; // calculate perpendicular distance to the boundary (by pythag)
				}

				if((s > h) && (t > h)) {	// print error if both horizontal and vertical distances greater than edge length
					printf("\nError! For node %i, both s & t > h in LocalInt, s=%lf, t=%lf", Nnum+1,s,t);
				}

				// Move node into known set
				trial[Nnum] = false;
				known[Nnum] = true;
			}
		}
	}
}

// ---- set of functions for SLP Level-set method ---- //

// obtain a boundary move vector from input of lam, s and move limits
void CHakLevelSet::get_delD(int n, int m, double *x, double *lam, double *s, double *up_lim, double *low_lim)
{
	// n = num boundary points
	// m = num functions
	// x = length n
	// s = length m x n (rows x colums)
	// lam = length m
	// limit - length n

	// get initial guess, x = s^T lam
	cblas_dgemv(CblasRowMajor, CblasTrans, m, n, 1.0, s, n, lam, 1, 0.0, x, 1);

	// modify x, wrt to the limits
	int i;
	int side = 0;
	for(i=0;i<n;i++)
	{
		side = (x[i] > up_lim[i] || x[i] < low_lim[i]) ? side+1 : side;
		x[i] = (x[i] > up_lim[i]) ? up_lim[i] : x[i];
		x[i] = (x[i] < low_lim[i]) ? low_lim[i] :  x[i];
	}

	double pct = (double)(side)/(double)(n);
	if(pct > 0.9){printf("\nWarning %12.4e pct side limit hit",pct);}
}

// function to obtain gradients by finite difference
void CHakLevelSet::get_slpGrad(int n, int m, int numVar, double *lam, double *s, double *c, double *up_lim,
				 double *low_lim, double *max_lam, double *min_lam, double *grad, int pinfo)
{
	double *x = (double *) malloc(n * sizeof(double));  // boundary move vector

	// define trial values of lam for finite difference
	double lam1,lam2;
	double *lam_temp = (double *) malloc(m * sizeof(double));
	double *del_lam = (double *) malloc(m * sizeof(double));

	// make initial guess
	int i,j;
	double inc = 1.0e-4;
	for(i=0;i<m;i++)
	{
		// set temp lam vector
		lam_temp[i] = lam[i];

		// set delta for lambda
		del_lam[i] = inc;
	}

	double *val1 = (double *) malloc(m * sizeof(double));
	double *val2 = (double *) malloc(m * sizeof(double));
	double delLam,ftemp;
	// evaluate each function for each value of lam
	for(i=0;i<m;i++) // for each lambda value (column in grad matrix)
	{
		lam1 = lam[i] + del_lam[i];
		lam2 = lam[i] - del_lam[i];

		// adjust for upper and lower limits
		lam1 = (lam1 > max_lam[i]) ? max_lam[i] : lam1;
		lam2 = (lam2 < min_lam[i]) ? min_lam[i] : lam2;

		// evaluate function for lam1 - fix other values
		lam_temp[i] = lam1;
		get_delD(n, m, x, lam_temp, s, up_lim, low_lim); // get boundary move vector

		for(j=0;j<m;j++) // for each function (row in grad matrix)
		{
			val1[j] = cblas_ddot(n, x, 1, &c[j*n], 1);
		}

		// evaluate function for lam2[i] - fix other values
		lam_temp[i] = lam2;
		get_delD(n, m, x, lam_temp, s, up_lim, low_lim); // get boundary move vector

		for(j=0;j<m;j++) // for each function (row in grad matrix)
		{
			val2[j] = cblas_ddot(n, x, 1, &c[j*n], 1);
		}

		// compute gradients for current lambda value
		delLam = lam1 - lam2;
		for(j=0;j<m;j++) // for each function (row in grad matrix)
		{
			ftemp = (val1[j] - val2[j]) / delLam;
			grad[numVar*j + i] = ftemp;

			if(pinfo==3) { printf("\ndlam%i / d%i = %12.4e",i,j,ftemp); }
		}

		// reset the lambda value
		lam_temp[i] = lam[i];
	}

	free(x);
	free(lam_temp);
	free(del_lam);
	free(val1);
	free(val2);
}

// sub-function for sorting in ascending order
int CHakLevelSet::dcmpfunc (const void * p1, const void * p2)
{
	double cd = ( *(double*)p1 - *(double*)p2 );
	int c = (cd < 0.0) ? -1 : 1;
	return c;
}

// function to get limits for lambda
void CHakLevelSet::getLamLim(int n, int m, double *inmax, double *inmin, double *s, double *c,
					double *up_lim, double *low_lim, int pinfo)
{
	int i,j,ind;
	int nr = m-1; // number of functions -1

	double denom_max,denom_min,ftemp;
	double sr_max, sr_min;

    // no constraints at all
    if(nr==0)
    {
        // find maximum and minimum limits
        inmax[0]=0.0; inmax[0]=0.0;
        for(i=0;i<n;i++)
        {
            if(fabs(s[i]) > 1.0e-6)
            {
                if(s[i] > 0.0)
                {
                    sr_max = up_lim[i] / s[i];
                    sr_min = low_lim[i] / s[i];
                }
                else
                {
                    sr_max = low_lim[i] / s[i];
                    sr_min = up_lim[i] / s[i];
                }
                inmax[0] = (sr_max > inmax[0]) ? sr_max : inmax[0];
                inmin[0] = (sr_min < inmin[0]) ? sr_min : inmin[0];
            }
        }
    }

    else
    {
        // first compute objective / constrant ratios
        double *rmax = (double *) malloc(nr*sizeof(double));
        double *rmin = (double *) malloc(nr*sizeof(double));

        double ob_max = cblas_ddot(n, c, 1, up_lim, 1);
        double ob_min = cblas_ddot(n, c, 1, low_lim, 1);

        for(i=1;i<m;i++)
        {
            j=i*n; // point to correct place in array
            ftemp = cblas_ddot(n, &c[j], 1, up_lim, 1);
            rmax[i-1] = (fabs(ftemp) > 1.0e-12) ? (ob_max / ftemp) : 1.0e12;

            ftemp = cblas_ddot(n, &c[j], 1, low_lim, 1);
            rmin[i-1] = (fabs(ftemp) > 1.0e-12) ? (ob_min / ftemp) : 1.0e12;
            if(pinfo==3){printf("\nFor const %i: rmax = %12.4e, rmin = %12.4e",i,rmax[i-1],rmin[i-1]);}
        }

        // then compute lam_max and lam_min for each point
        double *lammax = (double *) malloc(n*sizeof(double));
        double *lammin = (double *) malloc(n*sizeof(double));

        for(i=0;i<n;i++)
        {
            denom_max = s[i];
            denom_min = s[i];
            for(j=0;j<nr;j++)
            {
                ind = (j+1)*n + i; // indicate correct position in s

                sr_max = s[ind]*rmax[j];
                sr_min = s[ind]*rmin[j];

                // swap if necessary
                if(sr_max < sr_min)
                {
                    ftemp = sr_max;
                    sr_max = sr_min;
                    sr_min = ftemp;
                }

                denom_max += sr_max;
                denom_min += sr_min;
            }

            // need to choose denom correctly
            sr_max = up_lim[i] / denom_max;
            sr_min = up_lim[i] / denom_min;
            lammax[i] = (sr_min > sr_max) ? sr_min : sr_max;

            sr_max = low_lim[i] / denom_max;
            sr_min = low_lim[i] / denom_min;
            lammin[i] = (sr_min < sr_max) ? sr_min : sr_max;

            // final check on estimated lambda limits
            if(lammax[i] < lammin[i])
            {
                // swap if wrong way round
                sr_max = lammax[i];
                lammax[i] = lammin[i];
                lammin[i] = sr_max;
            }
        }

        // sort lammin & lammax in ascending order
        qsort(lammax, n, sizeof(double), dcmpfunc);
        qsort(lammin, n, sizeof(double), dcmpfunc);

        // choose 50% (median) into lammax & lammin
        i = (int)(0.5*n); // 50%
        j = (int)(0.5*n); // 50%

        inmax[0] = lammax[i];
        inmin[0] = lammin[j];

        // set limits for other lambda values
        for(i=1;i<m;i++)
        {
            // choose minimum absolute value
            ftemp = (fabs(rmax[i-1]) < fabs(rmin[i-1])) ? fabs(rmax[i-1]) : fabs(rmin[i-1]);

            inmax[i] = inmax[0] * ftemp;
            inmin[i] = inmin[0] * ftemp;

            if(inmax[i] < inmin[i])
            {
                ftemp = inmax[i];
                inmax[i] = inmin[i];
                inmin[i] = ftemp;
            }
        }

        if(pinfo==3) {
            for(i=0;i<m;i++) {
                printf("\nFor Lam %i: max = %12.4e , min = %12.4e",i+1,inmax[i],inmin[i]);
            }
        }

        // free memory
        free(rmax);
        free(rmin);
        free(lammax);
        free(lammin);
    }
}

// function to set velocity (or move dist) using SLP - filter method
int CHakLevelSet::SLPsubSol4(mesh *inMesh, levSet *levelset, boundary *bound_in, double alt, double *delCon, int numCon,
                double **sens, double *cA, int n, int *bound, double *lam_in, int *active, double *Vnorm, double *pred,
				int numAdd, double *add_sens, double *add_min, double *add_max, double *add_change, int pinfo)
{
	// read data
	double h = inMesh->h;
	double maxX = inMesh->maxX;
	double maxY = inMesh->maxY;
	int NumNodes = inMesh->NumNodes;
	Coord *NodeCoord = inMesh->NodeCoord;
    Coord *AuxNodes = bound_in->AuxNodes;

	int m = 1+numCon; // number of shape sensitivities (+1 for the objective)
	int numVar = m + numAdd; // total number of variables

	// n = number of boudnary points
	// cA = matrix of boundary integral coefficients (m x n, 1st row for objective)
	// sens = matrix of raw shape sensitvites (m x NumTot, 1st row for objective)
	// bound = boundary point numbers (length = n)
	// delCon = target change in constraints (length = numCon)

	int i,j,k,nd,err,count; // incementors
	double ftemp, ftemp2; // temp variable
	Coord Crd;
	double a1,del1,del2;

	double *sfg = (double *) malloc(m * n * sizeof(double)); // objective & constraint shape sens
	double *x = (double *) malloc(n * sizeof(double));  // boundary move vector
	double *u_low = (double *) malloc(n * sizeof(double)); // upper limit on x
	double *u_up = (double *) malloc(n * sizeof(double));  // lower limit on x
	double *b = (double *) malloc(numCon * sizeof(double)); // constraint change targets (sent to LP)
	double *lam = (double *) malloc(numVar * sizeof(double));	 // multiplers for each shape senstivity (+ add vars)
	double *maxS = (double *) calloc(m,sizeof(double)); // find abolute maximum sensitivty (of each function)
	double *maxC = (double *) calloc(m,sizeof(double)); // find abolute maximum integral coefficient (of each function)

	// condense the raw senstivity matrix
	for(i=0;i<n;i++)
	{
		nd = bound[i]; // node number

		for(j=0;j<m;j++) // for each function
		{
			k = (j*n)+i;
			sfg[k] = sens[j][nd]; // shape sens for function j
            //printf("\nsens[%i][%i]=%12.4e, cA=%12.4e",j,nd,sfg[k],cA[k]);

			// find maximums
			ftemp = fabs(sens[j][nd]);
			maxS[j] = (ftemp > maxS[j]) ? ftemp : maxS[j];
			ftemp = fabs(cA[k]);
			maxC[j] = (ftemp > maxC[j]) ? ftemp : maxC[j];
		}
	}

    // maximums for add sens (as well)
    for(j=0;j<m;j++) // for each function
	{
        k=j*numAdd;
        for(i=0;i<numAdd;i++)
        {
            ftemp = fabs(add_sens[k++]);
            maxC[j] = (ftemp > maxC[j]) ? ftemp : maxC[j];
        }
    }

	// normalize sfg (max magnitude = 1)
	for(j=0;j<m;j++) // for each function
	{
		k=(j*n);
		if(maxS[j] < 1.0e-20){maxS[j]=1.0;}
		if(maxC[j] < 1.0e-20){maxC[j]=1.0;}
		ftemp = 1.0 / maxS[j]; ftemp2 = 1.0 / maxC[j];
		if(pinfo==3){ printf("\nmaxS[%i]=%12.4e, maxC[%i]=%12.4e",j+1,maxS[j],j+1,maxC[j]); }
		for(i=0;i<n;i++)
		{
			sfg[k+i] *= ftemp; // scale shape sensitivities
			cA[k+i] *= ftemp2;  // also scale boundary integral coeffs
		}
		lam[j] = lam_in[j]; // read in start lambda values
        //lam[j] = 0.0; // EDIT

		k=j*numAdd;
		for(i=0;i<numAdd;i++)
		{
			add_sens[k++] *= ftemp2; // scale add var sensitivities
		}
	}

	for(j=1;j<m;j++)
	{
		delCon[j-1] /= maxC[j]; // need to also scale the constraint change targets
	}

	// set side constraints, special cosiderations for points near domain edge
	double max_delD = h * alt; // maximum boundary move value

	for(i=0;i<n;i++)
	{
		u_low[i] = -max_delD;
		u_up[i] = max_delD;

		j = bound[i]; // node number
		if(j < NumNodes) {
			Crd.x = NodeCoord[j].x;
			Crd.y = NodeCoord[j].y;
		}
		else {
			k = j - NumNodes;
			Crd.x = AuxNodes[k].x;
			Crd.y = AuxNodes[k].y;
		}

		// Now work out if that node lies within one grid spacing of the boundary of the domain
		if((Crd.x < h) || (Crd.x > (maxX - h)) || (Crd.y < h) || (Crd.y > (maxY - h)))
		{
			// find shortest distance from node to domain boundary
			a1 = maxX-Crd.x;
			del1 = (Crd.x < a1) ? Crd.x : a1;
			a1 = maxY - Crd.y;
			del2 = (Crd.y < a1) ? Crd.y : a1;
			a1 = (del1 < del2) ? del1 : del2;

			u_low[i] = (-a1 > u_low[i]) ? -a1 : u_low[i]; // modify lower limit
		}
	}

    // modify for aux nodes near fixed nodes
    if(levelset->fixed)
    {
        for(i=0;i<NumNodes;i++)
        {
            // if fixed, check for neighbouring aux nodes
            if(levelset->fixed[i])
            {
                nd = bound_in->na_conn_ind[i+1];
                for(j=bound_in->na_conn_ind[i+1];j<nd;j++)
                {
                    // compute dist from aux node to fixed node
                    k = bound_in->na_conn[j]; // aux node number
                    Crd.x = NodeCoord[j].x - AuxNodes[k].x; // x dist
                    Crd.y = NodeCoord[j].y - AuxNodes[k].y; // y dist
                    a1 = (Crd.x * Crd.x) + (Crd.y * Crd.y); // sqrd dist
                    a1 = sqrt(a1); // dist

                    // get variable number
                    count = k + NumNodes; // aux node num in bound
                    for(k=0;k<n;k++){ if(bound[k] == count){break;} }

                    // if fixed node is in - modify upper limit
                    if(levelset->lsf[i] > 0.0)
                    {
                        u_up[k] = (a1 < u_up[k]) ? a1 : u_up[k];
                    }
                    // otherwise - modify lower limit
                    else
                    {
                        u_low[k] = (-a1 > u_low[k]) ? -a1 : u_low[k];
                    }
                }
            }
        }
    }

	// estimate desired change in constraints to form b

	// first find max & min constrant changes (x c_mod)
	double c_mod = (numAdd>0) ? 0.1 : 0.2;
	err = 0;
	for(i=1;i<m;i++)
	{
		k=i*n; // point to correct place in cA array
		// del1 = max
		// del2 = min

		del1=0.0; del2=0.0;
		for(j=0;j<n;j++)
		{
			if(cA[k] > 0.0){
				del1 += cA[k]*u_up[j];
				del2 += cA[k]*u_low[j];
			}
			else {
				del1 += cA[k]*u_low[j];
				del2 += cA[k]*u_up[j];
			}
			k++;
		}

		del2 *= c_mod; // modify

		// add additional variable contributions
		if(numAdd>0)
		{
			k=i*numAdd;
			for(j=0;j<numAdd;j++)
			{
				if(add_sens[k] > 0.0){
					del1 += add_sens[k]*add_max[j];
					del2 += add_sens[k]*add_min[j]*c_mod;
				}
				else {
					del1 += add_sens[k]*add_min[j];
					del2 += add_sens[k]*add_max[j]*c_mod;
				}
				k++;
			}
		}

		if(delCon[i-1] < 0.0) // if reduction in constraint
		{
			b[i-1] = (delCon[i-1] < del2) ? del2 : delCon[i-1]; // set limit
			active[i-1] = 1; // active constraint
		}
		else // if increase (or no change) in constraint
		{
			if(delCon[i-1] > del1)
			{
				b[i-1] = del1;
				active[i-1] = 0; // not active
				err = 1; // need to reorganise some data
			}
			else
			{
				b[i-1] = delCon[i-1];
				active[i-1] = 2; // "nearly" active
			}
		}
		if(pinfo==3)
		{ printf("\ndel1=%12.4e, del2=%12.4e, delCon=%12.4e",del1,del2,delCon[i-1]);
		  printf("\nConstraint %i change target = %12.4e, %12.4e",i,b[i-1]*maxC[i],b[i-1]); }
	}

	// reduce arrays if some constraint not active
	int numAct = numCon;
	if(err == 1)
	{
		count = 0; // count active constraints
		for(i=0;i<numCon;i++)
		{
			if(active[i] != 0)
			{
				if(count != i)
				{
					// copy data
					k=(i+1)*n;
					nd=(count+1)*n;

					for(j=0;j<n;j++)
					{
						sfg[nd+j] = sfg[k+j]; // copy shape sensitivities
						cA[nd+j] = cA[k+j];  // copy boundary integral coeffs
					}

					k=(i+1)*numAdd;
					nd=(count+1)*numAdd;

					for(j=0;j<numAdd;j++)
					{
						add_sens[nd++] = add_sens[k++]; // copy add var sensitivities
					}

					// other data to copt accross
					maxS[count+1] = maxS[i+1];
					maxC[count+1] = maxC[i+1];
					b[count] = b[i];
					lam[count+1] = lam[i+1];
				}
				count++; // increase number of active constraint
			}
		}
		// other variables to set
		numAct = count;
		m = 1+numAct;
		numVar = m+numAdd;
	}

	printf("\nActive constraints = %i",numAct);

	// need default if numAct = 0 !!
	if(numAct == 0)
	{
		printf("\nNo active constraints");
		lam[0] = -0.5; // reduction in objective
		get_delD(n, m, x, lam, sfg, u_up, u_low);

		// use x to define Vnorm (at boundry points)
		for(i=0;i<n;i++)
		{
			j = bound[i]; // node number
			Vnorm[j] = x[i];
		}

		// pseudo move on add var
        for(i=0;i<numAdd;i++)
        {
            add_change[i] = (add_sens[i]>0.0) ? add_min[i] : add_max[i];
        }

		// save predicted value
		pred[0] = (cblas_ddot(n, x, 1, cA, 1)+cblas_ddot(numAdd, add_sens, 1, add_change, 1)) * maxC[0];
		printf("\nPredicted obj change = %12.4e", pred[0]);

		// clean-up and exit
		free(sfg);
		free(x);
		free(u_up);
		free(u_low);
		free(b);
		free(maxS);
		free(maxC);
		free(lam);
		return 0;
	}

	// if there are more than 10 variables, need to add slack variables for LP solve using primal-dual
    int numVarLP = numVar; // no slack
    if(numVar > 10)
    {
        // number of slack variables is equal to the number of active constriants (all inequalities)
        numVarLP = numVar + numAct; // total number of variables (inc slack) for LP solve
    }

	// get lambda limits
	double *min_lam = (double *) malloc(m * sizeof(double));
	double *max_lam = (double *) malloc(m * sizeof(double));
	getLamLim(n, m, max_lam, min_lam, sfg, cA, u_up, u_low, pinfo);

	// re-scale to make each lamba approx -1 -> +1
	for(j=0;j<m;j++) // for each function
	{
		k=(j*n);
		ftemp = (fabs(max_lam[j]) > fabs(min_lam[j])) ? fabs(max_lam[j]) : fabs(min_lam[j]);
		for(i=0;i<n;i++) { sfg[k++] *= ftemp; } // scale shape sensitivities

		min_lam[j] /= ftemp;
		max_lam[j] /= ftemp; // these should be magnitude near 1
	}

	// SLP loop to find "optimal" values for lam
	// lots of data arrays
	double *min_lam_step = (double *) calloc(numVarLP, sizeof(double));
	double *max_lam_step = (double *) malloc(numVarLP * sizeof(double)); // limits for current step size
	double *lim_lam = (double *) malloc(numVar * sizeof(double));      // limits for reduced step lengths in SLP
	double *dlam = (double *) malloc(m * numVarLP * sizeof(double));     // gradients of obj & constraint, wrt lam
	double *lam_trial = (double *) malloc(numVar * sizeof(double));    // trial values for lambda
	double *lam_step = (double *) malloc(numVarLP * sizeof(double));     // change in lambda
	double *lam_best =  (double *) malloc(numVar * sizeof(double));    // best solutionlambda
	double *b_step = (double *) malloc(numAct * sizeof(double));  // targets for constraint change, change !
	double *g_trial = (double *) malloc(numAct * sizeof(double)); // trial constraint values (using lam_trial)
	double *g_lam = (double *) malloc(numAct * sizeof(double));   // constraint values
	double *g_shift = (double *) malloc(numAct * sizeof(double)); // shift due to shift to lam_min=0 in LP
	double f_filter[200];
	double h_filter[200]; // filters for objective and constraints
	double f_best, h_best;

	double obj, act_obj, del_obj, pred_obj, pred_obj_trial, hcon, hcon_trial; // obj & constraint values
	int accept, redo; // flags for what is happening in algorithm
	int fcount = 1; // count number of filters

	// step size variables
	double step, step_max, step_min;
	step_min = 1.0e-4;
	step_max = 0.1;
	step = step_max; // initial step size

	// copy add var sensitivites & limits to dlam (will not change during SLP opt)
	//double *add_shift = calloc(numAct,sizeof(double)); // shift in constraints due to shift in add var
	if(numAdd > 0)
	{
		for(j=0;j<m;j++) // for each function
		{
			k=(j*numAdd); // pointer to add_sens
			nd=(j*numVarLP + m); // pointer to dlam

			for(i=0;i<numAdd;i++)
			{
				dlam[nd++] = add_sens[k++]; // copy across
			}
		}

		for(i=0;i<numAdd;i++)
		{
			lim_lam[i+m] = add_max[i] - add_min[i]; // copy across
			min_lam_step[i+m] = add_min[i];
		}
	}

	// if there are more than 10 variables, add "gradients" for slack variables
    if(numVar > 10)
    {
        // objective (does not get effetced by slack vars)
        for(j=numVar;j<numVarLP;j++){ dlam[j] = 0.0; }

        for(j=1;j<m;j++) // for each constraint (j = row in dLam)
		{
            nd = j*numVarLP + numVar; // start of row + jump to start of additional grads
            for(k=0;k<numAct;k++) // column in dLam
            {
                dlam[nd++] = (k+1 == j) ? 1.0 : 0.0; // identity matrix
            }
        }
    }

    // test
    /*lam_trial[0] = -1.5;
    lam_trial[1] = 0.0;
    for(i=0;i<300;i++)
    {
        // evaluate constraints & objective at lam_trail
        get_delD(n, m, x, lam_trial, sfg, u_up, u_low); // boundary movement vector
        obj = cblas_ddot(n, x, 1, cA, 1); // objective
        printf("\n%f , %lf",lam_trial[0],obj);
        lam_trial[0] += 0.01;
    }*/

	// choose initial lam for vel func
    if(numVar<5)
    {
		// minimise constraint violation to get start point for SLP
		get_lam0(n, m, numAct, lam, sfg, cA, b, min_lam, max_lam);
		err = con_min3(n, m, numAct, lam, sfg, &cA[n], b, min_lam, max_lam, u_up, u_low, pinfo);
    }
    else{ for(i=0;i<numVar;i++){lam[i] = 0.0;} } // initially zero

	for(i=0;i<m;i++){ lam_best[i] = lam[i]; }
    for(i=m;i<numVar;i++){ lam_best[i] = 0.0; lam[i]=0.0; } // initial change in add var is zero

	// evaluate constraints & objective at initial point (lam)
	get_delD(n, m, x, lam, sfg, u_up, u_low); // initial boundary movement vector
	cblas_dgemv(CblasRowMajor, CblasNoTrans, numAct, n, 1.0, &cA[n], n, x, 1, 0.0, g_lam, 1); // constraints
	obj = cblas_ddot(n, x, 1, cA, 1); // objective
	pred_obj = 1.0e+20; hcon=0.0; // initial values (to start algorithm)
	if(pinfo==3){printf("\nobj = %12.4e",obj);}
	f_best = 1.0e+20;

	// compute initial hcon value
	for(i=0;i<numAct;i++)
	{
		// normalize constraint violation (if target not small)
		if(pinfo==3){printf("\n%i g_lam = %12.4e",i+1,g_lam[i]);}
		ftemp = (g_lam[i]-b[i]) / fabs(b[i]); //*maxS[i+1];
		ftemp = (ftemp > 0.0) ? ftemp : 0.0; // inequality constraint
		if(pinfo==3){printf(" hcon = %f",ftemp);}
		hcon = (ftemp > hcon) ? ftemp : hcon; // choose maximum
	}
	f_filter[0] = -1.0e+12; h_filter[0] = (hcon < 1000.0) ? 1000.0 : hcon; // initial filters
	h_best = hcon;

	// get finite difference gradient approximations (at initial point)
	get_slpGrad(n, m, numVarLP, lam, sfg, cA, u_up, u_low, max_lam, min_lam, dlam, pinfo);

	int stop = 0; // flag to signal convergence
	int lp_info = (pinfo == 3) ? 2 : 0;
	count = 0; // counter for number of iterations
	redo = 0;
	// START loop
	do
	{
		if(pinfo==3){ printf("\nstep = %f",step); }

		// ------- solve LP sub-problem with step, lam ------ //

		// set limits for lambda step
		for(i=0;i<m;i++)
		{
			ftemp = lam[i] - step;
			min_lam_step[i] = (ftemp > min_lam[i]) ? -step : min_lam[i]-lam[i];
			ftemp = lam[i] + step;
			max_lam_step[i] = (ftemp < max_lam[i]) ? step : max_lam[i]-lam[i];

			// shift limit, so that lam_step min = 0: lam_new = lam + lam_step + min_lam_step
			lim_lam[i] = max_lam_step[i] - min_lam_step[i];
		}

		// modify constraint targets for LP sub-problem
		// shift in constraint for shift in lam_step
		cblas_dgemv(CblasRowMajor, CblasNoTrans, numAct, numVarLP, 1.0, &dlam[numVarLP], numVarLP, min_lam_step, 1, 0.0, g_shift, 1);

		for(i=0;i<numAct;i++)
		{
			// LP target = original - current - shift
			b_step[i] = b[i] - g_lam[i] - g_shift[i];
			if(pinfo==3){printf("\nLP constraint target %i = %f",i+1,b_step[i]);}
		}
		if(pinfo==3){printf("\n");}

		// sub-solve function for the trust region method (to get lam_step)
		if(numVar > 10)
        {
            i=0;
            ftemp = 0.5;
            do{
                err = getSolver().LPsolve(numVarLP, numAct, numVar, lam_step, dlam, &dlam[numVarLP], b_step, lim_lam, pinfo);
                if(err==-1){
                    for(j=0;j<numAct;j++){b_step[j] = ftemp*b[j] - g_lam[j] - g_shift[j];}
                    ftemp*=0.5;
                }
                i++;
            } while(err==-1 && i<5);
        }
		else {
			err = getSolver().lp_simplex(numVar, numAct, numVar, lam_step, dlam, &dlam[numVar], b_step, lim_lam, lp_info);
			if(err==-1 && numAdd==0){err = trust_sub(numVar, numAct, numVar, lam_step, dlam, &dlam[numVar], b_step, lim_lam, pinfo);}
		}
		for(i=0;i<m;i++)
		{
			lam_step[i] += min_lam_step[i]; // shift back
			if(pinfo==3){printf("\nLam_step%i = %12.4e",i+1,lam_step[i]);}
		}

		// -------- End of LP sub-solve phase -------- //

		// if the LP was feasible
		if(err != -1 && redo == 0)
		{
			// check change in lambda (convergence condition)
			ftemp = 0.0;
			for(i=0;i<m;i++)
			{
				ftemp = (fabs(lam_step[i]) > ftemp) ? fabs(lam_step[i]) : ftemp;
			}
			if(ftemp < 1.0e-6){stop=1;} // if small change, we have convergence!!
		}

		// if not converged yet
		if(stop != 1)
		{
			if(numVar<11 || err==0)
			{
				for(i=0;i<m;i++)
				{
					lam_trial[i] = lam[i] + lam_step[i]; // compute trial lambda values
					if(pinfo==3){printf("\nLam_trial %i = %f",i+1,lam_trial[i]);}
				}
				for(i=m;i<numVar;i++){ lam_step[i] += min_lam_step[i]; lam_trial[i] = lam_step[i];
					if(pinfo==3 && i<10){printf("\nLam_trial %i = %f",i+1,lam_trial[i]);}}

				get_delD(n, m, x, lam_trial, sfg, u_up, u_low);

				// compute actual & predicted objective change
				act_obj = cblas_ddot(n, x, 1, cA, 1); // actual new objective
				if(numAdd>0){ ftemp = cblas_ddot(numAdd, add_sens, 1, &lam_trial[m], 1);
								act_obj += ftemp; } // + for add variables
				del_obj = obj - act_obj; // actual reduction in obj value (+ve = reduction)
				pred_obj_trial = -cblas_ddot(m,lam_step,1,dlam,1); // predicted change (+ve = reduction)
                //printf("\npred_obj=%f, ftemp=%f",pred_obj_trial,ftemp);
				if(numAdd>0){ pred_obj_trial += cblas_ddot(numAdd, add_sens, 1, &lam[m], 1) - ftemp; }
				if(pinfo==3){printf("\nact_obj=%f, del_obj=%f, pred_obj=%f",act_obj,del_obj,pred_obj_trial);}

				// compute maximum magnitude of constraint violation (normalized)
				cblas_dgemv(CblasRowMajor, CblasNoTrans, numAct, n, 1.0, &cA[n], n, x, 1, 0.0, g_trial, 1);
				if(numAdd > 0)
				{
					// + for add variables
					for(i=0;i<numAct;i++)
					{
						k = (i+1)*numAdd; // point to correct location in add_sens
						g_trial[i] += cblas_ddot(numAdd, &add_sens[k], 1, &lam_trial[m], 1);
					}
				}

				hcon_trial = 0.0;
				for(i=0;i<numAct;i++)
				{
					// normalize constrant violation (if target not small)
					if(pinfo==3){printf("\n%i g_trial = %12.4e",i+1,g_trial[i]);}
					ftemp = (g_trial[i]-b[i]) / fabs(b[i]); //*maxS[i+1];
					ftemp = (ftemp > 0.0) ? ftemp : 1.0e-12; // inequality constraint
					if(pinfo==3){printf(" hcon = %12.4e",ftemp);}
					hcon_trial = (ftemp > hcon_trial) ? ftemp : hcon_trial; // choose maximum
				}
				if(hcon_trial < 1.0e-4){hcon_trial=1.0e-12;}
				ftemp2 = 1.0e-4*hcon*hcon; // squared x small
				if(pinfo==3){printf(" ftemp2 = %12.4e",ftemp2);}

				// check current values against the filter (and previous values)
				accept = 1;
				if(redo == 0)
				{
					if(pred_obj < ftemp2) // if previous an h-type iterate
					{
						// check acceptability to previous values
						// at least improvement in objective or constraint
						if(del_obj < (1.0e-4*hcon) && hcon_trial > 0.99*hcon){ accept = 0; }
						if(pinfo==3 && accept==0){printf("\nfail on prev filter");}
					}
				}
				if(accept == 1)
				{
					// check acceptability wrt the filter
					for(i=0;i<fcount;i++)
					{
						if(redo == 1){
							if(hcon_trial > 0.99*h_filter[i])
							{ accept = 0; break; }
						}
						// must dominate at least obj or constraints (maybe both)
						else if( act_obj > f_filter[i]-(1.0e-4*h_filter[i]) && hcon_trial > 0.99*h_filter[i])
						{ accept = 0; break; } // not acceptable
					}
					if(pinfo==3 && accept==0){printf("\nfail on current filter");}
				}

				// still need to check other conditions (trust region)
				if(redo == 0 && accept == 1)
				{
					redo=1;
					if(del_obj < 0.1*pred_obj_trial && pred_obj_trial >= ftemp2){ accept = 0; }
					if(pinfo==3 && accept==0){printf("\nfail on trust region");}
				}
			}
			else
			{
				accept = 0; // do not accept LPsolve that failed
			}

			// if not acceptable re-solve LP sub-problem
			if(accept == 0)
			{
				// also improve model by updating gradient using the new info
				step *= 0.5;

				/*if(numAdd > 0)
				{
					cblas_dgemv(CblasRowMajor, CblasNoTrans, numAct, n, 1.0, &cA[n], n, x, 1, 0.0, g_trial, 1);
				}
				for(i=0;i<m;i++) // row in dlam (function)
				{
					for(j=0;j<m;j++) // col in A
					{
						k = i*numVar + j;
						if(fabs(lam_step[j]) < 1.0e-4){ del1 = dlam[k]; }
						else
						{
							// pseudo gradient
							if(i>0) { del1 = (g_trial[i-1] - g_lam[i-1]) / lam_step[j]; }
							else { del1 = -del_obj / lam_step[j];}
						}
						if(pinfo==3){printf("\ndlam[%i,%i] %12.4e , %12.4e",i,j,dlam[k],del1);}
						dlam[k] += del1;
						dlam[k] *= 0.5; // average of linear + pseudo gradients
					}
				}*/
			}

			// otherwise accept the new point
			else
			{
				if(pinfo==3){ printf("\naccept point - pred_obj=%12.4e , ftemp2=%12.4e",pred_obj,ftemp2); }
				// if h-type iterate - add previous point to the filter
				if(pred_obj < ftemp2)
				{
					f_filter[fcount] = obj;
					h_filter[fcount++] = hcon;
					if(pinfo==3){ printf("\n%12.4e, %12.4e added to filter",obj,hcon); }
					if(fcount == 200){printf("\nWarn -filter reached 200!");}
				}

				// update obj and constraint violation values (& lam)
				pred_obj = pred_obj_trial;
				obj = act_obj;
				hcon = hcon_trial; // save for next iteration
				for(i=0;i<numVar;i++){ lam[i] = lam_trial[i]; }
				get_delD(n, m, x, lam, sfg, u_up, u_low); // initial boundary movement vector
				cblas_dgemv(CblasRowMajor, CblasNoTrans, numAct, n, 1.0, &cA[n], n, x, 1, 0.0, g_lam, 1);

				// update best point
				if( (hcon < 1.0e-4 && obj < f_best) || hcon < h_best)
				{
					f_best = obj;
					h_best = hcon;
					if(pinfo==3){ printf("\nupdate best lam"); }
					for(i=0;i<numVar;i++){ lam_best[i] = lam[i]; }
				}

				// get finite difference gradient approximations (at new point)
				get_slpGrad(n, m, numVarLP, lam, sfg, cA, u_up, u_low, max_lam, min_lam, dlam, pinfo);

				count++; // update count of accepted points

				// reset step to maximum
				step = step_max;
				redo = 0;
			}

			// if LP was not feasible
			if(step < step_min)
			{
				if(pinfo==3){ printf("\nbreak-out"); }
				count = 200;
			}
		}

		// reduce step max after some iterations
		if(count == 25){step_max = 0.05;}
		if(count == 50){step_max = 0.025;}
		if(count == 75){step_max = 0.0125;}

	} while (stop==0 && count < 100);
	// END loop

	if(h_best < 1.0e-4)
	{
		if(count == 200){ printf("\nSLP sub-solve stopped"); }
		else{printf("\nSLP sub-solve converged in %i itt",count);}
		err = 0;
		//if(count < 200)
		{
			lam_in[0] = lam_best[0];
			count = 0;
			for(i=0;i<numCon;i++)
			{
				if(active[i] != 0){ lam_in[i+1] = lam_best[count+1]; }
				else{ lam_in[i+1] = 0.0; }
			}
		}
	}
	else
	{
		err = 2;
		printf("\nSLP sub-solve did not converge!");
		for(i=0;i<=numCon;i++) { lam_in[i] = 0.0; } // re-set
	}

	// recover best solution
	for(i=0;i<m;i++){ lam[i] = lam_best[i]; printf("\nlam %i = %12.4e",i,lam[i]); }
	get_delD(n, m, x, lam, sfg, u_up, u_low);
	cblas_dgemv(CblasRowMajor, CblasNoTrans, numAct, n, 1.0, &cA[n], n, x, 1, 0.0, g_lam, 1);
	if(numAdd > 0)
	{
		// + for add variables
		for(i=0;i<numAct;i++)
		{
			k = (i+1)*numAdd; // point to correct location in add_sens
			g_lam[i] += cblas_ddot(numAdd, &add_sens[k], 1, &lam_best[m], 1);
		}
	}

	// use x to define Vnorm (at boundry points)
	for(i=0;i<n;i++)
	{
		j = bound[i]; // node number
		Vnorm[j] = x[i];
	}

	for(i=0;i<numAdd;i++)
	{
		add_change[i] = lam_best[i+m];
	}

	// save predicted values (used in realxation in main loop)
	// reset from potentially inactive constraints
	m = 1+numCon;
	pred[0] = f_best * maxC[0];
	printf("\nPredicted obj change = %12.4e", pred[0]);
	count = 0;
	for(i=1;i<m;i++)
	{
		if(active[i-1] != 0)
		{
			count++;
			pred[i] = g_lam[count-1] * maxC[count];
			printf("\nPredicted cnst %i change = %12.4e",i,pred[i]);
		}
	}

	// Clean up
	free(sfg);
	free(x);
	free(u_up);
	free(u_low);
	free(b);
	free(maxS);
	free(maxC);
	free(lam);
	free(lam_best);
	free(min_lam);
	free(max_lam);
	free(lim_lam);
	free(lam_trial);
	free(lam_step);
	free(g_lam);
	free(g_shift);
	free(g_trial);
	free(b_step);
	free(min_lam_step);
	free(max_lam_step);
	free(dlam);
	//free(add_shift);
	return err;
}

// sub-solve function for the trust region method
int CHakLevelSet::trust_sub(int n, int m, int nu, double *x, double *c, double *A, double *b, double *u, int pinfo)
{
	// n = num variables
	// m = num equality constraints
	// x = length n
	// c = length n
	// A = matrix (m x n)
	// b = length m
	// u = length n (upper limits on x) - assume 0 lower limit
	// assume n = m+1 (for now)

	int i,j,ind,ind2;
	char trans = 'T'; // transpose (for FORTRAN column major order)

	// first split the A matrix into 1st column and a square mxm matrix
	double *rhs = (double *) malloc(2*m*sizeof(double)); // first row = a1, second row = b
	double *Ar = (double *) malloc(m*m*sizeof(double));

	for(i=0;i<m;i++) // each row in A
	{
		ind = i*n; // indicator for A
		rhs[i] = A[ind++]; // first column
		rhs[i+m] = b[i]; // copy b into rhs

		ind2 = i*m; // indicator for Ar
		for(j=1;j<n;j++) // remaining columns
		{
			Ar[ind2++] = A[ind++]; // copy rest of row
		}
	}

	// next invert Ar with rhs (a1 & b)
	getSolver().dun_slvLPK(trans, m, 2, Ar, rhs);

	// check feasible limits on x[0]
	double x_max = u[0];
	double x_min = 0.0; // start with real limits
	double bc, ac, ftemp, ftemp2;

	for(i=1;i<n;i++)
	{
		ind = i-1;
		bc = rhs[m+ind]; ac = rhs[ind]; // linear coefficients

		// check x[0] max
		ftemp = bc/ac; // for x[i] = 0.0
		ftemp2 = (bc - u[i])/ac; //for x[i] = u[i]

		if(ftemp < ftemp2)
		{
			bc=ftemp; ftemp=ftemp2; ftemp2=bc; // swap so ftemp > ftemp2
		}

		x_max = (x_max < ftemp) ? x_max : ftemp; // choose smaller
		x_min = (x_min > ftemp2) ? x_min : ftemp2; // choose bigger
	}

	// infeasible - simply find solution closest to constraints
	int ret = 0;
	if(x_min > x_max)
	{
		if(pinfo==3){ printf("\nWarning x_min=%12.4e, x_max=%12.4e",x_min,x_max); }

		// try again !!!
		double *bd = (double *) malloc(m*sizeof(double));
		for(i=0;i<m;i++){ bd[i] = 0.5*u[i+1] - rhs[m+i]; }

		x[0] = cblas_ddot(m, rhs, 1, bd, 1);
		x[0] /= cblas_ddot(m, rhs, 1, rhs, 1);

		// then compute remaining variables
		for(i=1;i<n;i++)
		{
			ind = i-1;
			x[i] = rhs[m+ind] - x[0]*rhs[ind];
			if(pinfo==3){printf("\nx[%i] = %12.4e",i,x[i]);}
			if(x[i] < 0.0){ x[i] = 0.0; }
			else if(x[i] > u[i]){ x[i] = u[i]; }
			if(pinfo==3){printf(" after = %12.4e",x[i]);}
		}

		if(pinfo==3){printf("\nx[0] = %12.4e",x[0]);}
		if(x[0] < 0.0){ x[0] = 0.0; }
		else if (x[0] > u[0]){ x[0] = u[0]; }
		if(pinfo==3){printf("\tx[0] = %12.4e",x[0]);}

		ret = -1;
		free(bd);
	}
	// solve for 1st varible
	else
	{
		// compute gradient
		ftemp = c[0] - cblas_ddot(m, &c[1], 1, rhs, 1); // c[0] - cr * a1

		if(ftemp < 0.0){ x[0] = x_max; } // if -ve gradient, set x[0] to maximum
		else { x[0] = x_min; } // otherwise choose minimum

		if(pinfo==3){printf("\nx[0] = %12.4e",x[0]);}
		if(x[0] > u[0]){x[0] = u[0];}
		else if(x[0] < 0.0){x[0] = 0.0;}
		if(pinfo==3){printf("\tx[0] = %12.4e",x[0]);}

		// then compute remaining variables
		for(i=1;i<n;i++)
		{
			ind = i-1;
			x[i] = rhs[m+ind] - x[0]*rhs[ind];
			if(pinfo==3){printf("\nx[%i] = %12.4e",i,x[i]);}
			if(x[i] < 0.0){ x[i] = 0.0; }
			else if(x[i] > u[i]){ x[i] = u[i]; }
			if(pinfo==3){printf(" after = %12.4e",x[i]);}
		}
	}

	// finally return predicted constraint values (in b)
	if(pinfo==3){printf("\nB in:");
		for(i=0;i<m;i++) {printf(" %12.4e",b[i]);} }
	cblas_dgemv(CblasRowMajor, CblasNoTrans, m, n, 1.0, A, n, x, 1, 0.0, b, 1);
	if(pinfo==3){ printf("\nB out:");
		for(i=0;i<m;i++) {printf(" %12.4e",b[i]);} }

	// clear memory
	free(Ar); free(rhs);

	// QED :)
	return ret;
}

// function to minimise constraint violation using Newtons method (more efficient)
int CHakLevelSet::con_min3(int n, int m, int numCon, double *lam,  double *s, double *A, double *b,
			 double *lam_min, double *lam_max, double *up_lim, double *low_lim, int pinfo)
{
	int i,j,k,ind;

	// n = number boundary point
	// m = number of lambda values (1 + numCon)
	// A = constraint boundary integral coeffs
	// b = target constraint values

	// lots of data arrays
	double *dG_dlam = (double *) malloc(m*sizeof(double));
    //double dG_dlam[2];
	double *Hes_fd = (double *) malloc(m*m*sizeof(double));
    //double Hes_fd[4];
	double *dgdl_mat = (double *) malloc(numCon*m*sizeof(double)); // row is lam, col is g
	double *g_b = (double *) malloc(numCon*sizeof(double));
	double *g0 = (double *) malloc(numCon*sizeof(double));
	double *g1 = (double *) malloc(numCon*sizeof(double));
	double *g2 = (double *) malloc(numCon*sizeof(double));
	double *g3 = (double *) malloc(numCon*sizeof(double));
	double *g4 = (double *) malloc(numCon*sizeof(double));
	double *x = (double *) malloc(n*sizeof(double));
	double *lam_best = (double *) malloc(m*sizeof(double));

	double h_trial, step, ftemp;

	double h_lam = 1.0e+20;
	double h_best = 1.0e+20;
	int stop = 0;
	int count = 0;
	double inc = 1.0e-3; // step size
	double inc2 = 2.0*inc; // twice step size
	double incc = inc*inc; // step size squared
	do
	{
		// re-set gradient matrices
		for(i=0;i<m;i++){ dG_dlam[i] = 0.0; }
		k = m*m; for(i=0;i<k;i++){ Hes_fd[i] = 0.0; }

		// current constraint move values
		get_delD(n, m, x, lam, s, up_lim, low_lim);
		cblas_dgemv(CblasRowMajor, CblasNoTrans, numCon, n, 1.0, A, n, x, 1, 0.0, g0, 1);

		// current diff
		for(i=0;i<numCon;i++){ g_b[i] = (g0[i] - b[i]); }

		// compute sum of squared differences
		h_trial = 0.0;
		for(i=0;i<numCon;i++)
		{
			h_trial += (g_b[i]*g_b[i]) / (b[i]*b[i]); // normalize
		}
		h_trial = sqrt(h_trial);

		if(h_trial > h_best){step = 0.25;}
		else { step = 1.0; h_best = h_trial; for(i=0;i<m;i++){lam_best[i] = lam[i];} }
		h_lam = h_trial;
		if(h_lam < 1.0e-6){stop = 1;}
		if(pinfo==3){ printf("\nh_lam = %12.4e",h_lam); }

		if(stop == 0)
		{
			// compute dG / d_lam & FD Hessian
			for(i=0;i<m;i++) // for each lambda value
			{
				lam[i] += inc;
				get_delD(n, m, x, lam, s, up_lim, low_lim); // fwd diff
				cblas_dgemv(CblasRowMajor, CblasNoTrans, numCon, n, 1.0, A, n, x, 1, 0.0, g1, 1);

				lam[i] -= inc2;
				get_delD(n, m, x, lam, s, up_lim, low_lim); // back diff
				cblas_dgemv(CblasRowMajor, CblasNoTrans, numCon, n, 1.0, A, n, x, 1, 0.0, g2, 1);

				// central finite diff for 1st order gradient
				ind =(i*m)+i;
				k=(i*numCon);
				for(j=0;j<numCon;j++)
				{
					ftemp = (g1[j] - g2[j]) / inc2; // dg_j / dlam_i
					dG_dlam[i] += g_b[j] * ftemp;
					dgdl_mat[k++] = ftemp; // save for later (to complete Hessian)
				}

				lam[i] += inc; // re-set

				// compute diagonal Hessian entry
				for(j=0;j<numCon;j++)
				{
					Hes_fd[ind] += g_b[j] * (g1[j] + g2[j] - 2.0*g0[j]) / incc;
				}

				// off diagonal entries (symmetric) for row i
				for(k=i+1;k<m;k++)
				{
					lam[i] += inc;
					lam[k] += inc;
					get_delD(n, m, x, lam, s, up_lim, low_lim); // fwd-fwd diff
					cblas_dgemv(CblasRowMajor, CblasNoTrans, numCon, n, 1.0, A, n, x, 1, 0.0, g4, 1);

					lam[k] -= inc2;
					get_delD(n, m, x, lam, s, up_lim, low_lim); // fwd-back diff
					cblas_dgemv(CblasRowMajor, CblasNoTrans, numCon, n, 1.0, A, n, x, 1, 0.0, g3, 1);

					lam[i] -= inc2;
					get_delD(n, m, x, lam, s, up_lim, low_lim); // back-back diff
					cblas_dgemv(CblasRowMajor, CblasNoTrans, numCon, n, 1.0, A, n, x, 1, 0.0, g1, 1);

					lam[k] += inc2;
					get_delD(n, m, x, lam, s, up_lim, low_lim); // back-fwd diff
					cblas_dgemv(CblasRowMajor, CblasNoTrans, numCon, n, 1.0, A, n, x, 1, 0.0, g2, 1);

					lam[i] += inc;
					lam[k] -= inc; // re-set

					ind =(i*m)+k; // row i, col k
					for(j=0;j<numCon;j++)
					{
						// dg2_dj / (dlam_i dlam_k)
						Hes_fd[ind] += g_b[j] * (g4[j] - g3[j] - g2[j] + g1[j]) / (4.0*incc);
					}
				}
			}

			// complete Hessian
			for(i=0;i<m;i++)
			{
				for(j=i;j<m;j++)
				{
					ind =(i*m)+j; // row i, col j
					for(k=0;k<numCon;k++)
					{ Hes_fd[ind] += dgdl_mat[(i*numCon)+k] * dgdl_mat[(j*numCon)+k]; }
					if(j>i)
					{
						k = (j*m)+i;  // row j, col i
						Hes_fd[k] = Hes_fd[ind];
					}
				}
			}

			// check for zeros on diagonal
			/*for(i=0;i<m;i++)
			{
				ind =(i*m)+i; // row i, col i
				if(fabs(Hes_fd[ind]) < 1.0e-12 && fabs(dG_dlam[i]) < 1.0e-12){ Hes_fd[ind] = 1.0; }
			}*/

			// compute update step (will be in dG_dlam)
            //printf("\ndg_dlam: %12.4e , %12.4e",dG_dlam[0],dG_dlam[1]);
            //printf("\nHes: %12.4e , %12.4e , %12.4e , %12.4e", Hes_fd[0],Hes_fd[1],Hes_fd[2],Hes_fd[3]);
			i = getSolver().dsy_slvLPK(m, 1, Hes_fd, dG_dlam);
            if(i==-1)
            {
                printf("\nError returned from dsy_slvLPK");
                for(i=0;i<m;i++){dG_dlam[i] = 0.0;}
            }
			//dun_slvLPK('N', m, 1, Hes_fd, dG_dlam);

			stop=1;
			for(i=0;i<m;i++)
			{
				if(pinfo==3){ printf("\nlam_step[%i] = %12.4e",i,dG_dlam[i]); }
				if(fabs(dG_dlam[i]) > 1.0e-6){stop=0;}
			}
			for(i=0;i<m;i++)
			{
				lam[i] -= step*dG_dlam[i];
				if(lam[i] < lam_min[i]){lam[i] = lam_min[i];
					if(pinfo==3){printf("\nlam_min[%i]=%12.4e",i,lam_min[i]);} }
				else if(lam[i] > lam_max[i]){lam[i] = lam_max[i];
					if(pinfo==3){printf("\nlam_max[%i]=%12.4e",i,lam_max[i]);} }
				if(pinfo==3){printf("\nlam[%i]=%12.4e",i,lam[i]);}
			}
		}
		count++;

	} while(stop == 0 && count < 50);

	for(i=0;i<m;i++){lam[i] = lam_best[i];}

	free(dG_dlam);
	free(Hes_fd);
	free(dgdl_mat);
	free(g_b);
	free(g0);
	free(g1);
	free(g2);
	free(g3);
	free(g4);
	free(x);
	free(lam_best);

	return 0;
}

// sub-function to compute sum of constraint

// function to find an initial set of lambda values
void CHakLevelSet::get_lam0(int n, int m, int numCon, double *lam,  double *s, double *cA, double *b,
			  double *lam_min, double *lam_max)
{
	// first compute the matrix M = cA x s^T (m x m)
	double *M = (double *) malloc(m*m*sizeof(double));

	cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, m, m, n, 1.0, cA, n, s, n, 0.0, M, m);

	// invert M
	getSolver().din_LPK(M, m);

	int i,ind;
	double del_obj, ftemp;
	double *gc = (double *) malloc(m*sizeof(double));

	for(i=0;i<m;i++)
	{
		ind = i*m + 1; // start from 2nd column, row i
		gc[i] = cblas_ddot(numCon, &M[ind], 1, b, 1);
	}

	// use inverted matrix to find lowest possible value for obj change
	// loop though all lambda values
	for(i=0;i<m;i++)
	{
		ind = i*m; // 1st col, row i

		// if Min[0,i] > sub in lam_min[i], else sub in lam_max[i]
		if(M[ind] < 0.0){ ftemp = (lam_max[i] - gc[i]) / M[ind]; }
		else { ftemp = (lam_min[i] - gc[i]) / M[ind]; }

		// choose largest value for del_obj -> use with delCon to get initial lambda values
		if(i==0){ del_obj = ftemp; }
		else { del_obj = (ftemp > del_obj) ? ftemp : del_obj; }
	}

	for(i=0;i<m;i++)
	{
		ind = i*m; // 1st col, row i
		lam[i] = del_obj*M[ind] + gc[i];
	}

	free(M);
	free(gc);
}
