/*
   CHakMesh.cpp

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

#include "CHakMesh.h"
#include "CHakFiniteElement.h"
#include "CommonTypes.h"

//
// Constructor & Destructor
//
CHakMesh::CHakMesh() {
	// TODO Auto-generated constructor stub

}

CHakMesh::~CHakMesh() {
	// TODO Auto-generated destructor stub
}

//
// Implementation - Interfaces
//

/*
 Name:
 Description:
	-- This is originated from Numbering() function...
 	--
 Arguments:
 Return:
 Example:
*/



// Function to extract the structure from the lsf
// 1. determines the node and element status
// 2. discretizes the boundary
// 3. computes area ratios for AFG method
void CHakMesh::FindStruct(CHakLevelSet *m_levelset,CHakBoundary *m_boundary, double *alpha, double aMin)
{
	int n,m,o; //Incrementors
	int temp,temp2,num; // Temperary variables
	int count, count2, Acount; // More incrementors
	double ftemp;
	double lsf1, lsf2;
	double Ax, Ay;

	// read in data
	double h = this->m_lenEdge;
	double tol = this->m_tolerance;
	int elemX = this->m_elemX;
	int elemY = this->m_elemY;
	Elem **Number = this->m_pNumber;
	int NumElem = this->m_numElem;
	int NumNodes = this->m_numNodes;
	Coord *NodeCoord = this->m_pNodeCoord;
	double *lsf = m_levelset->m_pNodalLsf;
	Coord *AuxNodes = (Coord *) malloc(NumNodes * sizeof(Coord));
	Bseg *bptr = (Bseg *) malloc(NumElem * sizeof(Bseg));

	// Part 1: compute node & element status

	// temp arrays for node & element status
	short *NodeStat = (short *) malloc(NumNodes * sizeof(short));
	short *ElemStat = (short *) malloc(NumElem * sizeof(short));

	// Calculate node status
	for(n=0;n<NumNodes;n++)
	{
		// If signed distance function is 0.0 (within tollerance), then node lies on the boundary
		if(fabs(lsf[n]) < tol)	{
			NodeStat[n] = 2; // 2 = boundary node
		}
		// Otherwise if lsf = -ve node is OUT (0), +ve then node is IN (1)
		else	{
			NodeStat[n] = (lsf[n] > 0.0) ? 1 : 0;
		}
	}

	// Calculate Element status
	short sumO, sumB, sumI;	// Variables for summing IN, OUT and Boundary nodes
	int Lnodes[4];	// Array to temp store local node data

	// For all elements
	for(m=0;m<elemY;m++)
	{
		for(n=0;n<elemX;n++)
		{
			num = Number[n][m].n;	// Element number
			// Read in local node numbers
			Lnodes[0] = Number[n][m].a;
			Lnodes[1] = Number[n][m].b;
			Lnodes[2] = Number[n][m].c;
			Lnodes[3] = Number[n][m].d;

			//Set node sumuations back to zero
			sumO = 0; sumI = 0; sumB = 0;

			// Find number of OUT and Boundary nodes for the element
			for(o=0;o<4;o++)
			{
				sumO += (NodeStat[Lnodes[o]] == 0) ? 1 : 0;
				sumI += (NodeStat[Lnodes[o]] == 1) ? 1 : 0;
				sumB += (NodeStat[Lnodes[o]] == 2) ? 1 : 0;
			}

			if(sumO == 0)	{	// If No nodes are OUT then Element is IN
				ElemStat[num] = 10;
			}
			else if (sumI == 0)	{	/*If no nodes are IN then Element is OUT*/
				ElemStat[num] = 0;
			}
			else {
				ElemStat[num] = 1; /*Element status of an NIO element is initially 1*/
			}
		}
	}

	// Part 2. Discretize the boundary
	// Find & store all boundary sections and define Aux nodes
	int *na_count = (int *) calloc(NumNodes, sizeof(int)); // count number of aux nodes connected to an element
	int *na_conn = (int *) malloc(4*NumNodes * sizeof(int));
	int ncnt;
	count = 0;	// initialize count of auxiallry nodes
	count2 = 0; // initialize count of boundary segments
	int An[4]; // Array to store intersection node numbers

	// For all elements
	for(m=0;m<elemY;m++)
	{
		for(n=0;n<elemX;n++)
		{
			num = Number[n][m].n; //Element number
			if(ElemStat[num] != 0)
			{
				Acount = 0;
				// Read in local node numbers
				Lnodes[0] = Number[n][m].a;
				Lnodes[1] = Number[n][m].b;
				Lnodes[2] = Number[n][m].c;
				Lnodes[3] = Number[n][m].d;

				// Look at each edge in turn to determine if cut, or part of the boundary
				for(o=0;o<4;o++)
				{
					temp = (o == 3) ? 0 : (o + 1); // find next node round
					sumI = NodeStat[Lnodes[o]] + NodeStat[Lnodes[temp]]; // sum node status

					// If edge is cut (1 node is IN other is OUT) then update count
					if(sumI == 1)
					{
						// Compute co-ordintes of interection point (interpolation)
						lsf1 = lsf[Lnodes[o]];
						lsf2 = lsf[Lnodes[temp]];
						ftemp = (lsf1 + lsf2) / (lsf1 - lsf2);
						ftemp += 1.0;
						ftemp *= 0.5 * h;

						if( (ftemp > h) || (ftemp < 0.0) ) {
							printf("ERROR! Boundary distance = %lf, near node %i",ftemp,Lnodes[o]);
						}

						// co-ordinates of the aux node, depends on the cut edge
						switch(o) {
							case 0:
								Ax = NodeCoord[Lnodes[o]].x + ftemp;
								Ay = NodeCoord[Lnodes[o]].y;
								break;
							case 1:
								Ax = NodeCoord[Lnodes[o]].x;
								Ay = NodeCoord[Lnodes[o]].y + ftemp;
								break;
							case 2:
								Ax = NodeCoord[Lnodes[o]].x - ftemp;
								Ay = NodeCoord[Lnodes[o]].y;
								break;
							case 3:
								Ax = NodeCoord[Lnodes[o]].x;
								Ay = NodeCoord[Lnodes[o]].y - ftemp;
								break;
						}

						// Add the auxillary boundary node
						// need to check that the auxillary node doesn't already exist
						An[Acount] = -1;
						int end = na_count[Lnodes[o]];
						if(end > 0)
						{
							ncnt = 4*Lnodes[o];
							while(end > 0)
							{
								if( (fabs(Ax - AuxNodes[na_conn[ncnt]].x) < tol)
									&& (fabs(Ay - AuxNodes[na_conn[ncnt]].y) < tol) )
								{
									An[Acount] = na_conn[ncnt];
									end = 0;
								}
								ncnt++;
								end--;
							}
						}

						// If auxillary node does not exist, then create it
						if(An[Acount] == -1)
						{
							// add grid node -> aux node connectivity (for both nodes)
							ncnt = na_count[ Lnodes[o] ]; // number connected aux nodes so far
							ncnt += Lnodes[o]*4; // point to correct location
							na_conn[ncnt] = count; // add connectivity
							na_count[ Lnodes[o] ]++; // increase count

							ncnt = na_count[ Lnodes[temp] ]; // number connected aux nodes so far
							ncnt += Lnodes[temp]*4; // point to correct location
							na_conn[ncnt] = count; // add connectivity
							na_count[ Lnodes[temp] ]++;  // increase count

							AuxNodes[count].x = Ax;
							AuxNodes[count].y = Ay;
							An[Acount] = count;
							count++; // increase Aux node count
						}

						Acount++; //update count of boundary intersection points for this element
					}

					// If both nodes are on the boundary
					else if(sumI == 4)
					{
						// Add the elemet edge to boundary segment data
						bptr[count2].n1 = Lnodes[o];    // node number of point 1
						bptr[count2].n2 = Lnodes[temp]; // node number of point 2
						bptr[count2++].e = num;			// associated elememt number
					}
				}

				// For cut elements determine boundary segment(s)

				// if there are two edges cut, then a boundary segment must cross both
				if(Acount == 2)
				{
					bptr[count2].n1 = NumNodes + An[0];	// node number of point 1
					bptr[count2].n2 = NumNodes + An[1];	// node number of point 2
					bptr[count2++].e = num;				// associated elememt number
				}

				// if there is only one cut edge, then the boundary must also cross an element node
				else if(Acount == 1)
				{
					// find a node that is on the boundary & has an OUT neighbour
					for(o=0;o<4;o++)
					{
						// If node is on boudnary, check its neighbours
						if(NodeStat[Lnodes[o]] == 2)
						{
							temp = (o == 3) ? 0 : (o + 1);   // find next node round
							temp2 = (o == 0) ? 3 : (o - 1);  // find previous node

							// If a neighbour is an OUT node, then add boundary segment
							if( (NodeStat[Lnodes[temp]] == 0) || (NodeStat[Lnodes[temp2]] == 0) )
							{
								bptr[count2].n1 = NumNodes + An[0];	// node number of point 1
								bptr[count2].n2 = Lnodes[o];		// node number of point 2
								bptr[count2++].e = num;				// associated elememt number
							}
						}
					}
				}

				// if there are four cut edges, then determine which Aux node pairs form the boundary
				else if(Acount == 4)
				{
					// first determine lsf value at element centre - sign is only of interest
					ftemp = 0.0;
					for(o=0;o<4;o++) {
						ftemp += lsf[Lnodes[o]];
					}

					// Now look at status of node 1
					temp = NodeStat[Lnodes[0]];

					// Node pairs that correspond to boundary segments can easily be determined
					if( ( (temp == 1) && (ftemp > 0.0) ) || ( (temp == 0) && (ftemp < 0.0) ) )
					{
						bptr[count2].n1 = NumNodes + An[0];		// node number of point 1
						bptr[count2].n2 = NumNodes + An[1];		// node number of point 2
						bptr[count2++].e = num;					// associated elememt number

						bptr[count2].n1 = NumNodes + An[2];		// node number of point 1
						bptr[count2].n2 = NumNodes + An[3];		// node number of point 2
						bptr[count2++].e = num;					// associated elememt number
					}

					else
					{
						bptr[count2].n1 =  NumNodes + An[0];    // node number of point 1
						bptr[count2].n2 =  NumNodes + An[3];	// node number of point 2
						bptr[count2++].e = num;					// associated elememt number

						bptr[count2].n1 =  NumNodes + An[1];    // node number of point 1
						bptr[count2].n2 =  NumNodes + An[2];	// node number of point 2
						bptr[count2++].e = num;					// associated elememt number
					}

					// update element status to indicate if centre of element is IN or OUT
					//  4 = centre IN, 5 = centre OUT
					ElemStat[num] = (ftemp > 0.0) ? 4:5;

				}

				// if no edges are cut and element is not IN
				//  then boundary segment must cross diagonal
				else if( (Acount == 0) && (ElemStat[num] != 10) )
				{
					// find the two boudnary nodes
					for(o=0;o<4;o++)
					{
						if( NodeStat[Lnodes[o]] == 2)
						{
							An[Acount] = Lnodes[o];
							Acount++;
						}
					}
					bptr[count2].n1 = An[0];    // node number of point 1
					bptr[count2].n2 = An[1];	// node number of point 2
					bptr[count2++].e = num;		// associated elememt number
				}
			}
		}
	}

	// read back the totals for m_numAux & NumBound
	m_boundary->m_numAux = count;
	m_boundary->m_numBound = count2;

	// re-size auxillary node array and boundary segment data array to min size
    free(m_boundary->m_pAuxNodes);
	if(count>0){AuxNodes = (Coord *) realloc(AuxNodes, (count * sizeof(Coord)));}
    else{AuxNodes = (Coord *) realloc(AuxNodes, (sizeof(Coord)));}
    m_boundary->m_pAuxNodes = AuxNodes;

    free(m_boundary->Bound);
	if(count2>0){bptr = (Bseg *) realloc(bptr, (count2 * sizeof(Bseg)));}
    else{bptr = (Bseg *) realloc(bptr, (sizeof(Bseg)));}
    m_boundary->Bound = bptr;

	// consolodate na_conn data
	// each aux node between 2 grid nodes
    free(m_boundary->m_pConnAuxNode);
	if(count>0){m_boundary->m_pConnAuxNode = (int *) malloc(2*count*sizeof(int));}
    else{m_boundary->m_pConnAuxNode = (int *) malloc(sizeof(int));}
	count = 0;
	for(n=0;n<NumNodes;n++)
	{
		temp2 = na_count[n];
		m_boundary->m_pIndConnAuxNode[n]=count;
		if(temp2 > 0)
		{
			temp = 4*n; // point ot start of connected aux nodes
			while(temp2 > 0)
			{
				m_boundary->m_pConnAuxNode[count++] = na_conn[temp++];
				temp2--;
			}
		}
	}
	m_boundary->m_pIndConnAuxNode[NumNodes] = count; // end point
	free(na_count);
	free(na_conn);

	// Part 3. Compute area ratios for AFG method
	ComputeAreaRatio(alpha, NodeStat, ElemStat, m_boundary->m_pAuxNodes, count2, m_boundary->Bound, aMin);

	// clear memory
	free(NodeStat);
	free(ElemStat);
}

// function to compute area ratio for all elements
void CHakMesh::ComputeAreaRatio(double *alpha, short *NodeStat, short *ElemStat, Coord *AuxNodes, int NumBound, Bseg *Boundary, double aMin)
{
	int n,m;	// incrementors
	int temp,num; // temp variables
	double atemp;
	int Lnodes[4];	// Array to temp store local node data


	// read in data
	double h = this->m_lenEdge;
	double AreaElem = h * h; // element area
	int elemX = this->m_elemX;
	int elemY = this->m_elemY;
	Elem **Number = this->m_pNumber;
	int NumElem = this->m_numElem;
	int NumNodes = this->m_numNodes;
	Coord *NodeCoord = this->m_pNodeCoord;

	CHakFiniteElement fem;


	// For all elements
	for(m=0;m<elemY;m++)
	{
		for(n=0;n<elemX;n++)
		{
			num = Number[n][m].n;	// Element number
			// Read in local node numbers
			Lnodes[0] = Number[n][m].a;
			Lnodes[1] = Number[n][m].b;
			Lnodes[2] = Number[n][m].c;
			Lnodes[3] = Number[n][m].d;

			temp = ElemStat[num];

			// If element is IN
			if(temp == 10)
			{
				alpha[num] = 1.0;
			}
			// If element is OUT
			else if(temp == 0)
			{
				alpha[num] = aMin;
			}
			// If element is cut, but not status = 5

			else if(temp != 5)
			{
				atemp = fem.InArea(temp, num, Lnodes, NodeStat, NumNodes, NodeCoord, AuxNodes, NumBound, Boundary);
				alpha[num] = atemp / AreaElem;
				alpha[num] = (alpha[num] < aMin) ? aMin : alpha[num]; // enforce minimum
				//alpha[num] = (alpha[num] < 0.01) ? 0.01 : alpha[num];
			}
			// ElemStat=5 indicates that centre of element is outside structure
			else
			{
				atemp = AreaElem - fem.InArea(temp, num, Lnodes, NodeStat, NumNodes, NodeCoord, AuxNodes, NumBound, Boundary);
				alpha[num] = atemp / AreaElem;
				alpha[num] = (alpha[num] < aMin) ? aMin : alpha[num]; // enforce minimum
				//alpha[num] = (alpha[num] < 0.01) ? 0.01 : alpha[num];
			}
		}
	}
}


//
// Interfaces - OOD versions
//

// Function that numbers all elements and nodes in the FG domain
void CHakMesh::Numbering()
{
	int i, j;
	int num = 0; // first element number is zero

	// read data
	Elem **Number = this->m_pNumber; // pointer to Numbering array
	int m_elemX = this->m_elemX;
	int elemY = this->m_elemY;

	// Element Numbering loop
	for (j = 0; j < elemY; j++)
	{
		for (i = 0; i <m_elemX; i++)
		{
			Number[i][j].n = num++;
		}
	}

	// State first element node numbers
	num = 0; // first node number is zero
	Number[0][0].a = num;
	Number[0][0].b = ++num;
	Number[0][0].c = ++num;
	Number[0][0].d = ++num;

	// Number first X row (Y = 0)
	if (m_elemX > 1) // if there is more than one element in the x direction number first row
	{
		for (i = 1; i <m_elemX; i++)
		{
			Number[i][0].a = Number[i-1][0].b;
			Number[i][0].b = ++num;
			Number[i][0].c = ++num;
			Number[i][0].d = Number[i-1][0].c;
		}
	}

	// Number rest of FG mesh
	if (elemY > 1)
	{
		for (j = 1; j < elemY; j++)
		{
			for (i = 0; i < m_elemX; i++)
			{
				Number[i][j].a = Number[i][j-1].d;
				Number[i][j].b = Number[i][j-1].c;

				if (i == 0)
				{
					Number[i][j].c = ++num;
					Number[i][j].d = ++num;
				}
				else
				{
					Number[i][j].c = ++num;
					Number[i][j].d = Number[i-1][j].c;
				}
			}
		}
	}
}

// Node Co-ordinate calculation function
void CHakMesh::Coordinates()
{
	// data from inMesh
	int m_elemX = this->m_elemX;
	int elemY = this->m_elemY;
	double hx = this->m_lenEdge; 
	double hy = this->m_lenEdge;
	double tol = this->m_tolerance;
	Elem **Number = this->m_pNumber; // pointer to Numbering array
	Coord *NodeCoord = this->m_pNodeCoord; // pointer to node coords

	int i, j, num;

	for (j = 0; j < elemY; j++)
	{
		for (i = 0; i < m_elemX; i++)
		{
			num = Number[i][j].a;			
			if((NodeCoord[num].x + NodeCoord[num].y) < tol)
			{
				// If node co-ordinate hasn't been calculated already
				NodeCoord[num].x = i * hx;
				NodeCoord[num].y = j * hy;
			}
			
			num = Number[i][j].b;
			if((NodeCoord[num].x + NodeCoord[num].y) < tol)
			{
				NodeCoord[num].x = (i * hx) + hx;
				NodeCoord[num].y = j * hy;
			}

			num = Number[i][j].c;
			if((NodeCoord[num].x + NodeCoord[num].y) < tol)
			{
				NodeCoord[num].x = (i * hx) + hx;
				NodeCoord[num].y = (j * hy) + hy;
			}
			
			num = Number[i][j].d;
			if((NodeCoord[num].x + NodeCoord[num].y) < tol)
			{
				NodeCoord[num].x = i * hx;
				NodeCoord[num].y = (j * hy) + hy;
			}
		}
	}
}

// function that orders node numbers into a 2D based on their relative positions
void CHakMesh::NodeNums2D()
{
	// data from inMesh
	int NodeX = this->m_nodeX;
	int NodeY = this->m_nodeY;
	int NumNodes = this->m_numNodes;
	double h = this->m_lenEdge;
	double tol = this->m_tolerance;
	int **Nodes2 = this->m_pNodes2D; // pointer to Node numbering array
	Coord *NodeCoord = this->m_pNodeCoord; // pointer to node coords

	int i; // incrementor
	int X, Y; // node position variables

	// for all grid node determine position in ordered 2D array
	//  by their co-ordinates and element edge length
	for (i = 0; i < NumNodes; i++)
	{
		X = (int)(floor((NodeCoord[i].x / h)+tol) + 1);
		Y = (int)(floor((NodeCoord[i].y / h)+tol) + 1); // +1 for layer of ghost nodes
		Nodes2[X][Y] = i;
	}

	// now copy node numbers to create a layer of ghost nodes around the grid

	// fill in the four corners
	Nodes2[0][0] = Nodes2[1][1];
	Nodes2[0][NodeY-1] = Nodes2[1][NodeY-2];
	Nodes2[NodeX-1][0] = Nodes2[NodeX-2][1];
	Nodes2[NodeX-1][NodeY-1] = Nodes2[NodeX-2][NodeY-2];

	// now fill in top and bottom rows
	for (i = 1; i < NodeX - 1; i++)
	{
		Nodes2[i][0] = Nodes2[i][1];
		Nodes2[i][NodeY-1] = Nodes2[i][NodeY-2];
	}

	// Finally fill in left and right columns
	for (i = 1; i < NodeY - 1; i++)
	{
		Nodes2[0][i] = Nodes2[1][i];
		Nodes2[NodeX-1][i] = Nodes2[NodeX-2][i];
	}
}

// function to number bars elements
void CHakMesh::NumberingBarElements()
{
	// read data
	Elem **pNumber = this->m_pNumber; // pointer to Numbering array
	int m_elemX = this->m_elemX;
	int elemY = this->m_elemY;

	int i,j;
	int xend = m_elemX - 1;
	int yend = elemY - 1;
	int count = 0;

	// loop through all elements
	for (j = 0; j < elemY; j++)
	{
		for (i = 0; i < m_elemX; i++)
		{
			// use bottom edge as x-bars
			this->m_pBarNums[count].e = count; // element number (bar elems only)
			this->m_pBarNums[count].n1 = pNumber[i][j].a * 2; // x-dof of node 1
			this->m_pBarNums[count++].n2 = pNumber[i][j].b * 2; // x-dof of node 2

			// also use top edge for top row of elements
			if (j == yend)
			{
				this->m_pBarNums[count].e = count; // element number (bar elems only)
				this->m_pBarNums[count].n1 = pNumber[i][j].d * 2; // x-dof of node 1
				this->m_pBarNums[count++].n2 = pNumber[i][j].c * 2; // x-dof of node 2
			}
		}
	}

	// loop through all elements
	for (i = 0;i < m_elemX; i++)
	{
		for (j = 0;j < elemY; j++)
		{
			// use left edge as y-bars
			this->m_pBarNums[count].e = count; // element number (bar elems only)
			this->m_pBarNums[count].n1 = (pNumber[i][j].a * 2) +1; // y-dof of node 1
			this->m_pBarNums[count++].n2 = (pNumber[i][j].d * 2) +1; // y-dof of node 2

			// also use right edge for last column of elements
			if(i == xend)
			{
				this->m_pBarNums[count].e = count; // element number (bar elems only)
				this->m_pBarNums[count].n1 = (pNumber[i][j].b * 2) +1; // y-dof of node 1
				this->m_pBarNums[count++].n2 = (pNumber[i][j].c * 2) +1; // y-dof of node 2
			}
		}
	}

	// that is all!!
}

// function to find the nearest grid node number to a set of co-ordinates
int CHakMesh::CloseNode(double xp, double yp)
{
	// read data
	int NumNodes = this->m_numNodes;
	Coord *NodeCoord = this->m_pNodeCoord;
	double tol = this->m_tolerance;

	int n;
	int nodenum=-1;
	double dtemp, d2, dx, dy;	// variables to keep track of distnace closest node
	double xn, yn1;				// variables for current node co-ordinates
	dtemp =  (2.0 * this->m_lenEdge);
	dtemp *= dtemp;			// set initial distance to twice an elements edge length squared

	for (n = 0; n < NumNodes; n++) // For all nodes
	{
		xn = NodeCoord[n].x;
		yn1 = NodeCoord[n].y; // read in node co-ords

		dx = xp - xn;
		dy = yp - yn1; // find difference between node and target co-ords
		d2 = (dx * dx) + (dy * dy);	// Squared distance between node and point

		if (fabs(d2) < tol)
		{	// if co-ordinates are same as current node, that that is closest
			return(n);
		}
		if (d2 < dtemp)
		{	// if distance of current node is less than previous, then update
			nodenum = n;
			dtemp = d2;
		}
	}

	return(nodenum);
}
