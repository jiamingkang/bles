/*
 *  ABFG.c
 *
 *  Created by Peter Dunning
 *	Version 5.4 @ 23/04/2014
 *  Functions related to the area-ratio weighted fixed grid method
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "ABFG.h"
#include "EMatrix.h"
#include "Solve.h"
#include "Sens.h"

// function to find the nearest grid node number to a set of co-ordinates
int closeNode(mesh *inMesh, double xp, double yp)
{
	// read data
	int NumNodes = inMesh->NumNodes;
	Coord *NodeCoord = inMesh->NodeCoord;
	double tol = inMesh->tol;
	
	int n;
	int nodenum=-1;
	double dtemp, d2, dx, dy;	// variables to keep track of distnace closest node
	double xn, yn1;				// variables for current node co-ordinates
	dtemp =  (2.0 * inMesh->h);
	dtemp *= dtemp;			// set initial distance to twice an elements edge length squared
	
	for(n=0;n<NumNodes;n++) // For all nodes
	{
		xn = NodeCoord[n].x;
		yn1 = NodeCoord[n].y; // read in node co-ords
		
		dx = xp - xn;
		dy = yp - yn1; // find difference between node and target co-ords
		d2 = (dx * dx) + (dy * dy);	// Squared distance between node and point
		
		if(fabs(d2) < tol)	{	// if co-ordinates are same as current node, that that is closest
			return(n);
		}
		if(d2 < dtemp)	{	// if distance of current node is less than previous, then update
			nodenum = n;
			dtemp = d2;
		}
	}
	return(nodenum);
}

// Function to extract the structure from the lsf
// 1. determines the node and element status
// 2. discretizes the boundary
// 3. computes area ratios for AFG method
void Find_struct(mesh *inMesh, levSet *levelset, boundary *bound_in, double *alpha, double aMin)
{
	int n,m,o; //Incrementors
	int temp,temp2,num; // Temperary variables
	int count, count2, Acount; // More incrementors
	double ftemp;
	double lsf1, lsf2;
	double Ax, Ay;
	
	// read in data
	double h = inMesh->h;
	double tol = inMesh->tol;
	int elemX = inMesh->elemX;
	int elemY = inMesh->elemY;
	Elem **Number = inMesh->Number;
	int NumElem = inMesh->NumElem;
	int NumNodes = inMesh->NumNodes;
	Coord *NodeCoord = inMesh->NodeCoord;
	double *lsf = levelset->lsf;
	Coord *AuxNodes = malloc(NumNodes * sizeof(Coord));
	Bseg *bptr = malloc(NumElem * sizeof(Bseg));
	
	// Part 1: compute node & element status
	
	// temp arrays for node & element status
	short *NodeStat = malloc(NumNodes * sizeof(short));
	short *ElemStat = malloc(NumElem * sizeof(short));
	
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
	int *na_count = calloc(NumNodes, sizeof(int)); // count number of aux nodes connected to an element
	int *na_conn = malloc(4*NumNodes * sizeof(int));
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
	
	// read back the totals for NumAux & NumBound
	bound_in->NumAux = count;
	bound_in->NumBound = count2;
	
	// re-size auxillary node array and boundary segment data array to min size
    free(bound_in->AuxNodes);
	if(count>0){AuxNodes = realloc(AuxNodes, (count * sizeof(Coord)));}
    else{AuxNodes = realloc(AuxNodes, (sizeof(Coord)));}
    bound_in->AuxNodes = AuxNodes;
    
    free(bound_in->Bound);
	if(count2>0){bptr = realloc(bptr, (count2 * sizeof(Bseg)));}
    else{bptr = realloc(bptr, (sizeof(Bseg)));}
    bound_in->Bound = bptr;
	
	// consolodate na_conn data
	// each aux node between 2 grid nodes
    free(bound_in->na_conn);
	if(count>0){bound_in->na_conn = malloc(2*count*sizeof(int));}
    else{bound_in->na_conn = malloc(sizeof(int));}
	count = 0;
	for(n=0;n<NumNodes;n++)
	{
		temp2 = na_count[n];
		bound_in->na_conn_ind[n]=count;
		if(temp2 > 0)
		{
			temp = 4*n; // point ot start of connected aux nodes
			while(temp2 > 0)
			{
				bound_in->na_conn[count++] = na_conn[temp++];
				temp2--;
			}
		}
	}
	bound_in->na_conn_ind[NumNodes] = count; // end point
	free(na_count);
	free(na_conn);
	
	// Part 3. Compute area ratios for AFG method
	AFG_area(inMesh, alpha, NodeStat, ElemStat, bound_in->AuxNodes, count2, bound_in->Bound, aMin);
	
	// clear memory
	free(NodeStat);
	free(ElemStat);
}

// function to compute area ratio for all elements
void AFG_area(mesh *inMesh, double *alpha, short *NodeStat, short *ElemStat, Coord *AuxNodes,
				int NumBound, Bseg *Boundary, double aMin)
{
	int n,m;	// incrementors
	int temp,num; // temp variables
	double atemp;
	int Lnodes[4];	// Array to temp store local node data
	
	// read in data
	double h = inMesh->h;
	double AreaElem = h * h; // element area
	int elemX = inMesh->elemX;
	int elemY = inMesh->elemY;
	Elem **Number = inMesh->Number;
	int NumNodes = inMesh->NumNodes;
	Coord *NodeCoord = inMesh->NodeCoord;
	
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
				atemp = InArea(temp, num, Lnodes, NodeStat, NumNodes, NodeCoord, AuxNodes, NumBound, Boundary);
				alpha[num] = atemp / AreaElem;
				alpha[num] = (alpha[num] < aMin) ? aMin : alpha[num]; // enforce minimum
				//alpha[num] = (alpha[num] < 0.01) ? 0.01 : alpha[num];
			}
			// ElemStat=5 indicates that centre of element is outside structure
			else
			{
				atemp = AreaElem - InArea(temp, num, Lnodes, NodeStat, NumNodes, NodeCoord, AuxNodes, NumBound, Boundary);
				alpha[num] = atemp / AreaElem;
				alpha[num] = (alpha[num] < aMin) ? aMin : alpha[num]; // enforce minimum
				//alpha[num] = (alpha[num] < 0.01) ? 0.01 : alpha[num];
			}
		}
	}
}

// function to compute a cut element area
double InArea(int eStat, int eNum, int *Lnodes, short *NodeStat, int NumNodes, Coord *NodeCoord,
			  Coord *AuxNodes, int NumBound, Bseg *Boundary)
{
	int numPts = 0;
	int i, j, k, temp, temp2;
	Coord pts[6]; // Max number of points will be 6
	
	// need to change the algorithm if element status = 5
	short nLook = (eStat==5) ? 0:1; // compute area outside structure if ElemStat is 5
	
	// first look at all elemet nodes
	for(i=0;i<4;i++)
	{
		// if the node is IN (or OUT for ElemStat=5) add it to the pts array
		if( NodeStat[Lnodes[i]] == nLook )
		{
			pts[numPts].x = NodeCoord[Lnodes[i]].x;
			pts[numPts++].y = NodeCoord[Lnodes[i]].y;
		}
	
		// if node on the boundary, check it is part of the area polygon, but not a boundary segment
		else if( NodeStat[Lnodes[i]] == 2 )
		{
			// check both neighbouring nodes
			temp = (i == 3) ? 0 : (i + 1);   // find next node round
			temp2 = (i == 0) ? 3 : (i - 1);  // find previous node
			
			// if both neighbours are IN then add to pts array
			if( (NodeStat[Lnodes[temp]] == 1) && (NodeStat[Lnodes[temp2]] == 1) )
			{
				pts[numPts].x = NodeCoord[Lnodes[i]].x;
				pts[numPts++].y = NodeCoord[Lnodes[i]].y;
			}
		}
	}
	
	// now find and add boundary segment nodes
	temp2 = NumBound;
	for(i=0;i<temp2;i++)
	{
		if(Boundary[i].e == eNum)
		{
			if(Boundary[i].n1 >= NumNodes) {
				temp = Boundary[i].n1 - NumNodes;
				pts[numPts].x = AuxNodes[temp].x;
				pts[numPts++].y = AuxNodes[temp].y;
			}
			else {
				temp = Boundary[i].n1;
				pts[numPts].x = NodeCoord[temp].x;
				pts[numPts++].y = NodeCoord[temp].y;
			}
			
			if(Boundary[i].n2 >= NumNodes) {
				temp = Boundary[i].n2 - NumNodes;
				pts[numPts].x = AuxNodes[temp].x;
				pts[numPts++].y = AuxNodes[temp].y;
			}
			else {
				temp = Boundary[i].n2;
				pts[numPts].x = NodeCoord[temp].x;
				pts[numPts++].y = NodeCoord[temp].y;
			}
			
			// once a segment is found, only need to look at next one
			temp2 = (temp2 == NumBound) ? i+2:temp2; 
			temp2 = (temp2 > NumBound) ? NumBound:temp2;
		}
	}
	
	if(numPts < 3) {
		printf("ERROR! Only %i polygon points found for Element %i",numPts,eNum);
	}
	
	// Need to check if edges of the polygon cross - and untangle
	Coord ctemp; // variable to swap points
	Coord chkPts[4]; // Array to send data to LineCross function
	if(numPts > 3) // if only 3 points, then segments can't cross
	{
		do {
			temp = 0;
			for(i=0;i<(numPts-2);i++)
			{
				chkPts[0].x = pts[i].x;
				chkPts[0].y = pts[i].y;
				chkPts[1].x = pts[i+1].x;
				chkPts[1].y = pts[i+1].y; // current segment
				
				temp2 = (i==0) ? (numPts-1) : numPts;
				for(j=i+2;j<temp2;j++)
				{
					chkPts[2].x = pts[j].x;
					chkPts[2].y = pts[j].y;
					
					k = (j == (numPts-1)) ? 0:(j+1); // next point round
					
					chkPts[3].x = pts[k].x;
					chkPts[3].y = pts[k].y; // comparrison segment
					
					if(LineCross(chkPts) == 1)
					{
						// swap end point of 1st seg with start point of 2nd seg
						ctemp.x = pts[i+1].x;
						ctemp.y = pts[i+1].y;
						pts[i+1].x = pts[j].x;
						pts[i+1].y = pts[j].y;
						pts[j].x = ctemp.x;
						pts[j].y = ctemp.y;
						
						temp = 1; // indicate a line cross was found
						j=numPts;
						i=numPts; // break out of the for loops
					}
				}
			}
		} while(temp == 1); // Until no line crosses are found
	}
	
	// Use PolyArea to compute element area - return value
	return(PolyArea(numPts, pts));
}

// Function that calculates the area of any Polygon
// NB: vertices have to be numbered anti-clockwise
double PolyArea(int N,Coord *point)	
{
	int i,j;
	double area = 0.0;
	
	for (i=0;i<N;i++) // For all verticies of the polygon
	{
		j = (i == N-1) ? 0 : (i+1);	// Need to loop back to first point to complete algorithm
		area += point[i].x * point[j].y;
		area -= point[j].x * point[i].y;
	}
	
	area *= 0.5;
	return(fabs(area)); // return absolute value just in case
}

// function to determine if two lines (p1 -> p2 & p3 -> p4) cross
short LineCross(Coord *pts)
{
	int i;
	double tol = 1.0e-5;
	
	double dx1 = pts[0].x - pts[1].x;
	double dx3 = pts[2].x - pts[3].x;
	double dy1 = pts[0].y - pts[1].y;
	double dy3 = pts[2].y - pts[3].y;
	
	double DetA = det2(pts[0].x, pts[0].y, pts[1].x, pts[1].y);
	double DetC = det2(pts[2].x, pts[2].y, pts[3].x, pts[3].y);
	double DetDiv = det2(dx1, dy1, dx3, dy3);
	
	// If DetDiv is zero, then lines are parallel*/
	if(fabs(DetDiv) < tol) { return(0); }
	
	// compute intercept co-ordinates
	Coord Inter;
	Inter.x = det2(DetA, dx1, DetC, dx3) / DetDiv;
	Inter.y = det2(DetA, dy1, DetC, dy3) / DetDiv;
	
	// compute max and min co-ordinates that contain the lines
	Coord max, min;
	max.x = pts[0].x;
	min.x = pts[0].x;
	max.y = pts[0].y;
	min.y = pts[0].y; // initialize
	

	for(i=1;i<4;i++)
	{
		max.x = (pts[i].x > max.x) ? pts[i].x : max.x;
		min.x = (pts[i].x < min.x) ? pts[i].x : min.x;
		max.y = (pts[i].y > max.y) ? pts[i].y : max.y;
		min.y = (pts[i].y < min.y) ? pts[i].y : min.y; // update
	}

	// shrink bounding box, in-case lines cross at element corner
	max.x -= tol;
	max.y -= tol;
	min.x += tol;
	min.y += tol;
	
	// if intersection point is outside limits, then lines do not cross
	if( (Inter.x > max.x) || (Inter.x < min.x) || (Inter.y > max.y) || (Inter.y < min.y) )
	{ return(0); }
	
	// Otherwise the lines DO cross
	return(1);
}

// compute determinant of a 2x2 matrix
double det2(double a, double b, double c, double d)
{
	return( (a * d) - (b * c) );
}

// Assembles global stiffness (& maybe mass) matrix for AFG method in triplet format (for MA57 solver)
void AFG_Matrix(int mass, double **KE, double **ME, sp_mat *Kg, sp_mat *Mg, sp_mat *lump_mass, double *alpha,
					mesh *inMesh, isoMat *inMat, double aMin, double mMin)
{
	int n,m,o,num;
	int tnodes[4];	// temp array for node numbers
	int K_begin = 0;  // Begin variable used to track how many entries have been indexed so far
	int M_begin = 0;
	int mat_count = 0; // count for designable material
	
	double KEQuad[KE_SIZE]; // 2D Array for an Quadrilateral Stiffness Matrix
	double MEQuad[KE_SIZE]; // 2D Array for an Quadrilateral Mass Matrix
	double atemp; // area ratio variable
	double mtemp,e_fac,p_fac,e_rat,p_rat; // variables for designable material modification
	// material ratios
	if(inMesh->des_mat)
	{
		e_rat = inMat[inMesh->mat2].e / inMat[inMesh->mat1].e;
		p_rat = inMat[inMesh->mat2].rho / inMat[inMesh->mat1].rho;
	}
	
	// read data
	int elemX = inMesh->elemX;
	int elemY = inMesh->elemY;
	Elem **Number = inMesh->Number;
	double *kptr, *mptr; // pointer to correct base matrices
	
	// For All Elements
	for(m=0;m<elemY;m++)
	{	
		for(n=0;n<elemX;n++)
		{
			num = Number[n][m].n; // Element number
			atemp = alpha[num]; // Element area ratio

			if(atemp > 1.0e-12) // If element has area ratio > 0: calculate area weighted matrix
			{
				if(atemp > 1.0) {
					printf("\nERROR! alpha for element %i = %lf!",num,alpha[num]);
				}
				
				o = inMesh->mat_type[num]; // material number
				kptr = KE[o]; mptr = ME[o]; // point to correct matrices for elem material
				
				// Read in grid node numbers of current element
				tnodes[0] = Number[n][m].a;
				tnodes[1] = Number[n][m].b;
				tnodes[2] = Number[n][m].c;
				tnodes[3] = Number[n][m].d;
				
				// If element is IN, no need to multiply by alpha (unless using designable material)
				if(!inMesh->des_mat && atemp > 0.9999)
				{
					// Assemble area ratio weighted stiffness matrix into index arrays
					Assemble2(&K_begin, &M_begin, 4, NUM_DOF, tnodes, kptr, mptr, Kg, Mg, mass);
				}
				
				else
				{
					// modify matrices for designable materials
					if(inMesh->des_mat && num==inMesh->mat_elems[mat_count])
					{
						// material mix variable
						mtemp = inMesh->mat_vars[mat_count++]; 
						
						// modulus factor
						if(inMesh->mat_lin){e_fac = mtemp*e_rat + (1.0-mtemp);}
                        else{e_fac = HS_mat(mtemp, 0.5, &inMat[inMesh->mat1], &inMat[inMesh->mat2]) / inMat[inMesh->mat1].e;}
                        e_fac *= atemp;
						
						// density factor
                        p_fac =  mtemp*p_rat + (1.0-mtemp);
						p_fac *= (atemp == aMin) ? mMin : atemp;
					}
					else
					{
						e_fac = atemp;
						p_fac = (atemp == aMin) ? mMin : atemp;
					}
					
					// multiply In stiffness matrix by area ratio
					for(o=0;o<KE_SIZE;o++) {
						KEQuad[o] = kptr[o] * e_fac;
					}
					// multiply mass matrix by area ratio
					if(mass==1) {
						for(o=0;o<KE_SIZE;o++) {
							MEQuad[o] = mptr[o] * p_fac;
						}
					}
					
					// Assemble area ratio weighted stiffness matrix into index arrays
					Assemble2(&K_begin, &M_begin, 4, NUM_DOF, tnodes, KEQuad, MEQuad, Kg, Mg, mass);
				}
			}
            else if(inMesh->des_mat && num==inMesh->mat_elems[mat_count]){mat_count++;}
		}
	}
	
	// Add lumped masses
	m = lump_mass->ne;
	if(m>0)
	{
		for(n=0;n<m;n++)
		{
			Mg->irn[M_begin] = lump_mass->irn[n];
			Mg->jcn[M_begin] = lump_mass->jcn[n];
			Mg->A[M_begin++] = lump_mass->A[n];
		}
	}
	
	// Add bar elements
	if(inMesh->bars)
	{
		int numBar = inMesh->NumBars;
		double EoL = inMesh->bar_mat->e / inMesh->h; // E / L
		double kbar;
		
		for(n=0;n<numBar;n++)
		{
			m = inMesh->bar_nums[n].n1;
			o = inMesh->bar_nums[n].n2;
			kbar = EoL * inMesh->bar_areas[n]; // bar stiffness
			
			// swap if necessary (m < o)
			if(m > o){num=o; o=m; m=num;}
			
			// 3 entries
			Kg->irn[K_begin] = m; Kg->jcn[K_begin] = m; Kg->A[K_begin++] = kbar;
			Kg->irn[K_begin] = o; Kg->jcn[K_begin] = o; Kg->A[K_begin++] = kbar;
			Kg->irn[K_begin] = m; Kg->jcn[K_begin] = o; Kg->A[K_begin++] = -kbar;
		}
	}
	
	// Add designable bc stiffness
	if(inMesh->des_bc)
	{
		int numBC = inMesh->NumBC;
		
		for(o=0;o<numBC;o++)
		{
			// compute element indices
			m = floor((double)inMesh->BC_nums[o] / (double)elemX);
			n = inMesh->BC_nums[o] - m*elemX;
			
			// Read in grid node numbers of current element
			tnodes[0] = Number[n][m].a;
			tnodes[1] = Number[n][m].b;
			tnodes[2] = Number[n][m].c;
			tnodes[3] = Number[n][m].d;
			
			for(n=0;n<4;n++) // each node
			{
				num = tnodes[n]*NUM_DOF; // start of nodal dof
				for(m=0;m<NUM_DOF;m++)   // each dof
				{
					Kg->irn[K_begin] = num+m; 
					Kg->jcn[K_begin] = Kg->irn[K_begin]; 
					Kg->A[K_begin++] = pow(inMesh->K_bc[o],3.0)*inMesh->K0_bc; // penalize design variable
				}
			}
		}
	}
	
	// resize matrices
	Kg->ne = K_begin;
	Mg->ne = M_begin;
}

// function to compute gauss point coords
void Gauss_Coord(mesh *inMesh, Coord *gCoord)
{
	// read in mesh data
	double h = inMesh->h;
	int elemX = inMesh->elemX;
	int elemY = inMesh->elemY;
	Elem **Number = inMesh->Number;
	Coord *NodeCoord = inMesh->NodeCoord;
	double h2 = 0.5 * h;
	
	int n,m,i,num,node1;
	double cx,cy;
	
	// Array storing data for location of the 4 nodes in the square element
	static Coord Po[4] = {{-1.,-1.},{1.,-1.},{1.,1.},{-1.,1.}};
	double ga = 0.57735026919;
	
	// For All Elements
	for(m=0;m<elemY;m++)
	{	
		for(n=0;n<elemX;n++)
		{	
			num = Number[n][m].n * 4;
			node1 = Number[n][m].a;
			
			// work out global co-ordinates of element center
			// from co-ords of node 1 and element edge length
			cx = NodeCoord[node1].x + h2;
			cy = NodeCoord[node1].y + h2;
			
			// compute global gauss point co-ordinates and store
			for(i=0;i<4;i++)
			{
				gCoord[num + i].x = cx + (ga * Po[i].x * h2);
				gCoord[num + i].y = cy + (ga * Po[i].y * h2);
			}
		}
	}	
}

// calculate sensitivies using least squares of integration points for AFG method
void AFG_Sens(mesh *inMesh, boundary *bound_in, double *alpha, isoMat *inMat,  double *Nsens,
				double **prim, double **dual, int numDual, int numCase, double *wgt, Coord *gCoord, 
					double aMin, int mode, double *fact, bool sw, Coord *acc)
{
	// read in mesh data
	double h = inMesh->h;
	double thk = inMesh->t;
	int NumNodes = inMesh->NumNodes;
	int elemX = inMesh->elemX;
	int elemY = inMesh->elemY;
	int NumElem = inMesh->NumElem;
	Elem **Number = inMesh->Number;
	Coord *NodeCoord = inMesh->NodeCoord;
	double rad = 2.0 * h; // hard coded for two elements around a node

	int irad = 2; //ceil(rad/h); // search space
	rad *=rad; // input squared dist to Lsens
	Coord *AuxNodes = bound_in->AuxNodes;
	int numBound = bound_in->NumBound;
	Bseg *Boundary = bound_in->Bound;
	int Ntot = NumNodes + bound_in->NumAux;
	double alMin = (aMin < 1.0e-6) ? 1.0e-6 : aMin; // numerical tollerance for sens calc
	if(mode==1){ alMin = (aMin < 1.0e-2) ? 1.0e-2 : aMin; } // increase for eigenvalue sensitivities
	else if(mode==2){alMin = 0.0;} // for compliant mechanisms
	
	int i,n,m,o,num,ind,node,ex,ey,xMin,xMax,yMin,yMax;
	double atemp;
	//if(mode==2){w1 = -1.0/wgt[0]; w2 = wgt[1]*w1*w1;} // extra weights for complinat mech disp constraint
	int Gcount = 0;		// varible to track number of points evaluated
	int tnodes[4];
	int gpoints = NumElem * 4; // number of sensitivity values (4 per element)
	double *gSens = calloc(gpoints*numDual,sizeof(double));	// define memory for gauss point sensitivities
	isoMat *mptr; // pointer for material
	
	int mat_count = 0; // count for designable material
	double mtemp; // variables for designable material modification
	isoMat fgm;
	
	// Step 1. Compute gauss point senstivities
	
	// For All Elements
	for(m=0;m<elemY;m++)
	{	
		for(n=0;n<elemX;n++)
		{	
			num = Number[n][m].n; // Element number
			atemp = alpha[num];
			if(atemp > alMin) // If Element isn't OUT
			{	
				tnodes[0] = Number[n][m].a;
				tnodes[1] = Number[n][m].b;
				tnodes[2] = Number[n][m].c;
				tnodes[3] = Number[n][m].d;
				
				o = inMesh->mat_type[num]; // material number
				mptr = &inMat[o]; // point to correct material
                
                // modify matrices for designable materials
                if(inMesh->des_mat && num==inMesh->mat_elems[mat_count])
                {
                    // material mix variable
                    mtemp = inMesh->mat_vars[mat_count++];
                    
                    // modulus
                    if(inMesh->mat_lin){ fgm.e =  (1.0-mtemp)*inMat[inMesh->mat1].e + mtemp*inMat[inMesh->mat2].e; }
                    else{ fgm.e = HS_mat(mtemp, 0.5, &inMat[inMesh->mat1], &inMat[inMesh->mat2]); }
                    
                    // density
                    fgm.rho =  mtemp*inMat[inMesh->mat2].rho + (1.0-mtemp)*inMat[inMesh->mat1].rho;
                    
                    // material prop matrix
                    mtemp = fgm.e / inMat[inMesh->mat1].e;
                    for(i=0;i<9;i++) {
                        //fgm.mat[i] = mtemp*inMat[inMesh->mat1].mat[i] + (1.0-mtemp)*inMat[inMesh->mat2].mat[i];
                        fgm.mat[i] = mtemp*inMat[inMesh->mat1].mat[i];
                    }
                    
                    fgm.v = inMat[inMesh->mat1].v;
                    mptr = &fgm;
                }
                else
                {
                    mptr = &inMat[o]; // point to correct material (not FGM)
                }
				
				//calcualte sensitivity at Gauss points
				if(mode==1) // eigenvalue sensitivites
				{
					GaEigSens_Q4(tnodes, prim, dual, atemp, h, thk, mptr, Gcount, gSens, numCase, wgt);
				}
				else if(mode==2) // compliant mechanism sensitivites
				{
					double *gSens2 = calloc(8, sizeof(double));
					double *gSens3 = calloc(4, sizeof(double));	// temp memory for gauss point sensitivities
					
					GaSens_Q4(tnodes, prim, dual, atemp, h, thk, mptr, 0, gSens2, 1, 2, wgt, sw, acc);
					GaSens_Q4(tnodes, &dual[1], &dual[1], atemp, h, thk, mptr, 0, gSens3, 1, 1, wgt, sw, acc);
					
					// process sensitivities for compliant mech
					for(i=0;i<4;i++)
					{		
						ind = (Gcount+i)*3; // point to place in gSens
						gSens[ind] = -fact[0]*gSens2[(2*i)+1] - fact[1]*gSens3[i];
						gSens[ind+1] = gSens2[2*i]*fact[4] + fact[2]*gSens2[(2*i)+1] + fact[3]*gSens3[i];
						gSens[ind+2] = gSens3[i]*fact[5];
					}
					
					free(gSens2);
					free(gSens3);
				}
				else
				{ GaSens_Q4(tnodes, prim, dual, atemp, h, thk, mptr, Gcount, gSens, numCase, numDual, wgt, sw, acc); }
			}
            else if(inMesh->des_mat && num==inMesh->mat_elems[mat_count]){mat_count++;} // keep counting
			Gcount += 4; // update point count
		}
	}
	
	// Step 2. Compute node sensitivities using least squares method
	
	double *sens_temp = malloc(numDual * sizeof(double)); // array for sens of each dual state
	
	// calculate smoothed sensitivities for all boundary (or non-OUT) nodes
	bool *done = calloc(Ntot, sizeof(bool));
	double ftemp;
	
	// first search through each boundary segment
	for(n=0;n<numBound;n++)
	{
		ey = Boundary[n].e / elemX;
		ex = Boundary[n].e - elemX*ey; // element indices
		xMin = ex-irad; xMax=ex+irad+1;
		yMin = ey-irad; yMax=ey+irad+1; // search limits
		xMin = (xMin < 0) ? 0 : xMin;
		yMin = (yMin < 0) ? 0 : yMin;
		xMax = (xMax > elemX) ? elemX : xMax;
		yMax = (yMax > elemY) ? elemY : yMax; // adjust search limits
		
		node=Boundary[n].n1; // 1st node number
		if(!done[node]) // if sens not already computed
		{
			ftemp = 1.0;
			if(node < NumNodes) // If node is a grid node
			{
				num = Lsens(&NodeCoord[node], xMax, xMin, yMax, yMin, alMin, alpha, rad, gCoord, gSens, Number, 4, numDual, sens_temp);
			}		
			else // Otherwise node is an auxillary node
			{
				m = node - NumNodes; // position in AuxNodes array
				num = Lsens(&AuxNodes[m], xMax, xMin, yMax, yMin, alMin, alpha, rad, gCoord, gSens, Number, 4, numDual, sens_temp);
			}
			for(i=0;i<numDual;i++)
			{
				m = Ntot*i + node; // point to correct place in Nsens
				Nsens[m] = ftemp * sens_temp[i]; // multiply smoothed sensitivities by weight
			}
			done[node] = true;
		}
		node=Boundary[n].n2; // 2nd node number
		if(!done[node]) // if sens not already computed
		{
			ftemp = 1.0;
			if(node < NumNodes) // If node is a grid node
			{
				num = Lsens(&NodeCoord[node], xMax, xMin, yMax, yMin, alMin, alpha, rad, gCoord, gSens, Number, 4, numDual, sens_temp);
			}		
			else // Otherwise node is an auxillary node
			{
				m = node - NumNodes; // position in AuxNodes array
				num = Lsens(&AuxNodes[m], xMax, xMin, yMax, yMin, alMin, alpha, rad, gCoord, gSens, Number, 4, numDual, sens_temp);
			}
			for(i=0;i<numDual;i++)
			{
				m = Ntot*i + node; // point to correct place in Nsens
				Nsens[m] = ftemp * sens_temp[i]; // multiply smoothed sensitivities by weight
			}
			done[node] = true;
		}
	}
	
	// clear memory
	free(done);
	free(gSens);
	free(sens_temp);
}
