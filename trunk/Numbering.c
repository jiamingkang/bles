/*
 *  Numbering.c
 *
 *  Created by Peter Dunning
 *	Version 5.4 @ 23/04/2014
 *  Funcitons to number the fixed grid elements & nodes, compute node coords
 */

#include "Numbering.h"
#include <stdlib.h>
#include <math.h>

// Function that numbers all elements and nodes in the FG domain
void Numbering(mesh *inMesh)
{
	int i,j;
	int num = 0; // first element number is zero
	
	// read data
	Elem **Number = inMesh->Number; // pointer to Numbering array
	int elemX = inMesh->elemX;
	int elemY = inMesh->elemY;
	
	// Element Numbering loop 
	for(j=0;j<elemY;j++)
	{
		for(i=0;i<elemX;i++)
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
	if(elemX > 1) // if there is more than one element in the x direction number first row
	{
		for(i=1;i<elemX;i++)
		{
			Number[i][0].a = Number[i-1][0].b;
			Number[i][0].b = ++num;
			Number[i][0].c = ++num;
			Number[i][0].d = Number[i-1][0].c;
		}
	}
	
	// Number rest of FG mesh
	if(elemY > 1)
	{
		for(j=1;j<elemY;j++)
		{
			for(i=0;i<elemX;i++)
			{
				Number[i][j].a = Number[i][j-1].d;
				Number[i][j].b = Number[i][j-1].c;
				
				if(i==0)
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
void Coordinates(mesh *inMesh)
{
	// data from inMesh
	int elemX = inMesh->elemX;
	int elemY = inMesh->elemY;
	double hx = inMesh->h; double hy = inMesh->h;
	double tol = inMesh->tol;
	Elem **Number = inMesh->Number; // pointer to Numbering array
	Coord *NodeCoord = inMesh->NodeCoord; // pointer to node coords
	
	int i,j,num;
	
	for(j=0;j<elemY;j++)
	{
		for(i=0;i<elemX;i++)
		{
			num = Number[i][j].a;
			// If node co-ordinate hasn't been calculated already
			if((NodeCoord[num].x + NodeCoord[num].y) < tol){   
				NodeCoord[num].x = i * hx;
				NodeCoord[num].y = j * hy;
			}
			num = Number[i][j].b;
			if((NodeCoord[num].x + NodeCoord[num].y) < tol){
				NodeCoord[num].x = (i * hx) + hx;
				NodeCoord[num].y = j * hy;
			}
			num = Number[i][j].c;
			if((NodeCoord[num].x + NodeCoord[num].y) < tol){
				NodeCoord[num].x = (i * hx) + hx;
				NodeCoord[num].y = (j * hy) + hy;
			}
			num = Number[i][j].d;
			if((NodeCoord[num].x + NodeCoord[num].y) < tol){
				NodeCoord[num].x = i * hx;
				NodeCoord[num].y = (j * hy) + hy;
			}
		}
	}
}

// function that orders node numbers into a 2D based on their relative positions
void NodeNums2(mesh *inMesh)
{
	// data from inMesh
	int NodeX = inMesh->NodeX;
	int NodeY = inMesh->NodeY;
	int NumNodes = inMesh->NumNodes;
	double h = inMesh->h;
	double tol = inMesh->tol;
	int **Nodes2 = inMesh->Nodes2; // pointer to Node numbering array
	Coord *NodeCoord = inMesh->NodeCoord; // pointer to node coords

	
	int i; // incrementor
	int X,Y; // node position variables
	
	// for all grid node determine position in ordered 2D array
	//  by their co-ordinates and element edge length
	for(i=0;i<NumNodes;i++)
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
	for(i=1;i<NodeX-1;i++)
	{
		Nodes2[i][0] = Nodes2[i][1];
		Nodes2[i][NodeY-1] = Nodes2[i][NodeY-2];
	}
	
	// Finally fill in left and right columns
	for(i=1;i<NodeY-1;i++)
	{
		Nodes2[0][i] = Nodes2[1][i];
		Nodes2[NodeX-1][i] = Nodes2[NodeX-2][i];
	}
}

// function to number bars elements
void Bar_numbering(mesh *inMesh)
{
	// read data
	Elem **Number = inMesh->Number; // pointer to Numbering array
	int elemX = inMesh->elemX;
	int elemY = inMesh->elemY;
	
	int i,j;
	int xend = elemX-1;
	int yend = elemY-1;
	int count = 0;
	
	// loop through all elements
	for(j=0;j<elemY;j++)
	{
		for(i=0;i<elemX;i++)
		{
			// use bottom edge as x-bars
			inMesh->bar_nums[count].e = count; // element number (bar elems only)
			inMesh->bar_nums[count].n1 = Number[i][j].a * 2; // x-dof of node 1
			inMesh->bar_nums[count++].n2 = Number[i][j].b * 2; // x-dof of node 2
			
			// also use top edge for top row of elements
			if(j == yend)
			{
				inMesh->bar_nums[count].e = count; // element number (bar elems only)
				inMesh->bar_nums[count].n1 = Number[i][j].d * 2; // x-dof of node 1
				inMesh->bar_nums[count++].n2 = Number[i][j].c * 2; // x-dof of node 2
			}
		}
	}
	
	// loop through all elements
	for(i=0;i<elemX;i++)
	{
		for(j=0;j<elemY;j++)
		{
			// use left edge as y-bars
			inMesh->bar_nums[count].e = count; // element number (bar elems only)
			inMesh->bar_nums[count].n1 = (Number[i][j].a * 2) +1; // y-dof of node 1
			inMesh->bar_nums[count++].n2 = (Number[i][j].d * 2) +1; // y-dof of node 2
			
			// also use right edge for last column of elements
			if(i == xend)
			{
				inMesh->bar_nums[count].e = count; // element number (bar elems only)
				inMesh->bar_nums[count].n1 = (Number[i][j].b * 2) +1; // y-dof of node 1
				inMesh->bar_nums[count++].n2 = (Number[i][j].c * 2) +1; // y-dof of node 2
			}
		}
	}
	
	// that is all!!
}

// function to free memory for a sparse matrix
void free_sp_mat(sp_mat *m)
{
	if(m->irn){free(m->irn); m->irn=0;}
	if(m->jcn){free(m->jcn); m->jcn=0;}
	if(m->A){free(m->A); m->A=0;}
}

// function to create memory for a sparse matrix
void set_sp_mat(sp_mat *m)
{
	// first ensure memory if free
	free_sp_mat(m);
	
	// then set memory
	int len = m->ne;
	if(len > 0)
	{
		m->irn = malloc(len*sizeof(int));
		m->jcn = malloc(len*sizeof(int));
		m->A = malloc(len*sizeof(double));
	}
}

// function to remove dof from a sparse matrix
void rem_sp_mat(sp_mat *m, int *map, int inc)
{
	int i,j,k;
	
	// read data
	int numEnt = m->ne;
	int *irn = m->irn;
	int *jcn = m->jcn;
	double *A = m->A;
	
	int temp_ent = 0; // count number entries
	for(i=0;i<numEnt;i++)
	{
		if(irn[i] > -1)
		{
			j = map[irn[i]];
			k = map[jcn[i]];
			
			if( (j > -1) && (k > -1) ) // if dof not fixed
			{
				m->A[temp_ent] = A[i];
				m->irn[temp_ent] = j + inc;
				m->jcn[temp_ent++] = k + inc; // use inc = 1 for Fortran solver
			}
		}
	}
	
	// reduce entries
	m->ne = temp_ent;
}
