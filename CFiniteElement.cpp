/*
   CFiniteElement.cpp

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

#include <stdlib.h>
#include <math.h>

#include "CMaterial.h"
#include "CMathUtility.h"
#include "CFiniteElement.h"

CFiniteElement::CFiniteElement() {
	// TODO Auto-generated constructor stub

}

CFiniteElement::~CFiniteElement() {
	// TODO Auto-generated destructor stub
}

//
// Implementation - Interfaces
//

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
						if(inMesh->mat_lin)
						{
							e_fac = mtemp*e_rat + (1.0-mtemp);
						}
                        else
                        {
                        	// jeehanglee@gmail.com: temp code. Should be refactored....
                        	CMaterial material;
                        	e_fac = material.HS_mat(mtemp, 0.5, &inMat[inMesh->mat1], &inMat[inMesh->mat2]) / inMat[inMesh->mat1].e;
                        }
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


// function to compute a cut element area
double InArea(int eStat, int eNum, int *Lnodes, short *NodeStat, int NumNodes, Coord *NodeCoord,
			  Coord *AuxNodes, int NumBound, Bseg *Boundary)
{
	int numPts = 0;
	int i, j, k, temp, temp2;
	Coord pts[6]; // Max number of points will be 6

	// need to change the algorithm if element status = 5
	short nLook = (eStat==5) ? 0:1; // compute area outside structure if ElemStat is 5

	// jeehanglee@gmail.com: temp code. Refactoring required.
	CMathUtility mathUtil;

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

					if(mathUtil.LineCross(chkPts) == 1)
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
	return(mathUtil.PolyArea(numPts, pts));
}

// Array storing data for location of the 4 nodes in the square element
static Coord Q4[4] = {{-1.,-1.},{1.,-1.},{1.,1.},{-1.,1.}};


// The following 2 functions calculate the values in each 2 x 2 block in the Element Stiffness Matrix
// These functions integrate the matrix exactly. Trust me it works

// For odd (i + j = odd) Matrix enrites
double Aodd(int i,int j,double e,double g)
{
	double a = (Q4[i].y * Q4[j].x);
	double b = (Q4[i].x * Q4[j].y);

	double value = (0.25 * ((b * e) + (a * g)));

	return (value);
}

// For even (i + j = even) Matrix enrites
double Aeven(int i,int j,double e,double g)
{
	double a = 1.5 + (0.5 * (Q4[i].y * Q4[j].y));
	a *= (Q4[i].x * Q4[j].x);

	double b = 1.5 + (0.5 * (Q4[i].x * Q4[j].x));
	b *= (Q4[i].y * Q4[j].y);

	double value = (0.166666666667 * ((e * a) + (g * b)));

	return (value);
}

// Complete square IN element stiffness matix computation
// only for isotropic material
void KEMatrix(double *KE, isoMat *mat, double t)
{
	// multiply by constant thickness
	double c = mat->mat[0] * t;
	double d = mat->mat[1] * t;
	double g = mat->mat[8] * t; // symmetric isotropic material!

	int n,m; // Incrementors

	// Populate entire Matrix with the following
	int row,col,temp;

	// For each 2 x 2 block
	for(n=0;n<4;n++)
	{
		for(m=0;m<4;m++)
		{
			row = 2 * n;
			col = 2 * m;

			temp = (8 * row) + col;
			KE[temp] = Aeven(n,m,c,g);
			KE[temp+1] = Aodd(n,m,d,g);
			KE[temp+8] = Aodd(n,m,g,d);
			KE[temp+9] = Aeven(n,m,g,c);
		}
	}
}

// consistent mass matrix for a 2D plane element (just needs scaling)
static const double me_const[64]= {
	4,0,2,0,1,0,2,0,
	0,4,0,2,0,1,0,2,
	2,0,4,0,2,0,1,0,
	0,2,0,4,0,2,0,1,
	1,0,2,0,4,0,2,0,
	0,1,0,2,0,4,0,2,
	2,0,1,0,2,0,4,0,
	0,2,0,1,0,2,0,4};

// consistent mass matrix for a 2D plane element
void MEMatrix(double *ME, double rho, double area, double t)
{
	double mass = t*area*rho/36.0; // mass = volume x density
	int i;
	for(i=0;i<64;i++)
	{
		ME[i] = me_const[i]*mass;
	}
}

// Function that assemebles the element matricies into the global matrix (in triplet form)
// Symmetric matrix - upper triangle
void Assemble2(int *K_begin, int *M_begin, int nNod, int nDof, int *tnodes, double *KE, double *ME,
					sp_mat *Kg, sp_mat *Mg, int mass)
{
	int start,row,col,n,m,rn,cm,k_num,m_num,p,q,i,j,temp; // Incrementors
	k_num = *K_begin;
	if(mass==1){m_num = *M_begin;}
	int dof = nDof * nNod; // total dof for element
	double ftemp;

	for(i=0;i<nNod;i++)
	{
		for(j=i;j<nNod;j++)
		{
			p = (tnodes[i] <= tnodes[j]) ? i : j;
			q = (tnodes[j] <= tnodes[i]) ? i : j; // Ensure upper triangle entries are read in
			rn = nDof * p;
			cm = nDof * q; // element matrix row and column numbers -> for upper triangle
			row = tnodes[p] * nDof;
			col = tnodes[q] * nDof; // Global matrix row and column numbers -> for upper triangle
			start = 0;
			for(n=0;n<nDof;n++){
				for(m=start;m<nDof;m++){
					temp = (dof * (n + rn)) + m + cm; // index correct place in elem matrix
					Kg->irn[k_num] = row + n;
					Kg->jcn[k_num] = col + m;
					Kg->A[k_num] = KE[temp];
					if(mass==1)// also make mass matrix if mass=1
					{
						ftemp = ME[temp];
						if(fabs(ftemp) > 1.0e-12) // ignore zero entries
						{
							Mg->irn[m_num] = Kg->irn[k_num];
							Mg->jcn[m_num] = Kg->jcn[k_num];
							Mg->A[m_num++] = ftemp;
						}
					}
					k_num++;
				}
				start = (i == j) ? n+1 : 0;	// Ensures only upper triangle for diagonal blocks are read in
			}
		}
	}

	*K_begin = k_num; // read back number of entries added
	if(mass==1){*M_begin = m_num;}
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

