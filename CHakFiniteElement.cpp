/*
   CHakFiniteElement.cpp

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

#include "CHakMaterial.h"
#include "CHakMathUtility.h"
#include "CHakFiniteElement.h"

CHakFiniteElement::CHakFiniteElement() {
	// TODO Auto-generated constructor stub

}

CHakFiniteElement::~CHakFiniteElement() {
	// TODO Auto-generated destructor stub
}

//
// Implementation - Interfaces
//
// knai20@bath.ac.uk: OOD version

// Assembles global stiffness (& maybe mass) matrix for AFG method in triplet format (for MA57 solver)
void CHakFiniteElement::AFG_Matrix(int mass, double **KE, double **ME, CHakSparseMatrix *Kg, CHakSparseMatrix *Mg, CHakSparseMatrix *lump_mass, double *alpha,
					CHakMesh *inMesh, CHakMaterial *inMat, double aMin, double mMin)
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
	if(inMesh->m_bDesignableMaterial)
	{
		e_rat = inMat[inMesh->m_materialTwo].m_e / inMat[inMesh->m_materialOne].m_e;
		p_rat = inMat[inMesh->m_materialTwo].m_rho / inMat[inMesh->m_materialOne].m_rho;
	}

	// read data
	int elemX = inMesh->m_elemX;
	int elemY = inMesh->m_elemY;
	Elem **Number = inMesh->m_pNumber;
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

				o = inMesh->m_pMaterialType[num]; // material number
				kptr = KE[o]; mptr = ME[o]; // point to correct matrices for elem material

				// Read in grid node numbers of current element
				tnodes[0] = Number[n][m].a;
				tnodes[1] = Number[n][m].b;
				tnodes[2] = Number[n][m].c;
				tnodes[3] = Number[n][m].d;

				// If element is IN, no need to multiply by alpha (unless using designable material)
				if(!inMesh->m_bDesignableMaterial && atemp > 0.9999)
				{
					// Assemble area ratio weighted stiffness matrix into index arrays
					Assemble2(&K_begin, &M_begin, 4, NUM_DOF, tnodes, kptr, mptr, Kg, Mg, mass);
				}

				else
				{
					// modify matrices for designable materials
					if(inMesh->m_bDesignableMaterial && num==inMesh->m_pMaterialElems[mat_count])
					{
						// material mix variable
						mtemp = inMesh->m_pMaterialVars[mat_count++];

						// modulus factor
						if(inMesh->m_bMaterialLin)
						{
							e_fac = mtemp*e_rat + (1.0-mtemp);
						}
                        else
                        {
                        	// temp code. Should be refactored....
                        	CHakMaterial material;
                        	e_fac = material.HS_mat(mtemp, 0.5, &inMat[inMesh->m_materialOne], &inMat[inMesh->m_materialTwo]) / inMat[inMesh->m_materialOne].m_e;
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
            else if(inMesh->m_bDesignableMaterial && num==inMesh->m_pMaterialElems[mat_count]){mat_count++;}
		}
	}

	// Add lumped masses
	m = lump_mass->m_numEntry;
	if(m>0)
	{
		for(n=0;n<m;n++)
		{
			Mg->m_indRow[M_begin] = lump_mass->m_indRow[n];
			Mg->m_indCol[M_begin] = lump_mass->m_indCol[n];
			Mg->m_pMatEntry[M_begin++] = lump_mass->m_pMatEntry[n];
		}
	}

	// Add bar elements
	if(inMesh->m_bars)
	{
		int numBar = inMesh->m_numBars;
		double EoL = inMesh->m_pBarMaterialType->m_e / inMesh->m_lenEdge; // E / L
		double kbar;

		for(n=0;n<numBar;n++)
		{
			m = inMesh->m_pBarNums[n].n1;
			o = inMesh->m_pBarNums[n].n2;
			kbar = EoL * inMesh->m_pBarAreas[n]; // bar stiffness

			// swap if necessary (m < o)
			if(m > o){num=o; o=m; m=num;}

			// 3 entries
			Kg->m_indRow[K_begin] = m; Kg->m_indCol[K_begin] = m; Kg->m_pMatEntry[K_begin++] = kbar;
			Kg->m_indRow[K_begin] = o; Kg->m_indCol[K_begin] = o; Kg->m_pMatEntry[K_begin++] = kbar;
			Kg->m_indRow[K_begin] = m; Kg->m_indCol[K_begin] = o; Kg->m_pMatEntry[K_begin++] = -kbar;
		}
	}

	// Add designable bc stiffness
	if(inMesh->m_bDesignableBc)
	{
		int numBC = inMesh->m_numBc;

		for(o=0;o<numBC;o++)
		{
			// compute element indices
			m = floor((double)inMesh->m_pBcNums[o] / (double)elemX);
			n = inMesh->m_pBcNums[o] - m*elemX;

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
					Kg->m_indRow[K_begin] = num+m;
					Kg->m_indCol[K_begin] = Kg->m_indRow[K_begin];
					Kg->m_pMatEntry[K_begin++] = pow(inMesh->m_K_bc[o],3.0)*inMesh->m_K0_bc; // penalize design variable
				}
			}
		}
	}

	// resize matrices
	Kg->m_numEntry = K_begin;
	Mg->m_numEntry = M_begin;
}


// function to compute a cut element area
double CHakFiniteElement::InArea(int eStat, int eNum, int *Lnodes, short *NodeStat, int NumNodes, Coord *NodeCoord,
			  Coord *AuxNodes, int NumBound, Bseg *Boundary)
{
	int numPts = 0;
	int i, j, k, temp, temp2;
	Coord pts[6]; // Max number of points will be 6

	// need to change the algorithm if element status = 5
	short nLook = (eStat==5) ? 0:1; // compute area outside structure if ElemStat is 5

	// jeehanglee@gmail.com: temp code. Refactoring required.
	CHakMathUtility mathUtil;

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
double CHakFiniteElement::Aodd(int i,int j,double e,double g)
{
	double a = (Q4[i].y * Q4[j].x);
	double b = (Q4[i].x * Q4[j].y);

	double value = (0.25 * ((b * e) + (a * g)));

	return (value);
}

// For even (i + j = even) Matrix enrites
double CHakFiniteElement::Aeven(int i,int j,double e,double g)
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
void CHakFiniteElement::KEMatrix(double *KE, CHakMaterial *mat, double t)
{
	// multiply by constant thickness
	double c = mat->m_mat[0] * t;
	double d = mat->m_mat[1] * t;
	double g = mat->m_mat[8] * t; // symmetric isotropic material!

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
void CHakFiniteElement::MEMatrix(double *ME, double rho, double area, double t)
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
void CHakFiniteElement::Assemble2(int *K_begin, int *M_begin, int nNod, int nDof, int *tnodes, double *KE, double *ME,
					CHakSparseMatrix *Kg, CHakSparseMatrix *Mg, int mass)
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
					Kg->m_indRow[k_num] = row + n;
					Kg->m_indCol[k_num] = col + m;
					Kg->m_pMatEntry[k_num] = KE[temp];
					if(mass==1)// also make mass matrix if mass=1
					{
						ftemp = ME[temp];
						if(fabs(ftemp) > 1.0e-12) // ignore zero entries
						{
							Mg->m_indRow[m_num] = Kg->m_indRow[k_num];
							Mg->m_indCol[m_num] = Kg->m_indCol[k_num];
							Mg->m_pMatEntry[m_num++] = ftemp;
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
void CHakFiniteElement::free_sp_mat(CHakSparseMatrix *m)
{
	if(m->m_indRow){free(m->m_indRow); m->m_indRow=0;}
	if(m->m_indCol){free(m->m_indCol); m->m_indCol=0;}
	if(m->m_pMatEntry){free(m->m_pMatEntry); m->m_pMatEntry=0;}
}

// function to create memory for a sparse matrix
void CHakFiniteElement::set_sp_mat(CHakSparseMatrix *m)
{
	// first ensure memory if free
	free_sp_mat(m);

	// then set memory
	int len = m->m_numEntry;
	if(len > 0)
	{
		m->m_indRow = (int *) malloc(len*sizeof(int));
		m->m_indCol = (int *) malloc(len*sizeof(int));
		m->m_pMatEntry = (double *) malloc(len*sizeof(double));
	}
}

// function to remove dof from a sparse matrix
void CHakFiniteElement::rem_sp_mat(CHakSparseMatrix *m, int *map, int inc)
{
	int i,j,k;

	// read data
	int numEnt = m->m_numEntry;
	int *irn = m->m_indRow;
	int *jcn = m->m_indCol;
	double *A = m->m_pMatEntry;

	int temp_ent = 0; // count number entries
	for(i=0;i<numEnt;i++)
	{
		if(irn[i] > -1)
		{
			j = map[irn[i]];
			k = map[jcn[i]];

			if( (j > -1) && (k > -1) ) // if dof not fixed
			{
				m->m_pMatEntry[temp_ent] = A[i];
				m->m_indRow[temp_ent] = j + inc;
				m->m_indCol[temp_ent++] = k + inc; // use inc = 1 for Fortran solver
			}
		}
	}

	// reduce entries
	m->m_numEntry = temp_ent;
}


// function to compute self-weight load vector
void CHakFiniteElement::SelfWeight(CHakMesh *inMesh, CHakMaterial *inMat, double aMin, double mMin, double *alpha,
                    int freeDof, int *dofMap, int numCase, double *load_in, double *load_out, Coord *acc)
{
    int i,j,n,m,o,num,ind,temp,temp2;
    double atemp, fx,fy;
    int tnodes[4], Xdof[4], Ydof[4];
    double mtemp,p_fac,p_rat; // variables for designable material
    int mat_count = 0; // count for designable material
    // material ratio
	if(inMesh->m_bDesignableMaterial)
	{
		p_rat = inMat[inMesh->m_materialTwo].m_rho / inMat[inMesh->m_materialOne].m_rho;
	}

    // read data
    int totDof = inMesh->m_numNodes * NUM_DOF;
    int elemX = inMesh->m_elemX;
	int elemY = inMesh->m_elemY;
    Elem **Number = inMesh->m_pNumber;
    double vol = inMesh->m_thinkness * inMesh->m_lenEdge * inMesh->m_lenEdge; // element volume

    // copy load_in to load_temp (expanding)
    double *load_temp = (double *) calloc(totDof*numCase, sizeof(double));
    for(i=0;i<totDof;i++)
    {
        if(dofMap[i] != -1)	// if dof not fixed copy accross the value
        {
            for(j=0;j<numCase;j++)
            {
                temp = j * totDof; // to get to correct place in disp
                temp2 = j * freeDof; // to get to correct place in disp_temp
                load_temp[i+temp] = load_in[dofMap[i]+temp2];
            }
        }
        // else - already zero
    }

    // For All Elements
	for(m=0;m<elemY;m++)
	{
		for(n=0;n<elemX;n++)
		{
            num = Number[n][m].n; // Element number
            atemp = alpha[num]; // Element area ratio
            o = inMesh->m_pMaterialType[num]; // material number

            // Read in grid node numbers of current element
            tnodes[0] = Number[n][m].a;
            tnodes[1] = Number[n][m].b;
            tnodes[2] = Number[n][m].c;
            tnodes[3] = Number[n][m].d;

            // for each node, compute global dof
            for(i=0;i<4;i++)
            {
                j = NUM_DOF*tnodes[i];
                Xdof[i] = j; Ydof[i] = j+1;
            }

            // compute mass of the element
            if(inMesh->m_bDesignableMaterial)
            {
                // material mix variable
                mtemp = inMesh->m_pMaterialVars[mat_count++];

                // density factor
                p_fac =  mtemp*p_rat + (1.0-mtemp);
            }
            else { p_fac = 1.0; }

            p_fac *= (atemp == aMin) ? mMin : atemp; // total density factor
            p_fac *= 0.25 * inMat[o].m_rho * vol; // elem mass / 4 (nodes)

            if(atemp > 0.0)
            {
                // compute element self-weight load vector for each load case
                for(j=0;j<numCase;j++)
                {
                    ind = totDof*j; // start of current load case

                    // nodal forces (mass x acceleration)
                    fx = p_fac * acc[j].x;
                    fy = p_fac * acc[j].y;

                    // for each node add self-weight forces
                    for(i=0;i<4;i++)
                    {
                        load_temp[ind+Xdof[i]] += fx;
                        load_temp[ind+Ydof[i]] += fy;
                    }
                }
            }
        }
    }

    // Additional load for bar elements
    {
        // ???
    }

    // finally remove fixed dof from load_temp -> load_out
    for(i=0;i<totDof;i++) // for all dof
    {
        ind = dofMap[i];
        if(ind > -1) // if dof not fixed
        {
            for(j=0;j<numCase;j++) // for all load cases
            {
                n = totDof * j; // point to loaction in load
                m = freeDof * j; // point to loaction in load_in
                load_out[m+ind] = load_temp[n+i]; // copy load accross
            }
        }
    }

    free(load_temp);
}

