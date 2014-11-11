/*
 *  EMatrix.c
 *
 *  Created by Peter Dunning
 *	Version 5.4 @ 23/04/2014
 *  Functions to compute element matrices
 */

#include "EMatrix.h"
#include <stdlib.h>
#include <math.h>

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

// function to compute elastic modulus from two materials
double HS_mat(double alpha, double hs_int, isoMat *mat1, isoMat *mat2)
{
    double kmax, kmin, gmax, gmin, khs, ghs;
    double fact, ftemp, ftemp2;
    double amin = 1.0-alpha;
    
    // bulk modulus
    fact = amin*mat1->k + alpha*mat2->k;
    ftemp = mat2->k - mat1->k;
    ftemp *= ftemp;
    ftemp *= amin*alpha;
    
    kmax = fact - (ftemp / (amin*mat2->k + alpha*mat1->k + mat2->g));
    kmin = fact - (ftemp / (amin*mat2->k + alpha*mat1->k + mat1->g));
    
    khs = hs_int*kmax + (1.0-hs_int)*kmin;
    
    // shear modulus
    fact = amin*mat1->g + alpha*mat2->g;
    ftemp = mat2->g - mat1->g;
    ftemp *= ftemp;
    ftemp *= amin*alpha;
    
    ftemp2 = (mat2->k * mat2->g) / (mat2->k + 2.0*mat2->g);
    gmax = fact - (ftemp / (amin*mat2->g + alpha*mat1->g + ftemp2));
    
    ftemp2 = (mat1->k * mat1->g) / (mat1->k + 2.0*mat1->g);
    gmin = fact - (ftemp / (amin*mat2->g + alpha*mat1->g + ftemp2));
    
    ghs = hs_int*gmax + (1.0-hs_int)*gmin;
    
    // Elastic modulus
    ftemp = (9.0*khs*ghs)/(3.0*khs + ghs);
    return (ftemp);
}

// function to compute self-weight load vector
void self_weight(mesh *inMesh, isoMat *inMat, double aMin, double mMin, double *alpha,
                    int freeDof, int *dofMap, int numCase, double *load_in, double *load_out, Coord *acc)
{
    int i,j,n,m,o,num,ind,temp,temp2;
    double atemp, fx,fy;
    int tnodes[4], Xdof[4], Ydof[4];
    double mtemp,p_fac,p_rat; // variables for designable material
    int mat_count = 0; // count for designable material
    // material ratio
	if(inMesh->des_mat)
	{
		p_rat = inMat[inMesh->mat2].rho / inMat[inMesh->mat1].rho;
	}
    
    // read data
    int totDof = inMesh->NumNodes * NUM_DOF;
    int elemX = inMesh->elemX;
	int elemY = inMesh->elemY;
    Elem **Number = inMesh->Number;
    double vol = inMesh->t * inMesh->h * inMesh->h; // element volume
    
    // copy load_in to load_temp (expanding)
    double *load_temp = calloc(totDof*numCase, sizeof(double));
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
            o = inMesh->mat_type[num]; // material number
            
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
            if(inMesh->des_mat)
            {
                // material mix variable
                mtemp = inMesh->mat_vars[mat_count++];
                
                // density factor
                p_fac =  mtemp*p_rat + (1.0-mtemp);
            }
            else { p_fac = 1.0; }
            
            p_fac *= (atemp == aMin) ? mMin : atemp; // total density factor
            p_fac *= 0.25 * inMat[o].rho * vol; // elem mass / 4 (nodes)
            
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
