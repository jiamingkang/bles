/*
   CMaterial.cpp

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

#include "CMaterial.h"

CMaterial::CMaterial() {
	// TODO Auto-generated constructor stub

}

CMaterial::~CMaterial() {
	// TODO Auto-generated destructor stub
}

//
// Implementation - Interfaces
//

// function to compute elastic modulus from two materials
double CMaterial::HS_mat(double alpha, double hs_int, isoMat *mat1, isoMat *mat2)
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
void CMaterial::self_weight(mesh *inMesh, isoMat *inMat, double aMin, double mMin, double *alpha,
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
