/*
 * CHak2DTri3.cpp
 *
 *  Created on: 23 Feb 2015
 *      Author: knai20-admin
 */

#include "CHak2DTri3.h"

CHak2DTri_3::CHak2DTri_3() {
	// TODO Auto-generated constructor stub
	// constructor sets num nodes
    numNode=3;
	nodeSpace();
	Ke=0;
	Me=0;
    t = -1.0;
    volume = -1.0;

}

CHak2DTri_3::~CHak2DTri_3() {
	// TODO Auto-generated destructor stub
	// destructor (delete memory)
	delMat();
}


// compute the Tri3 elem stiffness matrix
void CHak2DTri_3::makeKe()
{
    int g,i,j,ij; // incrementors

    if(Ke){delete [] Ke;}
	Ke = new double [45]();

    // read in node coordinates
	Coord3D np[3]; getNcrd(np);

    // compute strain displacement matrix
    double factors[7]; // geometric factors, factors[6] = 2xarea

    // compute common geometric factors for the triangle
    tri3_geomFact(np, factors);

	double B[27];
    double BtD[27];
    double *D = Mat->getDsts(); // plane stress property matrix

	Coord gauss[3];
	gauss[0].x = 0.5;
	gauss[0].y = 0.0;
	gauss[1].x = 0.5;
	gauss[1].y = 0.5;
	gauss[2].x = 0.0;
	gauss[2].y = 0.5; // coord for gauss intergration rule

	double sub;

	// constant factor over entire matrix (for const thk)
	// includes weighting from gauss intergration
	double t_fact = t / (6.0*factors[6]);

    for(g=0;g<3;g++)
	{
        // compute the B matrix at the input point
        tri3_Bmat(gauss[g], factors, B);

		// compute stiffness matrix

        // first multiply B^T * D
		cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, 9, 3, 3, 1.0, B, 9, D, 3, 0.0, BtD, 3);

		// complete multiplication (for upper triangle only)
		for(i=0;i<9;i++) // rows
		{
			for(j=i;j<9;j++) // columns (upper triangle)
			{
				sub = 0.0;
				for(ij=0;ij<3;ij++)
				{
					sub += BtD[index2(i,ij,3)]*B[index2(ij,j,9)];
				}
				Ke[tri_ind(9,i,j)] += sub*t_fact;
			}
		}
	}

}

// function to compute element mass matrix (from OPT shape functions)
void CHak2DTri_3::makeMe()
{
    int g,i,j,ij; // incrementors

    if(Me){delete [] Me;}
	Me = new double [45]();

    if(volume < 0.0){ makeVol(); } // compute volume if not done already

    // total mass
    double mass = volume * Mat->getDens(); // mass = vol x density

	double e,n;		// non-dim coords
	double wMod,sub; // intermediate values
	double Nw[18];   // need 2 rows

    // read in node coordinates
	Coord3D np[3]; getNcrd(np);

    // compute common geometric factors for the triangle
    double factors[7]; // geometric factors, factors[6] = 2xarea
    tri3_geomFact(np, factors);

    // integration points and weights
	static double tri3_e_gauss[7] = {1.0/3.0,0.47014,0.05971,0.47014,0.10128,0.79742,0.10128};
	static double tri3_n_gauss[7] = {1.0/3.0,0.47014,0.47014,0.05971,0.10128,0.10128,0.79742};
	static double tri3_weight[7] = {0.225,0.13239,0.13239,0.13239,0.12593,0.12593,0.12593};

    for(g=0;g<7;g++)
	{
        e = tri3_e_gauss[g];
		n = tri3_n_gauss[g];

		tri3_Nw(e,n,factors,Nw);

		wMod = tri3_weight[g]*mass; // combined weight

        // compute upper triangle entries
		for(i=0;i<9;i++) // row
		{
			for(j=i;j<9;j++) // column
			{
                sub = 0.0;
				for(ij=0;ij<2;ij++)
				{
					sub += Nw[index2(ij,i,9)]*Nw[index2(ij,j,9)]; // N^T x N
				}
				Me[tri_ind(9,i,j)] += wMod*sub; // weight * N^T x N
			}
		}
	}
}

// function to create consistent mass matrix (from LST shape functions)
void CHak2DTri_3::makeMe2()
{
    int g,i,j,k,i2,j2,ind; // incrementors

    if(Me){delete [] Me;}
	Me = new double [45]();

    if(volume < 0.0){ makeVol(); } // compute volume if not done already

    // total mass
    double mass = volume * Mat->getDens(); // mass = vol x density

	double e,n,s1;		// non-dim coords
	double wMod,sub; // intermediate values
	double Nw[6];   // for LST shape functions

    // read in node coordinates
	Coord3D np[3]; getNcrd(np);

    // compute common geometric factors for the triangle
    double factors[7]; // geometric factors, factors[6] = 2xarea
    tri3_geomFact(np, factors);

    // integration points and weights
	static double tri3_e_gauss[7] = {1.0/3.0,0.47014,0.05971,0.47014,0.10128,0.79742,0.10128};
	static double tri3_n_gauss[7] = {1.0/3.0,0.47014,0.47014,0.05971,0.10128,0.10128,0.79742};
	static double tri3_weight[7] = {0.225,0.13239,0.13239,0.13239,0.12593,0.12593,0.12593};
	double *Mlst = new double[144](); // full LST mass matrix

    for(g=0;g<7;g++)
	{
        e = tri3_e_gauss[g];
		n = tri3_n_gauss[g];
		s1 = 1.0-e-n;

		Nw[0] = s1*(2.0*s1-1.0); Nw[1] = e*(2.0*e-1.0); Nw[2] = e*(2.0*e-1.0);
        Nw[3] = 4.0*s1*e; Nw[4] = 4.0*e*n; Nw[5] = 4.0*s1*n;

		wMod = tri3_weight[g]*mass; // combined weight

        // compute entries
		for(i=0;i<6;i++) // node row
		{
			for(j=0;j<6;j++) // node column
			{
                sub = Nw[i]*Nw[j];
                i2 = 2*i; j2=2*j; // dof for u
				Mlst[index2(i2,j2,12)] += wMod*sub; // weight * N^T x N
                i2++; j2++; // dof for v
				Mlst[index2(i2,j2,12)] += wMod*sub; // weight * N^T x N
			}
		}
	}

    // compute dof tranformation matrix
    double *LST_trans = new double[108](); // 12 x 9 matrix
    LST_trans[0] = 1.0; LST_trans[10] = 1.0; LST_trans[21] = 1.0; LST_trans[31] = 1.0; LST_trans[42] = 1.0; LST_trans[52] = 1.0;
    LST_trans[54] = 0.5; LST_trans[56] = 0.125*factors[2]; LST_trans[57] = 0.5; LST_trans[59] = -0.125*factors[2];
    LST_trans[64] = 0.5; LST_trans[65] = 0.125*factors[5]; LST_trans[67] = 0.5; LST_trans[68] = -0.125*factors[5];
    LST_trans[75] = 0.5; LST_trans[77] = 0.125*factors[0]; LST_trans[78] = 0.5; LST_trans[80] = -0.125*factors[0];
    LST_trans[85] = 0.5; LST_trans[86] = 0.125*factors[3]; LST_trans[88] = 0.5; LST_trans[89] = -0.125*factors[3];
    LST_trans[90] = 0.5; LST_trans[92] = -0.125*factors[1]; LST_trans[96] = 0.5; LST_trans[98] = 0.125*factors[1];
    LST_trans[100] = 0.5; LST_trans[101] = -0.125*factors[4]; LST_trans[106] = 0.5; LST_trans[107] = 0.125*factors[4];

    // apply transformation
    double mtemp[108]; // temp matrix used in multiplication
    double Me_temp[81];
    // whole matrix
    cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, 9, 12, 12, 1.0, LST_trans, 9, Mlst, 12, 0.0, mtemp, 12);
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 9, 9, 12, 1.0, mtemp, 12, LST_trans, 9, 0.0, Me_temp, 9);
    delete [] LST_trans;
    delete [] Mlst;

    // reduce matrix to upper triangle only (consistent mass matrix)
    for(i=0;i<9;i++) // row
    {
        for(j=i;j<9;j++) // col
        {
            ind = tri_ind(9,i,j); // upper triangle indicator
            k = index2(i,j,9); // full matrix indicator
            Me[ind] = Me_temp[k];
        }
    }
}


// function to compute element volume
void CHak2DTri_3::makeVol()
{
    Coord3D p[3]; getNcrd(p);
    triangle tri(p);
    volume = tri.area() * t; // volume = area x thickness
}

void CHak2DTri_3::getDisp(double *elDisp, double *globDisp)
{
	int i,nd;
	int ind = 0;

	for(i=0;i<3;i++) // for each node
	{
        if(dof_ind)
        {
            nd = dof_ind[getGlobal(i)];
            elDisp[ind++] = globDisp[nd];
            elDisp[ind++] = globDisp[nd+1];
            elDisp[ind++] = globDisp[nd+2]; // assume all elems have same dof
        }
        else
        {
            nd = getGlobal(i)*6;  // 1st dof
            elDisp[ind++] = globDisp[nd]; // u
            elDisp[ind++] = globDisp[nd+1]; // v
            elDisp[ind++] = globDisp[nd+5]; // rot_w
        }
	}
}

// compute Tri3 strain tensor at a (non-dim) point in the element
void CHak2DTri_3::strain(Coord p, double *elDisp, double *stn)
{
    // p must be non-dimensional

	double B[27];	// Strain displacement matrix

    // read in node coordinates
	Coord3D np[3]; getNcrd(np);

    // compute common geometric factors for the triangle
    double factors[7]; // geometric factors, factors[6] = 2xarea
    tri3_geomFact(np, factors);

	// compute isoparametric B matrix
    tri3_Bmat(p, factors, B);

	// multiply: strain = B x u
	cblas_dgemv(CblasRowMajor, CblasNoTrans, 3, 9, 1.0, B, 9, elDisp, 1, 0.0, stn, 1);

    // need to dividie by 2 x area (not in original B matrix)
    stn[0] /= factors[6];
    stn[1] /= factors[6];
    stn[2] /= factors[6];
}

// compute Tri3 element thickness sensitivity
void CHak2DTri_3::elemSens(int numDual, int numCase, double **disp_prim, double **disp_dual, double **fg,
                   double *wgt_fact, double *wgt_case, double *elSens, int eNum, double *Nf)
{
    int l,d;
    double const_fact;
    double ud[9];  // Primary element displacement array
    double dd[9];  // Dual element displacement array
    double Ku[9]; double fg_temp;
    double sens; // variable to keep running total of a sensitivity value
    double K_temp[81];
    double ftemp;

    // expand triangle storage
    int i,j;
    for(i=0;i<9;i++) // row
    {
        for(j=0;j<9;j++) // col
        {
            K_temp[index2(i,j,9)] = (j < i) ? Ke[tri_ind(9,j,i)] : Ke[tri_ind(9,i,j)];
        }
    }

    // for each load case
    for(l=0;l<numCase;l++)
    {
        // compute the constant sensitivty factor for the element
        const_fact = (wgt_fact[l] * Mat->getDens() * volume) / t; // wgt_fact x dW/dt

        // get the primary element displacements
        getDisp(ud, disp_prim[l]);

        // compute pKu
        // multiply: K x u
        cblas_dgemv(CblasRowMajor, CblasNoTrans, 9, 9, 1.0, K_temp, 9, ud, 1, 0.0, Ku, 1);
        fg_temp = cblas_ddot(9, ud, 1, fg[l], 1);

        // for each dual response p
        for(d=0;d<numDual;d++)
        {
            // add the 3 componnents (multipled by the load case weight)
            // -pKu + fg(u+p) + const_fact

            // get the dual element displacements
            getDisp(dd, disp_dual[index2(l,d,numDual)]);

            // multiply:p x Ku
            sens = -1.0 * cblas_ddot(9, dd, 1, Ku, 1);

            // multiply fg(u+p)
            ftemp = fg_temp + cblas_ddot(9, dd, 1, fg[l], 1);
            sens += Nf[l] * ftemp / t; // N x dfg/dx x (u+p)

            // complete sensitivity computation
            sens += const_fact;
            sens *= wgt_case[l]; // multiply by load case weight

            // add to overall sensitivity
            elSens[index2(eNum,d,numDual)] += sens; // add value for this load case
        }
    }
}

// compute Tri3 stress tensor at a point in the element
void CHak2DTri_3::stress(Coord p, double *elDisp, double *strs)
{
    // p must be non-dimensional

	double stn[3];

	// compute strain tensor (vector)
	strain(p, elDisp, stn); // NB: p should be in non-dim coordinates (i.e center = 0,0,0)

	// get material property matrix
    double *Em = Mat->getDsts();

	// multiply: stress = material matrix x strain
	cblas_dgemv(CblasRowMajor, CblasNoTrans, 3, 3, 1.0, Em, 3, stn, 1, 0.0, strs, 1);
}
