/*
 * CHak2DPlate3.cpp
 *
 *  Created on: 23 Feb 2015
 *      Author: knai20-admin
 */

#include "CHak2DPlate3.h"

CHak2DPlate_3::CHak2DPlate_3() {
	// TODO Auto-generated constructor stub
	// constructor sets num nodes
    numNode=3;
	nodeSpace();
    volume = -1.0;
    t = -1.0;
	Ke=0;
	Me=0;
}

CHak2DPlate_3::~CHak2DPlate_3() {
	// TODO Auto-generated destructor stub
	// destructor (delete memory)
	delMat();
}

// create DKT element stiffness matrix
void CHak2DPlate_3::makeKe()
{
    if(Ke){delete [] Ke;}
    Ke = new double [45];

	// read nodal coords (assume the plate lies in a 2D space for now)
    Coord p1,p2,p3;
	p1.x = gNodes[getGlobal(0)].x;
    p1.y = gNodes[getGlobal(0)].y;
	p2.x = gNodes[getGlobal(1)].x;
    p2.y = gNodes[getGlobal(1)].y;
    p3.x = gNodes[getGlobal(2)].x;
    p3.y = gNodes[getGlobal(2)].y;

    // Compute matrix!!
    dktKe_sub(Ke, p1, p2, p3, t, Mat);
}

// create DKT element consistent mass matrix
void CHak2DPlate_3::makeMe()
{
    if(Me){delete [] Me;}
    Me = new double [45];

	// read nodal coords (assume the plate lies in a 2D space for now)
    Coord p1,p2,p3;
	p1.x = gNodes[getGlobal(0)].x;
    p1.y = gNodes[getGlobal(0)].y;
	p2.x = gNodes[getGlobal(1)].x;
    p2.y = gNodes[getGlobal(1)].y;
    p3.x = gNodes[getGlobal(2)].x;
    p3.y = gNodes[getGlobal(2)].y;

    // Compute matrix!!
    dktME_sub(Me, p1, p2, p3, t, Mat);
}

// compute the element volume
void CHak2DPlate_3::makeVol()
{
    Coord3D p[3]; getNcrd(p);
    triangle tri(p);
    volume = tri.area() * t; // volume = area x thickness
}

void CHak2DPlate_3::getDisp(double *elDisp, double *globDisp)
{
	int i,j,nd;
	int ind = 0;

	for(i=0;i<3;i++) // for each node
	{
		nd = getGlobal(i)*3; // 1st dof
		for(j=0;j<3;j++)
		{
			elDisp[ind++] = globDisp[nd++];
		}
	}
}

// compute DKT strain tensor at a point in the element
void CHak2DPlate_3::strain(Coord p, double *elDisp, double *stn)
{
    int i;
    double dx;

    // initialize to zero
    for(i=0;i<3;i++)
    {
        stn[0] = 0.0;
    }

    // read nodal coords
    Coord3D np[3];
    Coord nCoord[3];
    getNcrd(np);
    for(i=0;i<3;i++)
    {
        nCoord[i].x = np[i].x;
        nCoord[i].y = np[i].y;
    }

    double tri_disp[9]; // need 9 displacements for stress calculation (3 dof per node)

    // compute strain displacement matrix
    double area;
	double Hxe[9];
	double Hye[9];
	double Hxn[9];
	double Hyn[9];
    double factors[16]; // geometric factors

    // compute common geometric factors for the triangle
    dkT_geomFact(nCoord, factors, &area);

    // convert point to non-dimensional coords
    Coord pnd;
    pnd.x = (p.x*factors[15] - factors[13]*p.y - nCoord[0].x*factors[15] + factors[13]*nCoord[0].y)/area;
	pnd.y = -(p.x*factors[14] - factors[12]*p.y - nCoord[0].x*factors[14] + factors[12]*nCoord[0].y)/area;

    // compute the H matrices at the input point
    dkt_Hmatrices(pnd,factors,Hxe,Hye,Hxn,Hyn);

    // compute strains at the input point
	double B[27]; // array for B matrix

    // strain displacement matrix
    for(i=0;i<9;i++) // for each column in B
    {
        B[i] = factors[15]*Hxe[i]-factors[14]*Hxn[i];		// first row
        B[i+9] = factors[12]*Hyn[i]-factors[13]*Hye[i];	// second row
        B[i+18] = factors[12]*Hxn[i]-factors[13]*Hxe[i]+factors[15]*Hye[i]-factors[14]*Hyn[i]; // third row
    }

    for(i=0;i<9;i++) // for each column in B
    {
        stn[0] += B[i] * tri_disp[i];
        stn[1] += B[i+9] * tri_disp[i];
        stn[2] += B[i+18] * tri_disp[i];
    }

    dx = 0.5*t/area; // distance from mid-plane to surfcae is t/2
    stn[0] *= dx;
    stn[1] *= dx;
    stn[2] *= dx;
}

void CHak2DPlate_3::stress(Coord p, double *elDisp, double *strs)
{
    int i;
    for(i=0;i<3;i++) // initilaize to zero
    {
        strs[i] = 0.0;
    }

    // get material property matrix
    double *D = Mat->getDsts();

    // compute strains
    double stn[3];
    strain(p, elDisp, stn);

    // complete multiplication: stress = D x strain
    for(i=0;i<3;i++)
    {
        strs[0] += D[i] * stn[i];
        strs[1] += D[i+3] * stn[i];
        strs[2] += D[i+6] * stn[i];
    }
}
