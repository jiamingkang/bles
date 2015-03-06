/*
 * CHak2DPlate4.cpp
 *
 *  Created on: 23 Feb 2015
 *      Author: knai20-admin
 */

#include "CHak2DPlate4.h"

CHak2DPlate_4::CHak2DPlate_4() {
	// TODO Auto-generated constructor stub
	// constructor sets num nodes
    numNode=4;
	nodeSpace();
    volume = -1.0;
    t = -1.0;
	Ke=0;
	Me=0;

}

CHak2DPlate_4::~CHak2DPlate_4() {
	// TODO Auto-generated destructor stub
	// destructor (delete memory)
	delMat();
}

// create a DKT2 element matrix
// NB: a quad is made from two tri elems across shortest (or specified) diagonal
void CHak2DPlate_4::makeKe()
{
	int i;
    if(Ke){delete [] Ke;}
    Ke = new double [78];
	for(i=0;i<78;i++)
	{
		Ke[i] = 0.0; // initialize to zero
	}

	double tri[45];	// pointer for a tri Ke matrices
	int nt[3];  // array of global node nums wrt local nodes

	// read nodal coords (assume the plate lies in a 2D space for now)
    Coord p1,p2,p3,p4;
	p1.x = gNodes[getGlobal(0)].x;
    p1.y = gNodes[getGlobal(0)].y;
	p2.x = gNodes[getGlobal(1)].x;
    p2.y = gNodes[getGlobal(1)].y;
    p3.x = gNodes[getGlobal(2)].x;
    p3.y = gNodes[getGlobal(2)].y;
    p4.x = gNodes[getGlobal(3)].x;
    p4.y = gNodes[getGlobal(3)].y;

	// compute length of diagonals
	double L13,L24;
	if( (diag != 1) && (diag != 2) && (diag !=4) )
	{
		L13 = dist2D(p1, p3);
		L24 = dist2D(p2, p4);
		diag = (L24 > L13) ? 2 : 1; // decide which diagonal to use
	}

	if(diag==2 || diag==4)
	{
		// first triangle is 1,2,4
		nt[0] = 0;
		nt[1] = 1;
		nt[2] = 3;

		// Compute matrix!!
		// input node coords material property class and thickness
		dktKe_sub(tri, p1, p2, p4, t, Mat);

		// add returned matrix to Ke
		addMatrix(3, 4, tri, Ke, 3, nt);

		// second triangle is 2,3,4
		nt[0] = 1;
		nt[1] = 2;

		// Compute matrix!!
		dktKe_sub(tri, p2, p3, p4, t, Mat);

		// add returned matrix to Ke
		addMatrix(3, 4, tri, Ke, 3, nt);
	}
	else if(diag==1 || diag==4)
	{
		// first triangle is 1,2,3
		nt[0] = 0;
		nt[1] = 1;
		nt[2] = 2;

		// Compute matrix!!
		dktKe_sub(tri, p1, p2, p3, t, Mat);

		// add returned matrix to Ke
		addMatrix(3, 4, tri, Ke, 3, nt);

		// second triangle is 1,3,4
		nt[1] = 2;
		nt[2] = 3;

		// Compute matrix!!
		dktKe_sub(tri, p1, p3, p4, t, Mat);

		// add returned matrix to Ke
		addMatrix(3, 4, tri, Ke, 3, nt);
	}
}

void CHak2DPlate_4::makeMe()
{
    int i;
    if(Me){delete [] Me;}
    Me = new double [78];
	for(i=0;i<78;i++)
	{
		Me[i] = 0.0; // initialize to zero
	}

	double tri[45];			// pointer for a tri Ke matrices
	int nt[3];

	// read nodal coords (assume the plate lies in a 2D space for now)
    Coord p1,p2,p3,p4;
	p1.x = gNodes[getGlobal(0)].x;
    p1.y = gNodes[getGlobal(0)].y;
	p2.x = gNodes[getGlobal(1)].x;
    p2.y = gNodes[getGlobal(1)].y;
    p3.x = gNodes[getGlobal(2)].x;
    p3.y = gNodes[getGlobal(2)].y;
    p4.x = gNodes[getGlobal(3)].x;
    p4.y = gNodes[getGlobal(3)].y;

	// compute length of diagonals
	double L13,L24;
	if( (diag != 1) && (diag != 2) && (diag !=4) )
	{
		L13 = dist2D(p1, p3);
		L24 = dist2D(p2, p4);
		diag = (L24 > L13) ? 2 : 1; // decide which diagonal to use
	}

	if(diag==2 || diag==4)
	{
		// first triangle is 1,2,4
		nt[0] = 0;
		nt[1] = 1;
		nt[2] = 3;

		// Compute matrix!!
		// input node coords density and thickness
		dktME_sub(tri, p1, p2, p4, t, Mat);

		// add returned matrix to Me
		addMatrix(3, 4, tri, Me, 3, nt);

		// second triangle is 2,3,4
		nt[0] = 1;
		nt[1] = 2;

		// Compute matrix!!
		dktME_sub(tri, p2, p3, p4, t, Mat);

		// add returned matrix to Me
		addMatrix(3, 4, tri, Me, 3, nt);
	}
	else if(diag==1 || diag==4)
	{
		// first triangle is 1,2,3
		nt[0] = 0;
		nt[1] = 1;
		nt[2] = 2;

		// Compute matrix!!
		dktME_sub(tri, p1, p2, p3, t, Mat);

		// add returned matrix to Me
		addMatrix(3, 4, tri, Me, 3, nt);

		// second triangle is 1,3,4
		nt[1] = 2;
		nt[2] = 3;

		// Compute matrix!!
		dktME_sub(tri, p1, p3, p4, t, Mat);

		// add returned matrix to Me
		addMatrix(3, 4, tri, Me, 3, nt);
	}
}

// compute the element volume
void CHak2DPlate_4::makeVol()
{
    Coord3D p[3]; getNcrd(p);
    triangle tri(p);
    volume = tri.area() * t; // volume = area x thickness
}

// compute DKT strain tensor at a point in the element
void CHak2DPlate_4::strain(Coord p, double *elDisp, double *stn)
{
    int i,temp;
    double dx;

    // initialize to zero
    for(i=0;i<3;i++)
    {
        stn[0] = 0.0;
    }

    // read nodal coords
    Coord3D np[4];
    Coord nCoord[4];
    getNcrd(np);
    for(i=0;i<4;i++)
    {
        nCoord[i].x = np[i].x;
        nCoord[i].y = np[i].y;
    }

    double tri_disp[9]; // need 9 displacements for stress calculation (3 dof per node)

    // work out which half of the DKT the point lies in
    int tri = findTri(p);

    if(diag==1)
    {
        for(i=0;i<3;i++)
        {
            tri_disp[i] = elDisp[i]; // always use node 1
        }

        if(tri==1)// use triangle 1,2,3
        {
            for(i=3;i<9;i++)
            {
                tri_disp[i] = elDisp[i];
            }
        }

        else // otherwise use triangle 1,3,4
        {
            nCoord[1] = nCoord[2];
            nCoord[2] = nCoord[3];

            for(i=3;i<9;i++)
            {
                temp = i+3; // skip node 2
                tri_disp[i] = elDisp[temp];
            }
        }
    }
    else if(diag==2)
    {
        for(i=6;i<9;i++)
        {
            temp=i+3;
            tri_disp[i] = elDisp[temp]; // always use node 4
        }

        if(tri==1) // use triangle 1,2,4
        {
            for(i=0;i<6;i++)
            {
                tri_disp[i] = elDisp[i];
            }
            nCoord[2] = nCoord[3];
        }
        else // otherwise use triangle 2,3,4
        {
            for(i=0;i<6;i++)
            {
                temp = i+3;
                tri_disp[i] = elDisp[temp];
            }
            nCoord[0] = nCoord[1];
            nCoord[1] = nCoord[2];
        }
    }

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
        B[i+9] = factors[12]*Hyn[i]-factors[13]*Hye[i];     // second row
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

void CHak2DPlate_4::stress(Coord p, double *elDisp, double *strs)
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

// fucntion to work out which triangle a point lies in
int CHak2DPlate_4::findTri(Coord p)
{
    double dx,dy,Lxy,y_rel;
    int tri;
    y_rel=0.0;
    if(diag == 1)
    {
        dx = gNodes[getGlobal(2)].x - gNodes[getGlobal(0)].x;
        dy = gNodes[getGlobal(2)].y - gNodes[getGlobal(0)].y;
        Lxy = sqrt(dx*dx + dy*dy);

        // turn dx * dy into cos and sin, respectively
        dx /= Lxy; // cos
        dy /= Lxy; // sin

        // compute y coord of node realtive to side as the x-axis
        y_rel = (p.y-gNodes[getGlobal(0)].y)*dx-(p.x-gNodes[getGlobal(0)].x)*dy;
    }

    else if(diag == 2)
    {
        dx = gNodes[getGlobal(3)].x - gNodes[getGlobal(1)].x;
        dy = gNodes[getGlobal(3)].y - gNodes[getGlobal(1)].y;
        Lxy = sqrt(dx*dx + dy*dy);

        // turn dx * dy into cos and sin, respectively
        dx /= Lxy; // cos
        dy /= Lxy; // sin

        // compute y coord of node realtive to side as the x-axis
        y_rel = (p.x-gNodes[getGlobal(1)].x)*dy - (p.y-gNodes[getGlobal(1)].y)*dx;
    }

    tri = (y_rel < 0.0) ? 1 : 2;  // pick which triangle to use
    return tri;
}

void CHak2DPlate_4::getDisp(double *elDisp, double *globDisp)
{
	int i,j,nd;
	int ind = 0;

	for(i=0;i<4;i++) // for each node
	{
		nd = getGlobal(i)*3; // 1st dof
		for(j=0;j<3;j++)
		{
			elDisp[ind++] = globDisp[nd++];
		}
	}
}
