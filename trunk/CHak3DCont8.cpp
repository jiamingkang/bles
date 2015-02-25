/*
 * CHak3DCont8.cpp
 *
 *  Created on: 23 Feb 2015
 *      Author: knai20-admin
 */

#include "CHak3DCont8.h"

CHak3DCont_8::CHak3DCont_8() {
	// TODO Auto-generated constructor stub
	// constructor sets num nodes
	numNode=8;
	nodeSpace();
    volume = -1.0; // indicates volume not computed
	Ke=0;
	Me=0;
    dof_ind=0;
}

CHak3DCont_8::~CHak3DCont_8() {
	// TODO Auto-generated destructor stub
	// destructor (delete memory)
	delMat();
}

// make C8 element stiffness matrix
void CHak3DCont_8::makeKe()
{
	int i,j,k,inc; //,err;
	double sum, Jdet;
	Coord3D gp;
	Coord3D np[8];

	if(Ke){delete [] Ke;}
	Ke = new double [300]();

	// read in node coordinates
	getNcrd(np);

	// get material property matrix
    double *Em = Mat->getD();

	double B[144];	// Strain displacement matrix (6 x 24)
	double Bt[144];	// Auxillary matrix

	// use 2 point Gauss rule (8 points total)
	for(i=0;i<8;i++)
	{
		gp.x = gauss2 * c8_pos[i].x;
		gp.y = gauss2 * c8_pos[i].y;
		gp.z = gauss2 * c8_pos[i].z;

		// compute isoparametric B matrix (6 x 24)
		Jdet = C8_isoBMat(gp, np, B);

        if(Jdet < 1.0E-12 || Jdet > 1.0E+15)
        {
            std::cout << "\nERROR in makeKe for C8 element, Jdet = " << Jdet;
        }

		// first multiply B^T * D
		cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, 24, 6, 6, 1.0, B, 24, Em, 6, 0.0, Bt, 6);

		// complete multiplication (for upper triangle only)
		for(j=0;j<24;j++) // rows
		{
			for(k=j;k<24;k++) // columns (upper triangle)
			{
				sum = 0.0;
				for(inc=0;inc<6;inc++)
				{
					sum += Bt[index2(j,inc,6)]*B[index2(inc,k,24)];
				}
				Ke[tri_ind(24,j,k)] += sum*Jdet; //*weight;
			}
		}
	}
}

// make a C8 element (consistent) mass matrix
void CHak3DCont_8::makeMe()
{
	if(Me){delete [] Me;}
	Me = new double [300]();

	C8_makeMe(this);
}

// function to get element displacement array
void CHak3DCont_8::getDisp(double *elDisp, double *globDisp)
{
	int i,j,nd;
	int ind = 0;

	for(i=0;i<8;i++) // for each node
	{
		if(dof_ind){ nd = dof_ind[getGlobal(i)]; }
        else{ nd = getGlobal(i)*3; } // 1st dof
		for(j=0;j<3;j++)
		{
			elDisp[ind++] = globDisp[nd++];
		}
	}
}

// compute C8 strain tensor at a point in the element
void CHak3DCont_8::strain(Coord3D p, double *elDisp, double *stn)
{
	Coord3D np[8];
	double B[144];	// Strain displacement matrix

	// get element nodal coords
	getNcrd(np);

	// compute isoparametric B matrix (6 x 24)
	C8_isoBMat(p, np, B); // NB: p should be in non-dim coordinates (i.e center = 0,0,0)

	// multiply: strain = B x u
	cblas_dgemv(CblasRowMajor, CblasNoTrans, 6, 24, 1.0, B, 24, elDisp, 1, 0.0, stn, 1);
}

// compute C8 stress tensor at a point in the element
void CHak3DCont_8::stress(Coord3D p, double *elDisp, double *strs)
{
	double stn[6];

	// compute strain tensor (vector)
	strain(p, elDisp, stn); // NB: p should be in non-dim coordinates (i.e center = 0,0,0)

	// gwt material property matrix
    double *Em = Mat->getD();

	// multiply: stress = material matrix x strain
	cblas_dgemv(CblasRowMajor, CblasNoTrans, 6, 6, 1.0, Em, 6, stn, 1, 0.0, strs, 1);
}

void CHak3DCont_8::makeVol()
{
    if(volume < 0.0)
    {
        Coord3D np[8];
        // read in node coordinates
        getNcrd(np);
        // compute volume
        volume = Hex_Volume(np);
    }
}

// add self-weight loading to a global loading vector
void CHak3DCont_8::addWeight(Coord3D vec, double mag, int *dof_ind, double *gload)
{
    C8_addWeight(this, &vec, mag, dof_ind, gload);
}

// compute C8 sensitivities at all gauss points, for eigenvalue problems
void CHak3DCont_8::shpSens_Eig(int numDual, int numEig, double **disp_prim, double **disp_dual,
                      double *eig, double *wgt_case, double *gSens, int pNum, double alpha)
{
    int i,j,p,l,d,ind;
    double ftemp;
    Coord3D np[8], gp;
	double B[144];	// Strain displacement matrix
    double ud[24];  // Primary element displacement array
    double dd[24];  // Dual element displacement array
    double stn[6];  // strain tensor
    double strs[6];  // stress tensor
    double sens; // variable to keep running total of a sensitivity value
    double ux[8], uy[8], uz[8];
    double px[8], py[8], pz[8];
    double ug[3];

	// get element nodal coords
	getNcrd(np);

    // get material property matrix
    double *Em = Mat->getD();
    double rho = Mat->getDens();

    // for each gauss point
    for(p=0;p<8;p++)
    {
        gp.x = gauss2 * c8_pos[p].x;
		gp.y = gauss2 * c8_pos[p].y;
		gp.z = gauss2 * c8_pos[p].z; // non-dim gauss point coords

        // compute isoparametric B matrix (6 x 24) & Ba (6 x 9)
        C8_isoBMat(gp, np, B); // NB: gp should be in non-dim coordinates (i.e center = 0,0,0)

        // for each eigenvalue
        for(l=0;l<numEig;l++)
        {
            // get the primary element displacements
            getDisp(ud, disp_prim[l]);

            for(i=0;i<8;i++)
            {
                j=i*3;
                ux[i] = ud[j];
                uy[i] = ud[j+1];
                uz[i] = ud[j+2];
            }
            ug[0] = C8_interpolate(gp,ux); // x disp at Gauss point
            ug[1] = C8_interpolate(gp,uy); // x disp at Gauss point
            ug[2] = C8_interpolate(gp,uz); // x disp at Gauss point

            // compute stress tensor Ee(u) = strs
            // multiply: strain = B x u
            cblas_dgemv(CblasRowMajor, CblasNoTrans, 6, 24, 1.0, B, 24, ud, 1, 0.0, stn, 1);

            // Finally multiply: stress = material matrix x strain
            cblas_dgemv(CblasRowMajor, CblasNoTrans, 6, 6, 1.0, Em, 6, stn, 1, 0.0, strs, 1);

            // for each dual response p
            for(d=0;d<numDual;d++)
            {
                // add the 3 componnents (multipled by the load case weight)
                // Ee(u)e(p) - fg(u+p) - const_fact

                // get the dual element displacements
                getDisp(dd, disp_dual[index2(l,d,numDual)]);

                // compute dual strain tensor e(p) = stn
                // multiply: strain = B x u
                cblas_dgemv(CblasRowMajor, CblasNoTrans, 6, 24, 1.0, B, 24, dd, 1, 0.0, stn, 1);

                // compute stress(u) x strain(p)
                sens = 0.0; // reset to zero
                for(i=0;i<6;i++)
                {
                    sens += strs[i]*stn[i];
                }
                sens *= alpha ; // need to multiply by the volume ratio (account for modified modulus)

                // interpolate, then multiply
                // need to extract u, v, w disp and forces in seperate vectors
                for(i=0;i<8;i++)
                {
                    j=i*3;
                    px[i] = dd[j];
                    py[i] = dd[j+1];
                    pz[i] = dd[j+2];
                }

                // now interpolate and add to sens
                ftemp =  ug[0]*C8_interpolate(gp,px) + ug[1]*C8_interpolate(gp,py) + ug[2]*C8_interpolate(gp,pz);
                sens -= eig[l]*ftemp*rho*alpha;

                // add to overall sensitivity
                ind=pNum+p;
                gSens[index2(ind,d,numDual)] += sens*wgt_case[l]; // add value for this load case
            }
        }
    }
}

// compute C8 sensitivities at all gauss points
void CHak3DCont_8::shpSens(int numDual, int numCase, double **disp_prim, double **disp_dual,
                  double *wgt_fact, double *wgt_case, double *gSens, int pNum, double *Nf, double alpha)
{
    int i,j,p,l,d,ind;
    double const_fact, ftemp;
    Coord3D np[8], gp;
	double B[144];	// Strain displacement matrix
    double Ba[54];  // Strain displacement matrix for internal dof
    double ud[24];  // Primary element displacement array
    double dd[24];  // Dual element displacement array
    double stn[6];  // strain tensor
    double strs[6];  // stress tensor
    double sens; // variable to keep running total of a sensitivity value
    double /*ux[8], uy[8],*/ uz[8];
    double /*px[8], py[8],*/ pz[8];
    double sw_fact;

	// get element nodal coords
	getNcrd(np);

    // compute material property matrix
    double *Em = Mat->getD();

    // for each gauss point
    for(p=0;p<8;p++)
    {
        gp.x = gauss2 * c8_pos[p].x;
		gp.y = gauss2 * c8_pos[p].y;
		gp.z = gauss2 * c8_pos[p].z; // non-dim gauss point coords

        // compute isoparametric B matrix (6 x 24) & Ba (6 x 9)
        C8_isoBMat(gp, np, B); // NB: gp should be in non-dim coordinates (i.e center = 0,0,0)
        C8M_isoBMat(gp, np, Ba);

        // for each load case
        for(l=0;l<numCase;l++)
        {
            // compute the constant sensitivty factor for the element
            const_fact = wgt_fact[l] * Mat->getDens(); // DO NOT multiply by the volume ratio

            // get the primary element displacements
            getDisp(ud, disp_prim[l]);

            for(i=0;i<8;i++)
            {
                j=i*3;
                //ux[i] = ud[j];
                //uy[i] = ud[j+1];
                uz[i] = ud[j+2];
            }

            // compute self weight factor fg * u
            ftemp = Mat->getDens() * Nf[l] * -9.81 * alpha; // mass times acceleration (in +ve z direction)
            sw_fact = ftemp * C8_interpolate(gp,uz); //( C8_interpolate(gp,ux) + C8_interpolate(gp,uy) + C8_interpolate(gp,uz) );;

            // compute stress tensor Ee(u) = strs
            // multiply: strain = B x u
            cblas_dgemv(CblasRowMajor, CblasNoTrans, 6, 24, 1.0, B, 24, ud, 1, 0.0, stn, 1);

            // Finally multiply: stress = material matrix x strain
            cblas_dgemv(CblasRowMajor, CblasNoTrans, 6, 6, 1.0, Em, 6, stn, 1, 0.0, strs, 1);

            // for each dual response p
            for(d=0;d<numDual;d++)
            {
                // add the 3 componnents (multipled by the load case weight)
                // Ee(u)e(p) - fg(u+p) - const_fact

                // get the dual element displacements
                getDisp(dd, disp_dual[index2(l,d,numDual)]);

                // compute dual strain tensor e(p) = stn
                // multiply: strain = B x u
                cblas_dgemv(CblasRowMajor, CblasNoTrans, 6, 24, 1.0, B, 24, dd, 1, 0.0, stn, 1);

                // compute stress(u) x strain(p)
                sens = 0.0; // reset to zero
                for(i=0;i<6;i++)
                {
                    sens += strs[i]*stn[i];
                }
                sens *= alpha ; // need to multiply by the volume ratio (account for modified modulus)

                // interpolate, then multiply
                // need to extract u, v, w disp and forces in seperate vectors
                for(i=0;i<8;i++)
                {
                    j=i*3;
                    //px[i] = dd[j];
                    //py[i] = dd[j+1];
                    pz[i] = dd[j+2];
                }

                // now interpolate and add to sens
                sens -= sw_fact;
                sens -= ftemp * C8_interpolate(gp,pz); //( C8_interpolate(gp,px) + C8_interpolate(gp,py) + C8_interpolate(gp,pz) );
                sens -= const_fact;
                sens *= wgt_case[l]; // multiply by load case weight

                // add to overall sensitivity
                ind=pNum+p;
                gSens[index2(ind,d,numDual)] += sens; // add value for this load case
            }
        }
    }
}
