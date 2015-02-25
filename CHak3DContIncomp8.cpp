/*
 * CHak3DContIncomp8.cpp
 *
 *  Created on: 23 Feb 2015
 *      Author: knai20-admin
 */

#include "CHak3DContIncomp8.h"

CHak3DContIncomp_8::CHak3DContIncomp_8() {
	// TODO Auto-generated constructor stub
	// constructor sets num nodes
	numNode=8;
	nodeSpace();
    volume = -1.0; // indicates volume not computed
	Ke=0;
	Me=0;
}

CHak3DContIncomp_8::~CHak3DContIncomp_8() {
	// TODO Auto-generated destructor stub
	// destructor (delete memory)
	delMat();
}


// make C8M element stiffness matrix
void C8M::makeKe()
{
	int i,j,k,inc; //,err;
	double sum, Jdet, Jdet2;
	Coord3D gp;
	Coord3D np[8];
	double Krc[216]; // (9 x 24)
	double Kcc[81];  // (9 x 9) auxillary stiffness matrices

    if(Ke){delete [] Ke;}
	Ke = new double [300]();

	if(Km){delete [] Km;}
	Km = new double [216]();

	// initialize matrices to zero
    for(i=0;i<216;i++) {
		Krc[i] = 0.0;
	}
	for(i=0;i<81;i++) {
		Kcc[i] = 0.0;
	}

	// read in node coordinates
	getNcrd(np);

	// get material property matrix
    double *Em = Mat->getD();

	double B[144];	// Strain displacement matrix
	double Ba[54];   // Strain displacement matrix for nodelss dof
	double Bt[144];	// Auxillary matrix

	// use 2 point Gauss rule (8 points total)
	for(i=0;i<8;i++)
	{
		gp.x = gauss2 * c8_pos[i].x;
		gp.y = gauss2 * c8_pos[i].y;
		gp.z = gauss2 * c8_pos[i].z;

		// compute isoparametric B matrix
		Jdet = C8_isoBMat(gp, np, B);

		// compute isoparmetric nodeless B matrix (6 x 9)
		Jdet2 = C8M_isoBMat(gp, np, Ba);

		cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, 24, 6, 6, 1.0, B, 24, Em, 6, 0.0, Bt, 6);

		// complete multiplication for Krr matrix (for upper triangle only)
		for(j=0;j<24;j++) // rows
		{
			for(k=j;k<24;k++) // columns (upper triangle)
			{
				sum = 0.0;
				for(inc=0;inc<6;inc++)
				{
					sum += Bt[index2(j,inc,6)]*B[index2(inc,k,24)];
				}
				Ke[tri_ind(24,j,k)] += Jdet * sum;
			}
		}

        // complete multiplication for Krc matrix
		for(j=0;j<24;j++) // rows
		{
			for(k=0;k<9;k++) // columns
			{
				sum = 0.0;
				for(inc=0;inc<6;inc++)
				{
					sum += Bt[index2(j,inc,6)]*Ba[index2(inc,k,9)];
				}
				Krc[index2(j,k,9)] += Jdet * sum; // store column major ordered for FORTRAN subroutine (or transposed)
			}
		}

		// compute Ba^T x Em
		cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, 9, 6, 6, 1.0, Ba, 9, Em, 6, 0.0, Bt, 6);

		// complete multiplication for Kcc matrix
		for(j=0;j<9;j++) // rows
		{
			for(k=j;k<9;k++) // columns (only need upper triangle)
			{
				sum = 0.0;
				for(inc=0;inc<6;inc++)
				{
					sum += Bt[index2(j,inc,6)]*Ba[index2(inc,k,9)];
				}
				Kcc[index2(j,k,9)] += Jdet2 * sum;
			}
		}
	}

	// solve Kcc^-1 x -Kcr = Km (9 x 24)
	// copy Km = -Kcr^T
	for(i=0;i<216;i++)
	{
		Km[i] = -Krc[i];
	}

	static char uplo = 'L';
	i = 9;
	j = 24;
	int *ipiv = new int [9];
	double *work = new double [9];
	dsysv_(&uplo, &i, &j, Kcc, &i, ipiv, Km, &i, work, &i, &k);
	//LAPACKE_dsysv(9, uplo, 24, 1, Kcc, 9, ipiv, Km, 9);
    if(k!=0)
    {
        std::cout << "\nError in dsysv info = " << k;
    }
	// NB Fortran routine will return a column major ordered matrix
	// Km = (-Kcc^-1 x Kcr)^T in row major order
	delete [] ipiv;
	delete [] work;

	// multiply Kcr^T x Km = Ka (24 x 24)
	// and modfiy Ke + Ka
	for(j=0;j<24;j++) // rows
	{
		for(k=j;k<24;k++) // columns (upper triangle)
		{
			sum = 0.0;
			for(inc=0;inc<9;inc++)
			{
				sum += Km[index2(k,inc,9)]*Krc[index2(j,inc,9)];
			}
			Ke[tri_ind(24,j,k)] += sum;
		}
	}
}

void C8M::makeMe()
{
	if(Me){delete [] Me;}
	Me = new double [300]();

	C8_makeMe(this);
}

// function to get element displacement array
void C8M::getDisp(double *elDisp, double *globDisp)
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

// function to get element complex displacement array
void C8M::getDispC(dcompLPK *elDisp, dcompLPK *globDisp)
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

// compute C8M strain tensor at a point in the element
void C8M::strain(Coord3D p, double *elDisp, double *stn)
{
	Coord3D np[8];
	double B[144];	// Strain displacement matrix
    double Ba[54];  // Strain displacement matrix for internal dof

	// get element nodal coords
	getNcrd(np);

	// compute isoparametric B matrix (6 x 24) & Ba (6 x 9)
	C8_isoBMat(p, np, B); // NB: p should be in non-dim coordinates (i.e center = 0,0,0)
    C8M_isoBMat(p, np, Ba);

	// multiply: strain = B x u
	cblas_dgemv(CblasRowMajor, CblasNoTrans, 6, 24, 1.0, B, 24, elDisp, 1, 0.0, stn, 1);

    double stn_a[6]; // strain tensor for internal dof
    double Ua[9]; // displacements for internal dof

    // multiply: Ua = Km x u
	cblas_dgemv(CblasRowMajor, CblasTrans, 24, 9, 1.0, Km, 9, elDisp, 1, 0.0, Ua, 1);

    // multiply: stn_a = Ba x Ua
	cblas_dgemv(CblasRowMajor, CblasNoTrans, 6, 9, 1.0, Ba, 9, Ua, 1, 0.0, stn_a, 1);

    // add the two strain tensors
    int i;
    for(i=0;i<6;i++)
    {
        stn[i] += stn_a[i];
    }
}

// compute C8M stress tensor at a point in the element
void C8M::stress(Coord3D p, double *elDisp, double *strs)
{
	double stn[6];

	// compute strain tensor (vector)
	strain(p, elDisp, stn); // NB: p should be in non-dim coordinates (i.e center = 0,0,0)

	// get material property matrix
    double *Em = Mat->getD();

	// multiply: stress = material matrix x strain
	cblas_dgemv(CblasRowMajor, CblasNoTrans, 6, 6, 1.0, Em, 6, stn, 1, 0.0, strs, 1);
}

void C8M::makeVol()
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
void C8M::addWeight(Coord3D vec, double mag, int *dof_ind, double *gload)
{
    C8_addWeight(this, &vec, mag, dof_ind, gload);
}


// compute C8M sensitivities at all gauss points, for eigenvalue problems
void C8M::shpSens_Eig(int numDual, int numEig, double **disp_prim, double **disp_dual,
                    double *eig, double *wgt_case, double *gSens, int pNum, double alpha)
{
    int i,j,p,l,d,ind;
    double ftemp;
    Coord3D np[8], gp, gp2;
	double B[144];	// Strain displacement matrix
    double Ba[54];  // Strain displacement matrix for internal dof
    double ud[24];  // Primary element displacement array
    double dd[24];  // Dual element displacement array
    double stn[6];  // strain tensor
    double strs[6];  // stress tensor
    double stn_a[6]; // strain tensor for internal dof
    double Ua[9];  // displacements for internal dof
    double sens; // variable to keep running total of a sensitivity value
    double ux[8], uy[8], uz[8];
    double px[8], py[8], pz[8];
    double ug[3];

	// get element nodal coords
	getNcrd(np);

    // get material property matrix
    double *Em = Mat->getD();
    double rho = Mat->getDens();
    double Malpha = (alpha < 1.01E-3) ? 1.0E-5 : alpha; // extra penalization for void mass

    // for each gauss point
    for(p=0;p<8;p++)
    {
        gp.x = gauss2 * c8_pos[p].x;
		gp.y = gauss2 * c8_pos[p].y;
		gp.z = gauss2 * c8_pos[p].z; // non-dim gauss point coords
        gp2.x = gp.x*gp.x; gp2.y = gp.y*gp.y; gp2.z = gp.z*gp.z;

        // compute isoparametric B matrix (6 x 24) & Ba (6 x 9)
        C8_isoBMat(gp, np, B); // NB: gp should be in non-dim coordinates (i.e center = 0,0,0)
        C8M_isoBMat(gp, np, Ba);

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

            // compute stress tensor Ee(u) = strs
            // multiply: strain = B x u
            cblas_dgemv(CblasRowMajor, CblasNoTrans, 6, 24, 1.0, B, 24, ud, 1, 0.0, stn, 1);
            // multiply: Ua = Km x u
            cblas_dgemv(CblasRowMajor, CblasTrans, 24, 9, 1.0, Km, 9, ud, 1, 0.0, Ua, 1);
            // multiply: stn_a = Ba x Ua
            cblas_dgemv(CblasRowMajor, CblasNoTrans, 6, 9, 1.0, Ba, 9, Ua, 1, 0.0, stn_a, 1);
            // add the two strain tensors;
            for(i=0;i<6;i++) { stn[i] += stn_a[i]; }
            // Finally multiply: stress = material matrix x strain
            cblas_dgemv(CblasRowMajor, CblasNoTrans, 6, 6, 1.0, Em, 6, stn, 1, 0.0, strs, 1);

            ug[0] = C8_interpolate(gp,ux) + (1.0-gp2.x)*Ua[0] + (1.0-gp2.y)*Ua[1] + (1.0-gp2.z)*Ua[2]; // x disp at Gauss point
            ug[1] = C8_interpolate(gp,uy) + (1.0-gp2.x)*Ua[3] + (1.0-gp2.y)*Ua[4] + (1.0-gp2.z)*Ua[5]; // y disp at Gauss point
            ug[2] = C8_interpolate(gp,uz) + (1.0-gp2.x)*Ua[6] + (1.0-gp2.y)*Ua[7] + (1.0-gp2.z)*Ua[8]; // z disp at Gauss point

            // for each dual response p
            for(d=0;d<numDual;d++)
            {
                // get the dual element displacements
                getDisp(dd, disp_dual[index2(l,d,numDual)]);

                // compute dual strain tensor e(p) = stn
                // multiply: strain = B x u
                cblas_dgemv(CblasRowMajor, CblasNoTrans, 6, 24, 1.0, B, 24, dd, 1, 0.0, stn, 1);
                // multiply: Ua = Km x u
                cblas_dgemv(CblasRowMajor, CblasTrans, 24, 9, 1.0, Km, 9, dd, 1, 0.0, Ua, 1);
                // multiply: stn_a = Ba x Ua
                cblas_dgemv(CblasRowMajor, CblasNoTrans, 6, 9, 1.0, Ba, 9, Ua, 1, 0.0, stn_a, 1);
                // add the two strain tensors;
                for(i=0;i<6;i++) { stn[i] += stn_a[i]; }

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
                ftemp =  ug[0]*(C8_interpolate(gp,px) + (1.0-gp2.x)*Ua[0] + (1.0-gp2.y)*Ua[1] + (1.0-gp2.z)*Ua[2]);
                ftemp += ug[1]*(C8_interpolate(gp,py) + (1.0-gp2.x)*Ua[3] + (1.0-gp2.y)*Ua[4] + (1.0-gp2.z)*Ua[5]);
                ftemp += ug[2]*(C8_interpolate(gp,pz) + (1.0-gp2.x)*Ua[6] + (1.0-gp2.y)*Ua[7] + (1.0-gp2.z)*Ua[8]);
                sens -= eig[l]*ftemp*rho*Malpha;

                // add to overall sensitivity
                ind=pNum+p;
                gSens[index2(ind,l,numEig)] += sens*wgt_case[l]; // add value for this load case
            }
        }
    }

}

// compute C8M sensitivities at all gauss points, for complex eigenvalue problems
void C8M::shpSens_CEig(int numDual, int numEig, dcompLPK **disp_prim, dcompLPK **disp_dual,
                      dcompLPK *eig, dcompLPK *wgt_case, double *gSens, int pNum, double K_fact, double alpha)
{
    int i,j,p,ind;
    double ftemp;
    Coord3D np[8], gp, gp2;
	double B[144];	// Strain displacement matrix
    double Ba[54];  // Strain displacement matrix for internal dof
    dcompLPK ud[24], dd[24], stemp, ctemp, mtemp;
    double ud_real[24]; double ud_imag[24];  // Primary element displacement array
    double dd_real[24]; double dd_imag[24];  // Dual element displacement array
    double stn_real[6]; double stn_imag[6];  // strain tensor
    double strs_real[6]; double strs_imag[6];  // stress tensor
    double stn_a_real[6]; double stn_a_imag[6]; // strain tensor for internal dof
    double Ua_real[9]; double Ua_imag[9];  // displacements for internal dof
    double sens; // variable to keep running total of a sensitivity value
    double ux_real[8], uy_real[8], uz_real[8];
    double ux_imag[8], uy_imag[8], uz_imag[8];
    double px_real[8], py_real[8], pz_real[8];
    double px_imag[8], py_imag[8], pz_imag[8];
    double ug_real[3], ug_imag[3];
    double dg_real[3], dg_imag[3];

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
        gp2.x = gp.x*gp.x; gp2.y = gp.y*gp.y; gp2.z = gp.z*gp.z;

        // compute isoparametric B matrix (6 x 24) & Ba (6 x 9)
        C8_isoBMat(gp, np, B); // NB: gp should be in non-dim coordinates (i.e center = 0,0,0)
        C8M_isoBMat(gp, np, Ba);

        // -- Stiffness sensitvity -- //

        // get 1st right hand eigenvector
        getDispC(ud, disp_prim[0]);
        for(i=0;i<24;i++)
        {
            ud_real[i] = ud[i].r; ud_imag[i] = ud[i].i;
        }

        // compute stress tensor Ee(u) = strs
        // multiply: strain = B x u
        cblas_dgemv(CblasRowMajor, CblasNoTrans, 6, 24, 1.0, B, 24, ud_real, 1, 0.0, stn_real, 1);
        cblas_dgemv(CblasRowMajor, CblasNoTrans, 6, 24, 1.0, B, 24, ud_imag, 1, 0.0, stn_imag, 1);
        // multiply: Ua = Km x u
        cblas_dgemv(CblasRowMajor, CblasTrans, 24, 9, 1.0, Km, 9, ud_real, 1, 0.0, Ua_real, 1);
        cblas_dgemv(CblasRowMajor, CblasTrans, 24, 9, 1.0, Km, 9, ud_imag, 1, 0.0, Ua_imag, 1);
        // multiply: stn_a = Ba x Ua
        cblas_dgemv(CblasRowMajor, CblasNoTrans, 6, 9, 1.0, Ba, 9, Ua_real, 1, 0.0, stn_a_real, 1);
        cblas_dgemv(CblasRowMajor, CblasNoTrans, 6, 9, 1.0, Ba, 9, Ua_imag, 1, 0.0, stn_a_imag, 1);
        // add the two strain tensors;
        for(i=0;i<6;i++) { stn_real[i] += stn_a_real[i];
                           stn_imag[i] += stn_a_imag[i];}
        // Finally multiply: stress = material matrix x strain
        cblas_dgemv(CblasRowMajor, CblasNoTrans, 6, 6, 1.0, Em, 6, stn_real, 1, 0.0, strs_real, 1);
        cblas_dgemv(CblasRowMajor, CblasNoTrans, 6, 6, 1.0, Em, 6, stn_imag, 1, 0.0, strs_imag, 1);

        // get 2nd left hand eigenvector
        getDispC(dd, disp_dual[1]);
        for(i=0;i<24;i++)
        {
            dd_real[i] = dd[i].r; dd_imag[i] = -dd[i].i; // complex conjugate
        }

        // compute dual strain tensor e(p) = stn
        // multiply: strain = B x u
        cblas_dgemv(CblasRowMajor, CblasNoTrans, 6, 24, 1.0, B, 24, dd_real, 1, 0.0, stn_real, 1);
        cblas_dgemv(CblasRowMajor, CblasNoTrans, 6, 24, 1.0, B, 24, dd_imag, 1, 0.0, stn_imag, 1);
        // multiply: Ua = Km x u
        cblas_dgemv(CblasRowMajor, CblasTrans, 24, 9, 1.0, Km, 9, dd_real, 1, 0.0, Ua_real, 1);
        cblas_dgemv(CblasRowMajor, CblasTrans, 24, 9, 1.0, Km, 9, dd_imag, 1, 0.0, Ua_imag, 1);
        // multiply: stn_a = Ba x Ua
        cblas_dgemv(CblasRowMajor, CblasNoTrans, 6, 9, 1.0, Ba, 9, Ua_real, 1, 0.0, stn_a_real, 1);
        cblas_dgemv(CblasRowMajor, CblasNoTrans, 6, 9, 1.0, Ba, 9, Ua_imag, 1, 0.0, stn_a_imag, 1);
        // add the two strain tensors;
        for(i=0;i<6;i++) { stn_real[i] += stn_a_real[i];
                            stn_imag[i] += stn_a_imag[i];}

        // compute Real & Imaginary part of stress(u) x strain(p) - stiffness sensitivity
        stemp.r = 0.0; stemp.i = 0.0; // reset to zero
        for(i=0;i<6;i++)
        {
            stemp.r += strs_real[i]*stn_real[i] - strs_imag[i]*stn_imag[i];
            stemp.i += strs_real[i]*stn_imag[i] + strs_imag[i]*stn_real[i];
        }
        stemp.r = stemp.r * alpha; // * K_fact;
        stemp.i = stemp.i * alpha; // * K_fact; // need to multiply by the volume ratio (account for modified modulus)

        // -- mass sensitvity -- //

        // interpolate, then multiply
        // need to extract u, v, w disp and forces in seperate vectors
        for(i=0;i<8;i++)
        {
            j=i*3;
            px_real[i] = dd[j].r;   px_imag[i] = -dd[j].i;
            py_real[i] = dd[j+1].r; py_imag[i] = -dd[j+1].i;
            pz_real[i] = dd[j+2].r; pz_imag[i] = -dd[j+2].i; // complex conjugate
        }

        // now interpolate
        dg_real[0] = C8_interpolate(gp,px_real) + (1.0-gp2.x)*Ua_real[0] + (1.0-gp2.y)*Ua_real[1] + (1.0-gp2.z)*Ua_real[2]; // x disp at Gauss point
        dg_imag[0] = C8_interpolate(gp,px_imag) + (1.0-gp2.x)*Ua_imag[0] + (1.0-gp2.y)*Ua_imag[1] + (1.0-gp2.z)*Ua_imag[2]; // x disp at Gauss point
        dg_real[1] = C8_interpolate(gp,py_real) + (1.0-gp2.x)*Ua_real[3] + (1.0-gp2.y)*Ua_real[4] + (1.0-gp2.z)*Ua_real[5]; // y disp at Gauss point
        dg_imag[1] = C8_interpolate(gp,py_imag) + (1.0-gp2.x)*Ua_imag[3] + (1.0-gp2.y)*Ua_imag[4] + (1.0-gp2.z)*Ua_imag[5]; // y disp at Gauss point
        dg_real[2] = C8_interpolate(gp,pz_real) + (1.0-gp2.x)*Ua_real[6] + (1.0-gp2.y)*Ua_real[7] + (1.0-gp2.z)*Ua_real[8]; // z disp at Gauss point
        dg_imag[2] = C8_interpolate(gp,pz_imag) + (1.0-gp2.x)*Ua_imag[6] + (1.0-gp2.y)*Ua_imag[7] + (1.0-gp2.z)*Ua_imag[8]; // z disp at Gauss point

        // get 2nd right hand eigenvector
        getDispC(ud, disp_prim[1]);
        for(i=0;i<24;i++)
        {
            ud_real[i] = ud[i].r; ud_imag[i] = ud[i].i;
        }

        for(i=0;i<8;i++)
        {
            j=i*3;
            ux_real[i] = ud[j].r;   ux_imag[i] = ud[j].i;
            uy_real[i] = ud[j+1].r; uy_imag[i] = ud[j+1].i;
            uz_real[i] = ud[j+2].r; uz_imag[i] = ud[j+2].i;
        }

        cblas_dgemv(CblasRowMajor, CblasTrans, 24, 9, 1.0, Km, 9, ud_real, 1, 0.0, Ua_real, 1);
        cblas_dgemv(CblasRowMajor, CblasTrans, 24, 9, 1.0, Km, 9, ud_imag, 1, 0.0, Ua_imag, 1);

        // now interpolate
        ug_real[0] = C8_interpolate(gp,ux_real) + (1.0-gp2.x)*Ua_real[0] + (1.0-gp2.y)*Ua_real[1] + (1.0-gp2.z)*Ua_real[2]; // x disp at Gauss point
        ug_imag[0] = C8_interpolate(gp,ux_imag) + (1.0-gp2.x)*Ua_imag[0] + (1.0-gp2.y)*Ua_imag[1] + (1.0-gp2.z)*Ua_imag[2]; // x disp at Gauss point
        ug_real[1] = C8_interpolate(gp,uy_real) + (1.0-gp2.x)*Ua_real[3] + (1.0-gp2.y)*Ua_real[4] + (1.0-gp2.z)*Ua_real[5]; // y disp at Gauss point
        ug_imag[1] = C8_interpolate(gp,uy_imag) + (1.0-gp2.x)*Ua_imag[3] + (1.0-gp2.y)*Ua_imag[4] + (1.0-gp2.z)*Ua_imag[5]; // y disp at Gauss point
        ug_real[2] = C8_interpolate(gp,uz_real) + (1.0-gp2.x)*Ua_real[6] + (1.0-gp2.y)*Ua_real[7] + (1.0-gp2.z)*Ua_real[8]; // z disp at Gauss point
        ug_imag[2] = C8_interpolate(gp,uz_imag) + (1.0-gp2.x)*Ua_imag[6] + (1.0-gp2.y)*Ua_imag[7] + (1.0-gp2.z)*Ua_imag[8]; // z disp at Gauss point

        // compute mass complex eigenvalue sentivity
        ctemp.r = dg_real[0]*ug_real[0] - dg_imag[0]*ug_imag[0];
        ctemp.r += dg_real[1]*ug_real[1] - dg_imag[1]*ug_imag[1];
        ctemp.r += dg_real[2]*ug_real[2] - dg_imag[2]*ug_imag[2];

        ctemp.i = dg_real[0]*ug_imag[0] + dg_imag[0]*ug_real[0];
        ctemp.i += dg_real[1]*ug_imag[1] + dg_imag[1]*ug_real[1];
        ctemp.i += dg_real[2]*ug_imag[2] + dg_imag[2]*ug_real[2];

        ftemp = alpha*rho;
        ctemp.r *= ftemp; ctemp.i *= ftemp;

        // multiply mass sensitivity by eigenvalue
        mtemp.r = ctemp.r * eig[0].r - ctemp.i * eig[0].i;
        mtemp.i = ctemp.r * eig[0].i + ctemp.i * eig[0].r;

        // finally add stiffness and mass sensitivities and multiply by the weight
        stemp.r += mtemp.r; stemp.i += mtemp.i;

        ctemp.r = (stemp.r*wgt_case[1].r - stemp.i*wgt_case[1].i);
        ctemp.i = (stemp.r*wgt_case[1].i + stemp.i*wgt_case[1].r);

        sens = ctemp.r * wgt_case[0].r + ctemp.i * wgt_case[0].i;

        // add to sensitivty array
        ind=pNum+p;
        gSens[ind] = -sens; // correct sign for maximization (i.e -ve of real gradient)
    }
}

// compute C8M sensitivities at all gauss points, for complex eigenvalue problems
void C8M::shpSens_CEig2(int numDual, int numEig, dcompLPK **eig_flut, double **eig_vec,
                       dcompLPK *eig, dcompLPK *wgt_case, double *gSens, int pNum, double alpha)
{
    // numEig is the number of modes used in the reduced model
    // eig_vec points to each globaleigenvector in turn
    // eig_flut points to the 4 eigenvectors at the flutter point
        // [0] = right_1 , [1] = right_2
        // [2] = left_1 , [3] = left_2

    int i,j,m1,m2,p,d,ind;
    double vec1[24], vec2[24]; // eigenvectors at element nodes


    double ftemp;
    Coord3D np[8], gp, gp2;
	double B[144];	// Strain displacement matrix
    double Ba[54];  // Strain displacement matrix for internal dof
    dcompLPK stemp, ctemp, mtemp;

    double stn[6];  // strain tensor
    double strs[6]; // stress tensor
    double stn_a[6]; // strain tensor for internal dof
    double Ua[9]; // displacements for internal dof
    double sens; // variable to keep running total of a sensitivity value
    double ux[8], uy[8], uz[8];
    double px[8], py[8], pz[8];
    double ug[3];

    // shape sens matrices
    i = numEig*numEig;
    dcompLPK *dK = new dcompLPK [i]();
    dcompLPK *dM = new dcompLPK [i]();
    dcompLPK *aux = new dcompLPK [numEig];
    dcompLPK one; one.r = 1.0; one.i = 0.0;
    dcompLPK zero; zero.r = 0.0; zero.i = 0.0;

	// get element nodal coords
	getNcrd(np);

    // get material property matrix
    double *Em = Mat->getD();
    double rho = Mat->getDens();
    double Malpha = (alpha < 1.01E-3) ? 1.0E-5 : alpha; // extra penalization for void mass

    // for each gauss point
    for(p=0;p<8;p++)
    {
        gp.x = gauss2 * c8_pos[p].x;
		gp.y = gauss2 * c8_pos[p].y;
		gp.z = gauss2 * c8_pos[p].z; // non-dim gauss point coords
        gp2.x = gp.x*gp.x; gp2.y = gp.y*gp.y; gp2.z = gp.z*gp.z;

        // compute isoparametric B matrix (6 x 24) & Ba (6 x 9)
        C8_isoBMat(gp, np, B); // NB: gp should be in non-dim coordinates (i.e center = 0,0,0)
        C8M_isoBMat(gp, np, Ba);
        stemp.r = 0.0; stemp.i = 0.0;
        mtemp.r = 0.0; mtemp.i = 0.0; // re-set

        // loop over each combination of modes
        for(m1=0;m1<numEig;m1++)
        {
            getDisp(vec1, eig_vec[m1]);

            for(i=0;i<8;i++)
            {
                j=i*3;
                ux[i] = vec1[j];
                uy[i] = vec1[j+1];
                uz[i] = vec1[j+2];
            }

            // compute stress tensor Ee(u) = strs
            // multiply: strain = B x u
            cblas_dgemv(CblasRowMajor, CblasNoTrans, 6, 24, 1.0, B, 24, vec1, 1, 0.0, stn, 1);
            // multiply: Ua = Km x u
            cblas_dgemv(CblasRowMajor, CblasTrans, 24, 9, 1.0, Km, 9, vec1, 1, 0.0, Ua, 1);
            // multiply: stn_a = Ba x Ua
            cblas_dgemv(CblasRowMajor, CblasNoTrans, 6, 9, 1.0, Ba, 9, Ua, 1, 0.0, stn_a, 1);
            // add the two strain tensors;
            for(i=0;i<6;i++) { stn[i] += stn_a[i]; }
            // Finally multiply: stress = material matrix x strain
            cblas_dgemv(CblasRowMajor, CblasNoTrans, 6, 6, 1.0, Em, 6, stn, 1, 0.0, strs, 1);

            ug[0] = C8_interpolate(gp,ux) + (1.0-gp2.x)*Ua[0] + (1.0-gp2.y)*Ua[1] + (1.0-gp2.z)*Ua[2]; // x disp at Gauss point
            ug[1] = C8_interpolate(gp,uy) + (1.0-gp2.x)*Ua[3] + (1.0-gp2.y)*Ua[4] + (1.0-gp2.z)*Ua[5]; // y disp at Gauss point
            ug[2] = C8_interpolate(gp,uz) + (1.0-gp2.x)*Ua[6] + (1.0-gp2.y)*Ua[7] + (1.0-gp2.z)*Ua[8]; // z disp at Gauss point

            for(m2=m1;m2<numEig;m2++) // symmetry
            {
                getDisp(vec2, eig_vec[m2]);

                // -- Stiffness sensitvity -- //

                // compute dual strain tensor e(p) = stn
                // multiply: strain = B x u
                cblas_dgemv(CblasRowMajor, CblasNoTrans, 6, 24, 1.0, B, 24, vec2, 1, 0.0, stn, 1);
                // multiply: Ua = Km x u
                cblas_dgemv(CblasRowMajor, CblasTrans, 24, 9, 1.0, Km, 9, vec2, 1, 0.0, Ua, 1);
                // multiply: stn_a = Ba x Ua
                cblas_dgemv(CblasRowMajor, CblasNoTrans, 6, 9, 1.0, Ba, 9, Ua, 1, 0.0, stn_a, 1);
                // add the two strain tensors;
                for(i=0;i<6;i++) { stn[i] += stn_a[i]; }

                // compute stress(u) x strain(p)
                ftemp = 0.0; // reset to zero
                for(i=0;i<6;i++) { ftemp += strs[i]*stn[i]; }
                ftemp *= alpha ; // need to multiply by the volume ratio (account for modified modulus)

                dK[index2(m1,m2,numEig)].r = ftemp;
                if(m1!=m2){ dK[index2(m2,m1,numEig)].r = ftemp; } // symmetry

                // -- mass sensitvity -- //

                // interpolate, then multiply
                // need to extract u, v, w disp in seperate vectors
                for(i=0;i<8;i++)
                {
                    j=i*3;
                    px[i] = vec2[j];
                    py[i] = vec2[j+1];
                    pz[i] = vec2[j+2];
                }

                // now interpolate and multiply by vec1 values
                ftemp =  ug[0]*(C8_interpolate(gp,px) + (1.0-gp2.x)*Ua[0] + (1.0-gp2.y)*Ua[1] + (1.0-gp2.z)*Ua[2]);
                ftemp += ug[1]*(C8_interpolate(gp,py) + (1.0-gp2.x)*Ua[3] + (1.0-gp2.y)*Ua[4] + (1.0-gp2.z)*Ua[5]);
                ftemp += ug[2]*(C8_interpolate(gp,pz) + (1.0-gp2.x)*Ua[6] + (1.0-gp2.y)*Ua[7] + (1.0-gp2.z)*Ua[8]);
                // multiply by denstiy and alpha
                ftemp *= rho*Malpha;

                dM[index2(m1,m2,numEig)].r = ftemp;
                if(m1!=m2){ dM[index2(m2,m1,numEig)].r = ftemp; } // symmetry
            }
        }

        // now loop over each dual (or matched point 1 for flutter = obj, many if flutter = cnstr)
        for(d=0;d<numDual;d++)
        {
            // L2^H * dK * R1
            //cblas_zgemv(CblasRowMajor, CblasNoTrans, numEig, numEig, &one, dK, numEig, eig_flut[d*4], 1, &zero, aux, 1); - CEDIT
            //cblas_zdotc_sub(numEig, eig_flut[d*4 + 3], 1, aux, 1, &stemp); - CEDIT

            // L2^H * dM * R2
            //cblas_zgemv(CblasRowMajor, CblasNoTrans, numEig, numEig, &one, dM, numEig, eig_flut[d*4 + 1], 1, &zero, aux, 1); - CEDIT
            //cblas_zdotc_sub(numEig, eig_flut[d*4 + 3], 1, aux, 1, &mtemp); - CEDIT

            // multiply mass sensitivity by eigenvalue
            ctemp.r = mtemp.r * eig[d].r - mtemp.i * eig[d].i;
            ctemp.i = mtemp.r * eig[d].i + mtemp.i * eig[d].r;

            // add to stiffness sensitivity
            stemp.r += ctemp.r; stemp.i += ctemp.i;

            // multiply by inverse of denominator
            ind = 2*d + 1;
            ctemp.r = (stemp.r*wgt_case[ind].r - stemp.i*wgt_case[ind].i);
            ctemp.i = (stemp.r*wgt_case[ind].i + stemp.i*wgt_case[ind].r);

            // multiply by weights
            ind--;
            sens = ctemp.r * wgt_case[ind].r + ctemp.i * wgt_case[ind].i;

            // add to sensitivty array
            ind=pNum+p;
            gSens[index2(ind,d,numDual)] = -sens; // sensitivity for this matched point
        }
    }

    delete [] dK;
    delete [] dM;
    delete [] aux;
}

// compute C8M sensitivities at all gauss points
void C8M::shpSens(int numDual, int numCase, double **disp_prim, double **disp_dual,
                  double *wgt_fact, double *wgt_case, double *gSens, int pNum, double *Nf, double alpha, int mode)
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
    double stn_a[6]; // strain tensor for internal dof
    double Ua[9];  // displacements for internal dof
    double sens; // variable to keep running total of a sensitivity value
    double /*ux[8], uy[8],*/ uz[8];
    double /*px[8], py[8],*/ pz[8];
    double sw_fact;

	// get element nodal coords
	getNcrd(np);

    // get material property matrix
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

            // compute stress tensor Ee(u) = strs
            // multiply: strain = B x u
            cblas_dgemv(CblasRowMajor, CblasNoTrans, 6, 24, 1.0, B, 24, ud, 1, 0.0, stn, 1);
            // multiply: Ua = Km x u
            cblas_dgemv(CblasRowMajor, CblasTrans, 24, 9, 1.0, Km, 9, ud, 1, 0.0, Ua, 1);
            // multiply: stn_a = Ba x Ua
            cblas_dgemv(CblasRowMajor, CblasNoTrans, 6, 9, 1.0, Ba, 9, Ua, 1, 0.0, stn_a, 1);
            // add the two strain tensors;
            for(i=0;i<6;i++) { stn[i] += stn_a[i]; }
            // Finally multiply: stress = material matrix x strain
            cblas_dgemv(CblasRowMajor, CblasNoTrans, 6, 6, 1.0, Em, 6, stn, 1, 0.0, strs, 1);

            // compute self weight factor fg * u
            ftemp = Mat->getDens() * Nf[l] * -9.81 * alpha; // mass times acceleration (in +ve z direction)
            sw_fact = (mode==1) ? 0.0 : ftemp * (C8_interpolate(gp,uz) + (1.0-gp.x*gp.x)*Ua[6] + (1.0-gp.y*gp.y)*Ua[7] + (1.0-gp.z*gp.z)*Ua[8]);

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
                // multiply: Ua = Km x u
                cblas_dgemv(CblasRowMajor, CblasTrans, 24, 9, 1.0, Km, 9, dd, 1, 0.0, Ua, 1);
                // multiply: stn_a = Ba x Ua
                cblas_dgemv(CblasRowMajor, CblasNoTrans, 6, 9, 1.0, Ba, 9, Ua, 1, 0.0, stn_a, 1);
                // add the two strain tensors;
                for(i=0;i<6;i++) { stn[i] += stn_a[i]; }

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
                sens -= ftemp * (C8_interpolate(gp,pz) + (1.0-gp.x*gp.x)*Ua[6] + (1.0-gp.y*gp.y)*Ua[7] + (1.0-gp.z*gp.z)*Ua[8]);
                sens -= const_fact;
                sens *= wgt_case[l]; // multiply by load case weight

                // add to overall sensitivity
                ind=pNum+p;
                gSens[index2(ind,d,numDual)] += sens; // add value for this load case
            }
        }
    }

}

// compute C8M element density sensitivity
void C8M::elemSens(int numDual, int numCase, double **disp_prim, double **disp_dual,
                  double *wgt_fact, double *wgt_case, double *elSens, int eNum, double *Nf, int mode)
{
    int l,d,ind;
    double const_fact;
    double ud[24];  // Primary element displacement array
    double dd[24];  // Dual element displacement array
    double Ku[24];  // strain tensor
    double sens;    // variable to keep running total of a sensitivity value
    double K_temp[576];
    double ftemp;

    // expand triangle storage
    int i,j;
    for(i=0;i<24;i++) // row
    {
        for(j=0;j<24;j++) // col
        {
            K_temp[index2(i,j,24)] = (j < i) ? Ke[tri_ind(24,j,i)] : Ke[tri_ind(24,i,j)];
        }
    }

    // for each load case
    for(l=0;l<numCase;l++)
    {
        // compute the constant sensitivty factor for the element
        const_fact = wgt_fact[l] * Mat->getDens() * volume; // DO NOT multiply by the volume ratio

        // get the primary element displacements
        getDisp(ud, disp_prim[l]);

        // for each dual response p
        for(d=0;d<numDual;d++)
        {
            // add the 3 componnents (multipled by the load case weight)
            // -pKu + fg(u+p) + const_fact

            // get the dual element displacements
            getDisp(dd, disp_dual[index2(l,d,numDual)]);

            // compute pKu
            // multiply: K x u
            cblas_dgemv(CblasRowMajor, CblasNoTrans, 24, 24, 1.0, K_temp, 24, ud, 1, 0.0, Ku, 1);

            // multiply:p x Ku
            sens = -1.0 * cblas_ddot(24, dd, 1, Ku, 1);

            // multiply fg(u+p)
            ftemp = Mat->getDens() * volume * Nf[l] * -9.81; // mass times acceleration (in +ve z direction)
            ftemp *= 0.125; // 8 nodes
            if(mode!=1)
            {
                for(i=0;i<8;i++)
                {
                    ind = i*3;
                    sens += ftemp * ud[ind+2];
                }
            }

            for(i=0;i<8;i++)
            {
                ind = i*3;
                sens += ftemp * dd[ind+2];
            }

            // complete sensitivity computation
            sens += const_fact;
            sens *= wgt_case[l]; // multiply by load case weight

            // add to overall sensitivity
            elSens[index2(eNum,d,numDual)] += sens; // add value for this load case
        }
    }

}

// compute C8M element density sensitivity - complex Eigenvalue
void C8M::elemSens_CEig(int numDual, int numEig, dcompLPK **disp_prim, dcompLPK **disp_dual,
                        dcompLPK *eig, dcompLPK *wgt_case, double *elSens, double K_fact, int eNum)
{
    dcompLPK ud[24], dd[24], stemp, ctemp, mtemp;
    double ud_real[24]; double ud_imag[24];  // Primary element displacement array
    double dd_real[24]; double dd_imag[24];  // Dual element displacement array
    double Ku_real[24]; double Ku_imag[24]; // strain tensor
    double sens; // variable to keep running total of a sensitivity value
    double K_temp[576];
    double M_temp[576];

    // expand triangle storage
    int i,j;
    for(i=0;i<24;i++) // row
    {
        for(j=0;j<24;j++) // col
        {
            K_temp[index2(i,j,24)] = (j < i) ? Ke[tri_ind(24,j,i)] : Ke[tri_ind(24,i,j)];
            M_temp[index2(i,j,24)] = (j < i) ? Me[tri_ind(24,j,i)] : Me[tri_ind(24,i,j)];
        }
    }

    // get 1st right hand eigenvector
    getDispC(ud, disp_prim[0]);
    for(i=0;i<24;i++)
    {
        ud_real[i] = ud[i].r; ud_imag[i] = ud[i].i;
    }


    // get 2nd left hand eigenvector
    getDispC(dd, disp_dual[1]);
    for(i=0;i<24;i++)
    {
        dd_real[i] = dd[i].r; dd_imag[i] = -dd[i].i; // complex conjugate
    }

    // compute pKu
    // multiply: K x u
    cblas_dgemv(CblasRowMajor, CblasNoTrans, 24, 24, 1.0, K_temp, 24, ud_real, 1, 0.0, Ku_real, 1);
    cblas_dgemv(CblasRowMajor, CblasNoTrans, 24, 24, 1.0, K_temp, 24, ud_imag, 1, 0.0, Ku_imag, 1);

    // compute Real & Imaginary part of p x Ku - stiffness sensitivity
    stemp.r = 0.0; stemp.i = 0.0; // reset to zero
    for(i=0;i<24;i++)
    {
        stemp.r += Ku_real[i]*dd_real[i] - Ku_imag[i]*dd_imag[i];
        stemp.i += Ku_real[i]*dd_imag[i] + dd_real[i]*Ku_imag[i];
    }

    // get the 2nd right hand eigenvector
    getDispC(ud, disp_prim[1]);
    for(i=0;i<24;i++)
    {
        ud_real[i] = ud[i].r; ud_imag[i] = ud[i].i;
    }

    // compute pMu
    // multiply: M x u
    cblas_dgemv(CblasRowMajor, CblasNoTrans, 24, 24, 1.0, M_temp, 24, ud_real, 1, 0.0, Ku_real, 1);
    cblas_dgemv(CblasRowMajor, CblasNoTrans, 24, 24, 1.0, M_temp, 24, ud_imag, 1, 0.0, Ku_imag, 1);

    // compute Real & Imaginary part of p x Mu - mass sensitivity
    ctemp.r = 0.0; ctemp.i = 0.0; // reset to zero
    for(i=0;i<24;i++)
    {
        ctemp.r += Ku_real[i]*dd_real[i] - Ku_imag[i]*dd_imag[i];
        ctemp.i += Ku_real[i]*dd_imag[i] + dd_real[i]*Ku_imag[i];
    }

    // multiply mass sensitivity by eigenvalue * factor
    mtemp.r = ctemp.r * eig[0].r - ctemp.i * eig[0].i;
    mtemp.i = ctemp.r * eig[0].i + ctemp.i * eig[0].r;

    // finally add stiffness and mass sensitivities and multiply by the weight
    stemp.r += mtemp.r; stemp.i += mtemp.i;

    ctemp.r = (stemp.r*wgt_case[1].r - stemp.i*wgt_case[1].i);
    ctemp.i = (stemp.r*wgt_case[1].i + stemp.i*wgt_case[1].r);

    sens = ctemp.r * wgt_case[0].r + ctemp.i * wgt_case[0].i;

    // add to overall sensitivity
    elSens[eNum] = -sens; // add value for this load case
}

// compute C8M element density sensitivity - complex Eigenvalue
void C8M::elemSens_CEig2(int numDual, int numEig, double **eig_vecs, dcompLPK **right, dcompLPK **left,
                        dcompLPK eig, dcompLPK *wgt_case, double *elSens, double K_fact, int eNum)
{
    int l;
    dcompLPK stemp, ctemp, mtemp;
    double eig_loc[360];
    double *mult = new double [24*numEig];
    double *dK_dx = new double[numEig*numEig];
    double *dM_dx = new double[numEig*numEig];
    double Ku_real[24]; double Ku_imag[24];
    double ud_real[24]; double ud_imag[24];  // Primary element displacement array
    double dd_real[24]; double dd_imag[24];  // Dual element displacement array
    double sens; // variable to keep running total of a sensitivity value
    double K_temp[576];
    double M_temp[576];

    // expand triangle storage
    int i,j;
    for(i=0;i<24;i++) // row
    {
        for(j=0;j<24;j++) // col
        {
            K_temp[index2(i,j,24)] = (j < i) ? Ke[tri_ind(24,j,i)] : Ke[tri_ind(24,i,j)];
            M_temp[index2(i,j,24)] = (j < i) ? Me[tri_ind(24,j,i)] : Me[tri_ind(24,i,j)];
        }
    }

    // for each eigenvalue - make transformation matrix
    for(l=0;l<numEig;l++)
    {
        getDisp(&eig_loc[24*l], eig_vecs[l]);
    }

    // matrix multiplication
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, 24, numEig, 24, 1.0, K_temp, 24, eig_loc, 24, 0.0, mult, numEig);
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, numEig, numEig, 24, 1.0, eig_loc, 24, mult, numEig, 0.0, dK_dx, numEig);

    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, 24, numEig, 24, 1.0, M_temp, 24, eig_loc, 24, 0.0, mult, numEig);
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, numEig, numEig, 24, 1.0, eig_loc, 24, mult, numEig, 0.0, dM_dx, numEig);

    // get 1st right hand eigenvector
    for(i=0;i<numEig;i++)
    {
        ud_real[i] = right[0][i].r; ud_imag[i] = right[0][i].i;
    }

    // get 2nd left hand eigenvector
    for(i=0;i<numEig;i++)
    {
        dd_real[i] = left[1][i].r; dd_imag[i] = -left[1][i].i; // complex conjugate
    }

    // compute p(dK_dx)u
    // multiply: K x u
    cblas_dgemv(CblasRowMajor, CblasNoTrans, numEig, numEig, 1.0, dK_dx, numEig, ud_real, 1, 0.0, Ku_real, 1);
    cblas_dgemv(CblasRowMajor, CblasNoTrans, numEig, numEig, 1.0, dK_dx, numEig, ud_imag, 1, 0.0, Ku_imag, 1);

    // compute Real & Imaginary part of p x Ku - stiffness sensitivity
    stemp.r = 0.0; stemp.i = 0.0; // reset to zero
    for(i=0;i<numEig;i++)
    {
        stemp.r += Ku_real[i]*dd_real[i] - Ku_imag[i]*dd_imag[i];
        stemp.i += Ku_real[i]*dd_imag[i] + dd_real[i]*Ku_imag[i];
    }

    // get the 2nd right hand eigenvector
    for(i=0;i<24;i++)
    {
       ud_real[i] = right[1][i].r; ud_imag[i] = right[1][i].i;
    }

    // compute pMu
    // multiply: M x u
    cblas_dgemv(CblasRowMajor, CblasNoTrans, numEig, numEig, 1.0, dM_dx, numEig, ud_real, 1, 0.0, Ku_real, 1);
    cblas_dgemv(CblasRowMajor, CblasNoTrans, numEig, numEig, 1.0, dM_dx, numEig, ud_imag, 1, 0.0, Ku_imag, 1);

    // compute Real & Imaginary part of p x Mu - mass sensitivity
    ctemp.r = 0.0; ctemp.i = 0.0; // reset to zero
    for(i=0;i<numEig;i++)
    {
        ctemp.r += Ku_real[i]*dd_real[i] - Ku_imag[i]*dd_imag[i];
        ctemp.i += Ku_real[i]*dd_imag[i] + dd_real[i]*Ku_imag[i];
    }

    // multiply mass sensitivity by eigenvalue * factor
    mtemp.r = ctemp.r * eig.r - ctemp.i * eig.i;
    mtemp.i = ctemp.r * eig.i + ctemp.i * eig.r;

    // finally add stiffness and mass sensitivities and multiply by the weight
    stemp.r += mtemp.r; stemp.i += mtemp.i;

    ctemp.r = (stemp.r*wgt_case[1].r - stemp.i*wgt_case[1].i);
    ctemp.i = (stemp.r*wgt_case[1].i + stemp.i*wgt_case[1].r);

    sens = ctemp.r * wgt_case[0].r + ctemp.i * wgt_case[0].i;

    // add to overall sensitivity
    elSens[eNum] = -sens; // add value for this load case

    //delete [] eig_loc;
    delete [] mult;
    delete [] dK_dx;
    delete [] dM_dx;
}

