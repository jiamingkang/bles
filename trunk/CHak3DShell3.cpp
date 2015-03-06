/*
 * CHak3DShell3.cpp
 *
 *  Created on: 23 Feb 2015
 *      Author: knai20-admin
 */

#include "CHak3DShell3.h"

CHak3DShell_3::CHak3DShell_3() {
	// TODO Auto-generated constructor stub
	// constructor sets num nodes
	numNode=3;
	nodeSpace();
    t = -1.0;
    volume = -1.0; // indicates volume not computed
	Ke=0;
	Me=0;
    Kstress=0;
    tr_yes=false;

    // create a plate element (using local coord system)
    plate = new DKT;
    plate->setGNode(local);

    // create a membrane element
    mem = new Tri3;
    mem->setGNode(local);

    // set local node numbers
    int i;
    for(i=0;i<3;i++)
    {
        plate->addNode(i,i);
        mem->addNode(i,i);
    }

}

CHak3DShell_3::~CHak3DShell_3() {
	// TODO Auto-generated destructor stub
	// destructor (delete memory)
	delMat();
    delKs();
}

// compute transformation matrix and local coords
void CHak3DShell_3::makeTrans()
{
    // get global node coords
    Coord3D input[3]; getNcrd(input);

    // create local axes
    // x-axis vector from vertex 1 to vertex 2
    Coord3D x_local;
    x_local.x = input[1].x - input[0].x;
    x_local.y = input[1].y - input[0].y;
    x_local.z = input[1].z - input[0].z;
    x_local = vectorNorm(x_local);

    // z-axis normal to the plane
    Coord3D z_local = PlaneNorm3(input);
    z_local = vectorNorm(z_local);
    //z_local.x *= -1.0;  z_local.y *= -1.0;  z_local.z *= -1.0;

    // y-axis perpendicualr to x & z axes
    Coord3D y_local = CrossProd3(x_local, z_local);
    y_local = vectorNorm(y_local);

    // the direction cosine matrix is simply the values of the local unit vectors
    // convention is that u_glob = [Trans] x u_local
    Trans[0] = x_local.x; Trans[1] = y_local.x; Trans[2] = z_local.x;
    Trans[3] = x_local.y; Trans[4] = y_local.y; Trans[5] = z_local.y;
    Trans[6] = x_local.z; Trans[7] = y_local.z; Trans[8] = z_local.z;

    // use the resulting tranformation matrix to make all z coords of the input the same
    // i.e. the triangle will lie on a x,y plane in the global coordiante system
    local[0].x = x_local.x * input[0].x + x_local.y * input[0].y + x_local.z * input[0].z;
    local[0].y = y_local.x * input[0].x + y_local.y * input[0].y + y_local.z * input[0].z;

    local[1].x = x_local.x * input[1].x + x_local.y * input[1].y + x_local.z * input[1].z;
    local[1].y = y_local.x * input[1].x + y_local.y * input[1].y + y_local.z * input[1].z;

    local[2].x = x_local.x * input[2].x + x_local.y * input[2].y + x_local.z * input[2].z;
    local[2].y = y_local.x * input[2].x + y_local.y * input[2].y + y_local.z * input[2].z;

    local[0].z = 0.0; local[1].z = 0.0; local[2].z = 0.0; // flat in x-y plane

    tr_yes = true;
}

void CHak3DShell_3::makeVol()
{
    // geometric factors
    double y31 = local[2].y - local[0].y;
    double y12 = local[0].y - local[1].y;
    double x13 = local[0].x - local[2].x;
    double x21 = local[1].x - local[0].x;
    double a = x21*y31 - x13*y12; // 2 x area

    // compute volume
    volume = 0.5*a*t; // area x thickness
}

// compute the stiffness matrix (from DKT & Tri3 elements)
void CHak3DShell_3::makeKe()
{
    int i,j,k,ind;

    if(Ke){delete [] Ke;}
    Ke = new double [171](); // upper triangle of shell matrix

    if(t < 1.0e-12){ std::cout << "\n Error! thickness = " << t; }

    if(!tr_yes){makeTrans();} // make transfer matrix if required

    // make sub-matrices
    plate->setMat(Mat); mem->setMat(Mat); // ensure same material is used
    plate->makeKe();
    mem->makeKe();

    // transform them to gloabl coords (separately)
    double *Km = new double [324]();
    double *Kb = new double [324]();
    double *Tg = new double [324]();

    // first make the complete transformation matrix
    shell3_expandTrans(Tg, Trans);

    // expand local element matrices
    double ftemp, ftemp2;
    double *Kmtemp = mem->getKe(); double *Kbtemp = plate->getKe();

    //int rowb, rowm, colb, colm;
    for(i=0;i<9;i++) // row in local matrix
    {
        for(j=i;j<9;j++) // col
        {
            ind = tri_ind(9,i,j); // indicator to entry in local matrix
            ftemp = Kmtemp[ind];
            ftemp2 = Kbtemp[ind]; // local matrix entries

            Km[index2(dofmap_mem[i],dofmap_mem[j],18)] = ftemp;
            Kb[index2(dofmap_bnd[i],dofmap_bnd[j],18)] = ftemp2; // copy accross
            if(i != j){ Km[index2(dofmap_mem[j],dofmap_mem[i],18)] = ftemp;
                        Kb[index2(dofmap_bnd[j],dofmap_bnd[i],18)] = ftemp2; } // symetric matrix
        }
    }

    // reduce separate matrices to upper triangle only
    double mtemp[324]; // temp matrix used in multiplication

    // apply transformation (separately)
    // plate matrix
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, 18, 18, 18, 1.0, Kb, 18, Tg, 18, 0.0, mtemp, 18);
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 18, 18, 18, 1.0, Tg, 18, mtemp, 18, 0.0, Kb, 18);
    // membrane matrix
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, 18, 18, 18, 1.0, Km, 18, Tg, 18, 0.0, mtemp, 18);
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 18, 18, 18, 1.0, Tg, 18, mtemp, 18, 0.0, Km, 18);

    // reduce separate matrices to upper triangle only
    for(i=0;i<18;i++) // row
    {
        for(j=i;j<18;j++) // col
        {
            ind = tri_ind(18,i,j); // upper triangle indicator
            k = index2(i,j,18); // full matrix indicator
            Kbend[ind] = Kb[k];
            Kmem[ind] = Km[k];
            Ke[ind] = Kb[k] + Km[k]; // also make full matrix
        }
    }

    delete [] Km;
    delete [] Kb;
    delete [] Tg;

    // also delete original local matrices
    plate->delMat();
    mem->delMat();
}

// function to update data due to a thickness change
void CHak3DShell_3::change_thk(double tin)
{
    int i;
    double fact = tin / t; // new / old
    double fact3 = fact*fact*fact; // cubed

    t = tin; // update thickness
    volume *= fact; // update volume

    // update stiffness matrix (bending & membrane scale differently)
    if(Ke)
    {
        for(i=0;i<171;i++)
        {
            Kbend[i] *= fact3;
            Kmem[i] *= fact;
            Ke[i] = Kbend[i] + Kmem[i];
        }
    }

    // update mass matrix (scales linearly with thickness)
    if(Me)
    {
        for(i=0;i<171;i++) { Me[i] *= fact; }
    }
}

void CHak3DShell_3::makeMe()
{
    int i,j,k,ind;

    if(Me){delete [] Me;}
    Me = new double [171](); // upper triangle of shell matrix

    if(!tr_yes){makeTrans();} // make transfer matrix if required

    // make sub-matrices
    plate->setMat(Mat); mem->setMat(Mat); // ensure same material is used
    plate->makeMe();
    mem->makeMe();

    // transform them to gloabl coords
    double *Me_temp = new double [324]();
    double *Tg = new double [324]();

    // first make the complete transformation matrix
    shell3_expandTrans(Tg, Trans);

    // expand local element matrices
    double ftemp, ftemp2;
    double *Mmtemp = mem->getMe(); double *Mbtemp = plate->getMe();

    for(i=0;i<9;i++) // row in local matrix
    {
        for(j=i;j<9;j++) // col
        {
            ind = tri_ind(9,i,j); // indicator to entry in local matrix
            ftemp = Mmtemp[ind];
            ftemp2 = Mbtemp[ind]; // local matrix entries
            Me_temp[index2(dofmap_mem[i],dofmap_mem[j],18)] = ftemp;
            Me_temp[index2(dofmap_bnd[i],dofmap_bnd[j],18)] = ftemp2; // copy accross
            if(i != j){ Me_temp[index2(dofmap_mem[j],dofmap_mem[i],18)] = ftemp;
                Me_temp[index2(dofmap_bnd[j],dofmap_bnd[i],18)] = ftemp2; } // symetric matrix
        }
    }

    // apply transformation
    double mtemp[324]; // temp matrix used in multiplication
    // whole matrix
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, 18, 18, 18, 1.0, Me_temp, 18, Tg, 18, 0.0, mtemp, 18);
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 18, 18, 18, 1.0, Tg, 18, mtemp, 18, 0.0, Me_temp, 18);

    // reduce matrix to upper triangle only (consistent mass matrix)
    for(i=0;i<18;i++) // row
    {
        for(j=i;j<18;j++) // col
        {
            ind = tri_ind(18,i,j); // upper triangle indicator
            k = index2(i,j,18); // full matrix indicator
            Me[ind] = Me_temp[k];
        }
    }

    /*
    // reduce to diagonal matrix using HRZ lunping procedure
    double sumU=0.0; double sumV=0.0; double sumW=0.0;
    for(i=0;i<3;i++) // node number
    {
        j=6*i; // u dof
        sumU += Me[tri_ind(18,j,j)]; j++;
        sumV += Me[tri_ind(18,j,j)]; j++;
        sumW += Me[tri_ind(18,j,j)];
    }

    // scaling factors
    if(volume<0.0){makeVol();}
    double mass = volume * Mat->getDens();
    sumU = mass / sumU; sumV = mass / sumV; sumW = mass / sumW;

    for(i=0;i<3;i++) // node number
    {
        j=6*i;
        Me[tri_ind(18,j,j)] *= sumU; j++; // u dof
        Me[tri_ind(18,j,j)] *= sumV; j++; // v dof
        Me[tri_ind(18,j,j)] *= sumW; j++; // w dof
        Me[tri_ind(18,j,j)] *= sumV*sumW; j++; // rx dof
        Me[tri_ind(18,j,j)] *= sumU*sumW; j++; // ry dof
        Me[tri_ind(18,j,j)] *= sumU*sumV; j++; // rz dof
    }

    // zero out off-diagonal entries
    for(i=0;i<18;i++) // row
    {
        for(j=i+1;j<18;j++) // col
        {
            ind = tri_ind(18,i,j); // upper triangle indicator
            Me[ind] = 0.0;
        }
    }
    */

    delete [] Me_temp;
    delete [] Tg;

    // also delete original local matrices
    plate->delMat();
    mem->delMat();
}

// create stress stiffness matrix
void CHak3DShell_3::makeKs(double *disp)
{
    int i,j,ind,k,g;
    double strs[3]; // membrane stresses
    double Nmem[4]; // membrane forces
    double mtemp[324]; // aux matrix used in mutiplication
    double Ks_sub[9]; // 3x3 sub matrix for stress stiffness
    double *Ks_full = new double[324](); // 18x18 full order stress stiffness matrix

    Coord gauss[3];
	gauss[0].x = 0.5;
	gauss[0].y = 0.0;
	gauss[1].x = 0.5;
	gauss[1].y = 0.5;
	gauss[2].x = 0.0;
	gauss[2].y = 0.5; // coord for gauss intergration rule

    // compute isoparametic matrix for linear interpolation of {dw/dx , dw/dy}

    // geometric factors
    double factors[7];
    tri3_geomFact(local, factors);
    double a = factors[6]*6.0; // total weight factor for numerical integration

    double *G = factors;
    //G[0] = -factors[2]-factors[1]; G[1] = factors[1]; G[2] = factors[2];
    //G[3] = -factors[4]-factors[5]; G[4] = factors[4]; G[5] = factors[5];

    // get membrane displacement field - from shell displacement field (into local coords)
    double elDisp[18], memDisp[9];
    getDisp(elDisp, disp);

    if(!tr_yes){makeTrans();} // make transfer matrix if required
    double *Tg = new double [324]();
    shell3_expandTrans(Tg, Trans);

    // displacements in local coords
    double dispLoc[18]; // all dof for now
    cblas_dgemv(CblasRowMajor, CblasTrans, 18, 18, 1.0, Tg, 18, elDisp, 1, 0.0, dispLoc, 1);

    j=0;
    for(i=0;i<3;i++)
    {
        ind=6*i;
        memDisp[j++] = dispLoc[ind]; // x dof
        memDisp[j++] = dispLoc[ind+1]; // y dof
        memDisp[j++] = dispLoc[ind+5]; // wz dof
    }

    // for each gauss point
    for(g=0;g<3;g++)
    {
        // compute membrane forces
        // get stresses from membrane element
        mem->stress(gauss[g], memDisp, strs);

        // multiply by thickness to get forces / unit length
        Nmem[0] = strs[0]*t; Nmem[1] = strs[2]*t; Nmem[2] = Nmem[1]; Nmem[3] = strs[1]*t;

        // matrix multiplications
        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 2, 3, 2, 1.0, Nmem, 2, G, 3, 0.0, mtemp, 3);
        cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, 3, 3, 2, 1.0, G, 3, mtemp, 3, 0.0, Ks_sub, 3);

        // sum into expanded full order matrix (& apply factors)
        Ks_full[38] += Ks_sub[0]/a;  Ks_full[44] += Ks_sub[1]/a;  Ks_full[50] += Ks_sub[2]/a;
        Ks_full[146] += Ks_sub[3]/a; Ks_full[152] += Ks_sub[4]/a; Ks_full[158] += Ks_sub[5]/a;
        Ks_full[254] += Ks_sub[6]/a; Ks_full[260] += Ks_sub[7]/a; Ks_full[266] += Ks_sub[8]/a;
    }

    // transfer matrix
    double Ks_trans[324];

    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, 18, 18, 18, 1.0, Ks_full, 18, Tg, 18, 0.0, mtemp, 18);
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 18, 18, 18, 1.0, Tg, 18, mtemp, 18, 0.0, Ks_trans, 18);

    // reduce matrix to upper triangle only
    if(Kstress){delete [] Kstress;}
    Kstress = new double [171]();
    for(i=0;i<18;i++) // row
    {
        for(j=i;j<18;j++) // col
        {
            ind = tri_ind(18,i,j); // upper triangle indicator
            k = index2(i,j,18); // full matrix indicator
            Kstress[ind] = Ks_trans[k];
        }
    }

    delete [] Ks_full;
    delete [] Tg;
}

void CHak3DShell_3::getDisp(double *elDisp, double *globDisp)
{
	int i,j,nd;
	int ind = 0;

	for(i=0;i<3;i++) // for each node
	{
        if(dof_ind){ nd = dof_ind[getGlobal(i)]; }
        else{ nd = getGlobal(i)*6; } // 1st dof
		for(j=0;j<6;j++)
		{
			elDisp[ind++] = globDisp[nd++];
		}
	}
}

// compute Shell3 element thickness sensitivity
void CHak3DShell_3::elemSens(int numDual, int numCase, double **disp_prim, double **disp_dual,
                    double *wgt_fact, double *wgt_case, double *elSens, int eNum, double *Nf, int mode)
{
    // mode = 1: buckling
    // else, compliance

    int i,j,l,d,ind,ind2;
    double const_fact;
    double ud[18];  // Primary element displacement array
    double dd[18];  // Dual element displacement array
    double Ku[18];
    double sens; // variable to keep running total of a sensitivity value
    double dK_dt[324];
    double ftemp;

    // expand triangle storage
    // factor, to compute dK/dt
    double t3 = 3.0/t; // 3/t
    for(i=0;i<18;i++) // row
    {
        for(j=0;j<18;j++) // col
        {
            ind = index2(i,j,18);
            ind2 = (j < i) ? tri_ind(18,j,i) : tri_ind(18,i,j);
            dK_dt[ind] = (Kmem[ind2]/t) + (Kbend[ind2]*t3);
        }
    }

    // for each load case
    for(l=0;l<numCase;l++)
    {
        // compute the constant sensitivty factor for the element
        const_fact = (wgt_fact[l] * Mat->getDens() * volume) / t; // wgt_fact x dW/dt

        // get the primary element displacements
        getDisp(ud, disp_prim[l]);

        // compute p[dK/dt]u
        // multiply: dK/dt x u
        cblas_dgemv(CblasRowMajor, CblasNoTrans, 18, 18, 1.0, dK_dt, 18, ud, 1, 0.0, Ku, 1);

        // for each dual response p
        for(d=0;d<numDual;d++)
        {
            // add the 3 componnents (multipled by the load case weight)
            // -p[dK/dt]u + dfg/dt(u+p) + const_fact

            // get the dual element displacements
            getDisp(dd, disp_dual[index2(l,d,numDual)]);

            // multiply:p x Ku
            sens = -1.0 * cblas_ddot(18, dd, 1, Ku, 1);

            // multiply N x dfg/dt x (u+p)
            ftemp = Mat->getDens() * (volume/t) * Nf[l] * -9.81; // mass times acceleration (in +ve z direction)
            ftemp *= 0.333333333; // 3 nodes
            for(i=0;i<3;i++)
            {
                ind = i*6; // 6 dof / node
                ind += 2; // 3rd dof (z - direction)
                if(mode==1){sens += ftemp * dd[ind];}
                else{sens += ftemp * (ud[ind] + dd[ind]);}
            }

            // complete sensitivity computation
            sens += const_fact;
            sens *= wgt_case[l]; // multiply by load case weight

            // add to overall sensitivity
            elSens[index2(eNum,d,numDual)] += sens; // add value for this load case
        }
    }
}

// compute Shell3 element thickness sensitivity w.r.t an eigenvalue
void CHak3DShell_3::EigSens(int numDual, int numEig, double **disp_prim, double **disp_dual,
                     double *eig, double *wgt_case, double *elSens, int eNum, int mode)
{
    // mode = 1 -> use stress stiffness matrix (if available)
    // else -> use mass matrix

    int i,j,l,d,ind,ind2;
    double ud[18];  // Primary element displacement array
    double dd[18];  // Dual element displacement array
    double Ku[18], Mu[18];
    double sens; // variable to keep running total of a sensitivity value
    double dK_dt[324];
    double dM_dt[324];

    // expand triangle storage
    // factor, to compute dK/dt
    double t3 = 3.0/t; // 3/t
    for(i=0;i<18;i++) // row
    {
        for(j=0;j<18;j++) // col
        {
            ind = index2(i,j,18);
            ind2 = (j < i) ? tri_ind(18,j,i) : tri_ind(18,i,j);
            dK_dt[ind] = Kmem[ind2]/t + Kbend[ind2]*t3;
            if(mode==1 && Kstress){dM_dt[ind] = Kstress[ind2];}
            else if(mode!=1){dM_dt[ind] = Me[ind2];} // bring 1/t in later
        }
    }

    // for each eigenvalue
    for(l=0;l<numEig;l++)
    {
        // get the primary element displacements
        getDisp(ud, disp_prim[l]);

        // compute p[dK/dt]u & p[dM/dt]u
        // multiply: dK/dt x u & dM/dt x u
        cblas_dgemv(CblasRowMajor, CblasNoTrans, 18, 18, 1.0, dK_dt, 18, ud, 1, 0.0, Ku, 1);
        if(mode!=1 || Kstress){cblas_dgemv(CblasRowMajor, CblasNoTrans, 18, 18, 1.0, dM_dt, 18, ud, 1, 0.0, Mu, 1);}

        // for each dual response p
        for(d=0;d<numDual;d++)
        {
            // get the dual element displacements
            getDisp(dd, disp_dual[index2(l,d,numDual)]);

            // compute:p x Ku - eig x p x
            sens = cblas_ddot(18, dd, 1, Ku, 1);
            if(mode!=1 || Kstress){sens -= (eig[l]/t) * cblas_ddot(18, dd, 1, Mu, 1);}
            sens *= wgt_case[l]; // multiply by load case weight

            // add to overall sensitivity
            //elSens[index2(eNum,l,numEig)] += sens; // add value for this load case
            elSens[eNum] += sens;
        }
    }
}

// compute Shell3 flutter sensitivity - complex Eigenvalue
void CHak3DShell_3::CEigSens(int numDual, int numEig, dcompLPK **eig_flut, double **eig_vec,
                           dcompLPK *eig, dcompLPK *wgt_case, double *elSens, int eNum)
{
    int i,j,l,d,ind,ind2;
    dcompLPK stemp, ctemp, mtemp;
    double *eig_loc = new double [18*numEig];
    double sens; // variable to keep running total of a sensitivity value
    double K_temp[324];
    double M_temp[324];

    // sens matrices
    i = numEig*numEig;
    dcompLPK *dK_dt = new dcompLPK [i];
    dcompLPK *dM_dt = new dcompLPK [i];
    double *dK_temp = new double [i]();
    double *dM_temp = new double [i]();
    double *mult = new double [18*numEig];
    dcompLPK *aux = new dcompLPK [numEig];
    dcompLPK one; one.r = 1.0; one.i = 0.0;
    dcompLPK zero; zero.r = 0.0; zero.i = 0.0;

    // expand triangle storage
    // factor, to compute dK/dt
    double t3 = 3.0/t; // 3/t
    for(i=0;i<18;i++) // row
    {
        for(j=0;j<18;j++) // col
        {
            ind = index2(i,j,18);
            ind2 = (j < i) ? tri_ind(18,j,i) : tri_ind(18,i,j);
            K_temp[ind] = Kmem[ind2]/t + Kbend[ind2]*t3;
            M_temp[ind] = Me[ind2]; // bring 1/t in later
        }
    }

    // for each eigenvalue - make ROM transformation matrix
    for(l=0;l<numEig;l++)
    {
        getDisp(&eig_loc[18*l], eig_vec[l]);
    }

    // matrix multiplication
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, 18, numEig, 18, 1.0, K_temp, 18, eig_loc, 18, 0.0, mult, numEig);
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, numEig, numEig, 18, 1.0, eig_loc, 18, mult, numEig, 0.0, dK_temp, numEig);

    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, 18, numEig, 18, 1.0, M_temp, 18, eig_loc, 18, 0.0, mult, numEig);
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, numEig, numEig, 18, 1.0, eig_loc, 18, mult, numEig, 0.0, dM_temp, numEig);

    // make dK & dM complex for next section
    ind = numEig*numEig;
    for(i=0;i<ind;i++)
    {
        dK_dt[i].r = dK_temp[i]; dK_dt[i].i = 0.0;
        dM_dt[i].r = dM_temp[i]; dM_dt[i].i = 0.0;
    }

    delete [] dK_temp;
    delete [] dM_temp;
    delete [] eig_loc;
    delete [] mult; // redundant

    // now loop over each dual (or matched point 1 for flutter = obj, many if flutter = cnstr)
    for(d=0;d<numDual;d++)
    {
        // L2^H * dK * R1
        //cblas_zgemv(CblasRowMajor, CblasNoTrans, numEig, numEig, &one, dK_dt, numEig, eig_flut[d*4], 1, &zero, aux, 1); -CEDIT
        //cblas_zdotc_sub(numEig, eig_flut[d*4 + 3], 1, aux, 1, &stemp); -CEDIT

        // L2^H * dM * R2
        //cblas_zgemv(CblasRowMajor, CblasNoTrans, numEig, numEig, &one, dM_dt, numEig, eig_flut[d*4 + 1], 1, &zero, aux, 1); -CEDIT
        //cblas_zdotc_sub(numEig, eig_flut[d*4 + 3], 1, aux, 1, &mtemp); -CEDIT

        // multiply mass sensitivity by eigenvalue
        ctemp.r = mtemp.r * eig[d].r - mtemp.i * eig[d].i;
        ctemp.i = mtemp.r * eig[d].i + mtemp.i * eig[d].r;

        // add to stiffness sensitivity
        stemp.r += ctemp.r / t; stemp.i += ctemp.i / t; // bring back 1/t factor for mass sensitivity

        // multiply by inverse of denominator
        ind = 2*d + 1;
        ctemp.r = (stemp.r*wgt_case[ind].r - stemp.i*wgt_case[ind].i);
        ctemp.i = (stemp.r*wgt_case[ind].i + stemp.i*wgt_case[ind].r);

        // multiply by weights
        ind--;
        sens = ctemp.r * wgt_case[ind].r + ctemp.i * wgt_case[ind].i;

        // add to sensitivty array
        elSens[index2(eNum,d,numDual)] -= sens; // sensitivity for this matched point
    }

    delete [] aux;
    delete [] dK_dt;
    delete [] dM_dt;
}

// function to compute the contibution of an element to dks_du (buckling)
void CHak3DShell_3::DispSens(int numModes, double **disp_prim, double *fact, double *dk_du)
{
    int i,j,m,g,nd,ind,ind2;

    // get material property matrix
    double *Em = Mat->getDsts();

    Coord gauss[3];
	gauss[0].x = 0.5;
	gauss[0].y = 0.0;
	gauss[1].x = 0.5;
	gauss[1].y = 0.5;
	gauss[2].x = 0.0;
	gauss[2].y = 0.5; // coord for gauss intergration rule

    // geometric factors
    // geometric factors
    double factors[7];
    tri3_geomFact(local, factors);
    double wgt = t/(factors[6]*factors[6]*6.0); // total weight factor for numerical integration (including factor for B matrix & thk)

    // compute gx & gy (constant over the element)
    double *gx, *gy, ax, ay, axy;
    gx=factors; gy=&factors[3];

    if(!tr_yes){makeTrans();} // make transfer matrix if required
    double Tg[324];
    for(i=0;i<324;i++){Tg[i]=0.0;}
    shell3_expandTrans(Tg, Trans);

    // data arrays
    double elDisp[18], dispLoc[18], w_disp[3], B[27], B18[54], EB[54], dN_du[54], du_gauss[18];

    // for each mode
    for(m=0;m<numModes;m++)
    {
        // get displacements in local elemenet coord system
        getDisp(elDisp, disp_prim[m]);
        // displacements in local coords
        cblas_dgemv(CblasRowMajor, CblasTrans, 18, 18, 1.0, Tg, 18, elDisp, 1, 0.0, dispLoc, 1);
        // reduce to just out of plane dof
        w_disp[0] = dispLoc[2]; w_disp[1] = dispLoc[8]; w_disp[2] = dispLoc[14];

        // compute g x u terms -> a values
        ax = cblas_ddot(3, w_disp, 1, gx, 1);
        ay = cblas_ddot(3, w_disp, 1, gy, 1);
        axy = 2.0 * ax * ay;
        ax *= ax; ay *= ay;

        // for each integration point
        for(g=0;g<3;g++)
        {
            // compute strain-displacement matrix for membrane element
            tri3_Bmat(gauss[g], factors, B);

            // expand B for all dof
            for(i=0;i<3;i++) // each row
            {
                ind=18*i; ind2=9*i; // start of row
                // node 1
                B18[ind] = B[ind2]; B18[ind+1] = B[ind2+1]; B18[ind+2] = 0.0;
                B18[ind+3] = 0.0;   B18[ind+4] = 0.0;       B18[ind+5] = B[ind2+2];
                // node 2
                B18[ind+6] = B[ind2+3]; B18[ind+7] = B[ind2+4]; B18[ind+8] = 0.0;
                B18[ind+9] = 0.0;       B18[ind+10] = 0.0;      B18[ind+11] = B[ind2+5];
                // node 3
                B18[ind+12] = B[ind2+6]; B18[ind+13] = B[ind2+7]; B18[ind+14] = 0.0;
                B18[ind+15] = 0.0;       B18[ind+16] = 0.0;       B18[ind+17] = B[ind2+8];
            }

            // multiply E x B x T
            cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 3, 18, 3, 1.0, Em, 3, B18, 18, 0.0, EB, 18);
            cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, 3, 18, 18, 1.0, EB, 18, Tg, 18, 0.0, dN_du, 18);

            // compute contribution of dk_du from this gauss point
            for(i=0;i<18;i++)
            {
                du_gauss[i] = dN_du[i]*ax + dN_du[i+18]*ay + dN_du[i+36]*axy;
                du_gauss[i] *= wgt*fact[m]; // fact = eig / u^T Ks u
            }

            // sum dN/du x a values into gloabl vector dk_du
            for(i=0;i<3;i++)
            {
                ind = 6*i; // 1st dof in du_gauss
                nd = getGlobal(i); // global node number
                if(!dof_ind){nd*=6;}
                else{nd = dof_ind[nd];} // 1st degree of freedom (x)
                // for all dof
                for(j=0;j<6;j++)
                {
                    dk_du[nd++] += du_gauss[ind++];
                }
            }
        }
    }
}

// function to compute the contibution of an element to dvm_du (stress)
void CHak3DShell_3::DispSensVM(int numModes, double **disp_prim, double *fact, double *dvm_du)
{
    int i,j,m,nd,ind,ind2;
    // compute matrix ExBxT at element centre
    double *Em = Mat->getDsts(); // get material property matrix

    Coord cent; cent.x = 0.333333333333; cent.y = 0.333333333333;
    double B[27];
    double factors[7];
    tri3_geomFact(local, factors);
    tri3_Bmat(cent, factors, B);

    double EB[27], EBT[54], B18[54], elDisp[18], du[18], scale;
    scale = 1.0/(factors[6]*factors[6]); // scaling on B matrix (squared)

    if(!tr_yes){makeTrans();} // make transfer matrix if required
    double Tg[324];
    for(i=0;i<324;i++){Tg[i]=0.0;}
    shell3_expandTrans(Tg, Trans);

    // expand B for all dof
    for(i=0;i<3;i++) // each row
    {
        ind=18*i; ind2=9*i; // start of row
        // node 1
        B18[ind] = B[ind2]; B18[ind+1] = B[ind2+1]; B18[ind+2] = 0.0;
        B18[ind+3] = 0.0;   B18[ind+4] = 0.0;       B18[ind+5] = B[ind2+2];
        // node 2
        B18[ind+6] = B[ind2+3]; B18[ind+7] = B[ind2+4]; B18[ind+8] = 0.0;
        B18[ind+9] = 0.0;       B18[ind+10] = 0.0;      B18[ind+11] = B[ind2+5];
        // node 3
        B18[ind+12] = B[ind2+6]; B18[ind+13] = B[ind2+7]; B18[ind+14] = 0.0;
        B18[ind+15] = 0.0;       B18[ind+16] = 0.0;       B18[ind+17] = B[ind2+8];
    }

    // multiply E x B x T
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 3, 18, 3, 1.0, Em, 3, B18, 18, 0.0, EB, 18);
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, 3, 18, 18, 1.0, EB, 18, Tg, 18, 0.0, EBT, 18);

    // multiply (EBT)^T V EBT
    double aux[54];
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 3, 18, 3, scale, VM_matrix, 3, EBT, 18, 0.0, aux, 18);
    double final[324];
    cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, 18, 18, 3, 1.0, EBT, 18, aux, 18, 0.0, final, 18);

    // sum for each mode
    for(m=0;m<numModes;m++)
    {
        // get global displacements
        getDisp(elDisp, disp_prim[m]);

        // multiply with final
        cblas_dgemv(CblasRowMajor, CblasNoTrans, 18, 18, fact[m], final, 18, elDisp, 1, 0.0, du, 1);

        // add to dvm_du
        for(i=0;i<3;i++)
        {
            ind = 6*i; // 1st dof in du
            nd = getGlobal(i); // global node number
            if(!dof_ind){nd*=6;}
            else{nd = dof_ind[nd];} // 1st degree of freedom (x)
            // for all dof
            for(j=0;j<6;j++)
            {
                dvm_du[nd++] += du[ind++];
            }
        }
    }
}

// compute Shell3 strain tensor at a (non-dim) point in the element
void CHak3DShell_3::strain(Coord3D p, double *elDisp, double *stn)
{
    // p must be non-dimensional

    // transform elDisp into local membrane coords
    // make the complete transformation matrix
    double *Tg = new double [324]();
    shell3_expandTrans(Tg, Trans);

    // displacements in local coords
    double dispLoc[18]; // all dof for now
    cblas_dgemv(CblasRowMajor, CblasTrans, 18, 18, 1.0, Tg, 18, elDisp, 1, 0.0, dispLoc, 1);

    double dispLST[9]; // dof for membrane part
    dispLST[0] = dispLoc[0]; dispLST[1] = dispLoc[1]; dispLST[2] = dispLoc[5];
    dispLST[3] = dispLoc[6]; dispLST[4] = dispLoc[7]; dispLST[5] = dispLoc[11];
    dispLST[6] = dispLoc[12]; dispLST[7] = dispLoc[13]; dispLST[8] = dispLoc[17];

    // reduce p (non-dim coord) to 2D coord
    Coord pin; pin.x=p.x; pin.y=p.y;

	// assume in-plane strains from bending are small (ok for wing skins)
    mem->strain(pin, dispLST, stn);

    delete [] Tg;
    // leave in local coords (for now!)
}

// compute Shell3 stress tensor at a point in the element
void CHak3DShell_3::stress(Coord3D p, double *elDisp, double *strs)
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

// compute von Mises stress (at an input point)
double CHak3DShell_3::vonMises(Coord3D p, double *disp)
{
    // get displaement vector (global coords)
    double elDisp[18];
    getDisp(elDisp, disp);

    // get stresses
    double strs[3];
    stress(p,elDisp,strs);

    // compute & return von Mises
    double vm = strs[0]*strs[0] + strs[1]*strs[1] - strs[0]*strs[1] + 3.0*strs[2]*strs[2];
    return sqrt(vm);
}

// add self-weight loading to a global loading vector
void CHak3DShell_3::addWeight(Coord3D vec, double mag, int *dof_ind, double *gload)
{
    // first compute total element force (in each component)
    double force = Mat->getDens() * volume * mag; // force = mass x accel

    // split equally between each node
    force *= 0.333333333333;
    Coord3D eload;
    eload = vec * force; // force x vec

    // add force x vec for each node into gloabl load vector (gload)
    int i,nd;
    for(i=0;i<3;i++)
    {
        nd = getGlobal(i); // global node number
        if(!dof_ind){nd*=6;}
        else{nd = dof_ind[nd];} // 1st degree of freedom (x)
        gload[nd++] += eload.x;
        gload[nd++] += eload.y; // 2nd degree of freedom (y)
        gload[nd] += eload.z; // 3rd degree of freedom (z)
    }
}
