/*
 * CHakHomogenize.cpp
 *
 *  Created on: 22 Feb 2015
 *      Author: knai20-admin
 */

#include "CHakHomogenize.h"

CHakHomogenize::CHakHomogenize() {
	// TODO Auto-generated constructor stub

}

CHakHomogenize::~CHakHomogenize() {
	// TODO Auto-generated destructor stub
}


// function to homogenize the material properties of the microstructure - using symm periodic BCs
void CHakHomogenize::homogenize(mesh *inMesh, isoMat *mat1, sp_mat *K1, double *alpha, int **fixed, int freeDof, isoMat *homMat, double *disp)
{
    // read data
    int numDof = inMesh->NumNodes * NUM_DOF;
    int elemX = inMesh->elemX;
    int elemY = inMesh->elemY;
    Elem **Number = inMesh->Number;
    double area = (double)(elemX * elemY);

    int i,j,n,m,o,o2,num,temp,temp2;
    double atemp;
    int tdof[8];

    // --- compute load vector for each unit strain load case --- //

    // for now, assume all elements made from material 1
    double e11 = mat1->mat[0] * 0.5;
    double e22 = mat1->mat[4] * 0.5;
    double e12 = mat1->mat[1] * 0.5;
    double e66 = mat1->mat[8] * 0.5;
    double fb[24] = {-e11, -e12, e11, -e12, e11, e12, -e11, e12,
                     -e12, -e22, e12, -e22, e12, e22, -e12, e22,
                     -e66, -e66, -e66, e66, e66, e66, e66, -e66};

    double *hom_force = calloc(3*numDof, sizeof(double));

    // For All Elements
    for(m=0;m<elemY;m++)
    {
        for(n=0;n<elemX;n++)
        {
            num = Number[n][m].n; // Element number
            atemp = alpha[num]; // Element area ratio

            // Read in grid node dof of current element
            tdof[0] = NUM_DOF * Number[n][m].a;
            tdof[1] = tdof[0] + 1;
            tdof[2] = NUM_DOF * Number[n][m].b;
            tdof[3] = tdof[2] + 1;
            tdof[4] = NUM_DOF * Number[n][m].c;
            tdof[5] = tdof[4] + 1;
            tdof[6] = NUM_DOF * Number[n][m].d;
            tdof[7] = tdof[6] + 1;

            // for each load case, add loading contribution
            for(j=0;j<3;j++)
            {
                o = numDof*j;
                o2 = 8*j;
                for(i=0;i<8;i++)
                {
                    hom_force[tdof[i]+o] += fb[i+o2]*atemp;
                }
            }
        }
    }

    // --- solve each load case --- //

    // need to copy the stiffness matrix (2 types of BC)
    sp_mat K2; K2.irn = 0; K2.jcn = 0; K2.A = 0;
    copy_sp_mat(K1, &K2);

    // remove fixed dof
    rem_sp_mat(K1, fixed[0], 1);
    rem_sp_mat(&K2, fixed[1], 1);

    // remove fixed dof from hom_force
    double *hom_disp = malloc(3*numDof * sizeof(double));

    for(i=0;i<numDof;i++) // for all dof
    {
        // first set of BCs
        num = fixed[0][i];
        if(num > -1) // if dof not fixed
        {
            for(j=0;j<2;j++) // for first 2 load cases
            {
                n = numDof * j; // point to location in hom_force
                m = freeDof * j; // point to location in hom_disp
                hom_disp[m+num] = hom_force[n+i]; // copy load accross
            }
        }

        // second set of BCs
        num = fixed[1][i];
        if(num > -1) // if dof not fixed
        {
            n = numDof * 2; // point to location in hom_force
            m = numDof * 2; // point to location in hom_disp (ensure info correctly returned from FE_Solve)
            hom_disp[m+num] = hom_force[n+i]; // copy load accross
        }
    }

    // solve unit load cases
    FE_Solve(K1, fixed[0], hom_disp, freeDof, numDof, 2); // first 2
    FE_Solve(&K2, fixed[1], &hom_disp[2*numDof], freeDof, numDof, 1); // last 1

    //OutDispVTK(inMesh, 3, hom_disp, 0, 0, 0, "hom_disp");

    // --- compute the homogenized stiffness tensor --- //

    double stn[9];
    double Edisp[24];
    double Ehom[9];
    for(i=0;i<8;i++){homMat->mat[i] = 0.0;} // reset

    // For All Elements
    for(m=0;m<elemY;m++)
    {
        for(n=0;n<elemX;n++)
        {
            num = Number[n][m].n; // Element number
            atemp = alpha[num]; // Element area ratio

            if(atemp > 0.0)
            {
                // Read in grid node dof of current element
                tdof[0] = NUM_DOF * Number[n][m].a;
                tdof[1] = tdof[0] + 1;
                tdof[2] = NUM_DOF * Number[n][m].b;
                tdof[3] = tdof[2] + 1;
                tdof[4] = NUM_DOF * Number[n][m].c;
                tdof[5] = tdof[4] + 1;
                tdof[6] = NUM_DOF * Number[n][m].d;
                tdof[7] = tdof[6] + 1;

                // get element displacement arrays
                for(j=0;j<3;j++) // each case
                {
                    temp = j*8;
                    temp2 = j*numDof;
                    for(i=0;i<8;i++) // each dof
                    {
                        Edisp[temp + i] =  hom_disp[temp2 + tdof[i]];
                    }
                }

                // get element (centre) strains for each load case: stn = B x disp
                cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, 3, 3, 8, 0.5, BHom, 8, Edisp, 8, 0.0, stn, 3);
                // 1/2 factor for proper B matrix

                // multiply material property matrix with strains
                cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 3, 3, 3, 1.0, mat1->mat, 3, stn, 3, 0.0, Ehom, 3);

                // subtract from elem mat matrix & multiply by area ratio and sum into homogenized values
                // do not include some entires as shear test case BCs mess them up
                for(i=0;i<5;i++){j=Esymm_map[i]; homMat->mat[j] += atemp * (mat1->mat[j] - Ehom[j]);}
            }
        }
    }

    // also compute other homogenized properties
    for(i=0;i<9;i++){homMat->mat[i] /= area;}
    double E1, E2, v12, v21, G, B;
    E1 = homMat->mat[0] - ((homMat->mat[1] * homMat->mat[1]) / homMat->mat[4]);
    E2 = homMat->mat[4] - ((homMat->mat[1] * homMat->mat[1]) / homMat->mat[0]);
    v12 = homMat->mat[1] / homMat->mat[0];
    v21 = homMat->mat[1] / homMat->mat[4];
    G = homMat->mat[8];
    B = 0.5*(0.5*(homMat->mat[0] + homMat->mat[4]) + homMat->mat[1]);

    // save data
    homMat->g = G;
    homMat->k = B;
    homMat->e = E1;
    homMat->v = v12;

    printf("\nE[]:");
    for(i=0;i<9;i+=3)
    {
        printf("\n%12.4e\t%12.4e\t%12.4e",homMat->mat[i],homMat->mat[i+1],homMat->mat[i+2]);
    }

    printf("\nE1=%12.4e\nE2=%12.4e\nv12=%12.4e\nv21=%12.4e\nG=%12.4e\nB=%12.4e",E1, E2, v12, v21, G, B);

    // process displacement fields for sensitivity computation
    int NumNodes = inMesh->NumNodes;
    Coord * NodeCoord = inMesh->NodeCoord;

    // unit strain in x, y and shear
    temp = 2*numDof;
    for(i=0;i<NumNodes;i++)
    {
        j = i*2;
        disp[j] = NodeCoord[i].x;
        disp[j+1] = 0.0;
        disp[numDof + j] = 0.0;
        disp[numDof + j + 1] = NodeCoord[i].y;
        disp[temp + j] = NodeCoord[i].y;
        disp[temp + j + 1] = NodeCoord[i].x;
    }

    // difference in disp fields
    temp = 3 * numDof;
    for(i=0;i<temp;i++) { disp[i] -= hom_disp[i]; }
}

// function to homogenize the material properties of the microstructure - using constraint equn periodic BCs
void CHakHomogenize::homogenize2(mesh *inMesh, isoMat *mat1, sp_mat *K1, double *alpha, int **fixed, int freeDof, sp_mat *CE,
                 isoMat *homMat, double *disp, int pinfo)
{
    // read data
    int numDof = inMesh->NumNodes * NUM_DOF;
    int elemX = inMesh->elemX;
    int elemY = inMesh->elemY;
    Elem **Number = inMesh->Number;
    double area = (double)(elemX * elemY);

    int i,j,n,m,o,o2,num,temp,temp2;
    double atemp;
    int tdof[8];

    // --- compute load vector for each unit strain load case --- //

    // for now, assume all elements made from material 1
    double e11 = mat1->mat[0] * 0.5;
    double e22 = mat1->mat[4] * 0.5;
    double e12 = mat1->mat[1] * 0.5;
    double e66 = mat1->mat[8] * 0.5;
    double fb[24] = {-e11, -e12, e11, -e12, e11, e12, -e11, e12,
        -e12, -e22, e12, -e22, e12, e22, -e12, e22,
        -e66, -e66, -e66, e66, e66, e66, e66, -e66};

    double *hom_force = calloc(3*numDof, sizeof(double));

    // For All Elements
    for(m=0;m<elemY;m++)
    {
        for(n=0;n<elemX;n++)
        {
            num = Number[n][m].n; // Element number
            atemp = alpha[num]; // Element area ratio

            // Read in grid node dof of current element
            tdof[0] = NUM_DOF * Number[n][m].a;
            tdof[1] = tdof[0] + 1;
            tdof[2] = NUM_DOF * Number[n][m].b;
            tdof[3] = tdof[2] + 1;
            tdof[4] = NUM_DOF * Number[n][m].c;
            tdof[5] = tdof[4] + 1;
            tdof[6] = NUM_DOF * Number[n][m].d;
            tdof[7] = tdof[6] + 1;

            // for each load case, add loading contribution
            for(j=0;j<3;j++)
            {
                o = numDof*j;
                o2 = 8*j;
                for(i=0;i<8;i++)
                {
                    hom_force[tdof[i]+o] += fb[i+o2]*atemp;
                }
            }
        }
    }

    // --- solve each load case --- //

    // add constraint equations
    add_sp_mat(K1, CE);

    // remove fixed dof
    rem_sp_mat(K1, fixed[0], 1);

    // remove fixed dof from hom_force
    int con_dof = 2*(elemY+elemX+2); // extra dof (Lagrange mult) for constraint equns
    int con_free = freeDof + con_dof;
    int con_numDof = numDof +con_dof;
    double *hom_disp = calloc(3*con_numDof,sizeof(double));

    for(i=0;i<numDof;i++) // for all dof
    {
        // first set of BCs
        num = fixed[0][i];
        if(num > -1) // if dof not fixed
        {
            for(j=0;j<3;j++) // for all 3 load cases
            {
                n = numDof * j; // point to location in hom_force
                m = con_free * j; // point to location in hom_disp
                hom_disp[m+num] = hom_force[n+i]; // copy load accross
            }
        }
    }
    free(hom_force);
    // solve unit load cases
    FE_Solve(K1, fixed[0], hom_disp, con_free, numDof, 3); // All 3 at once

    //OutDispVTK(inMesh, 3, hom_disp, 0, 0, 0, "hom_disp");

    // --- compute the homogenized stiffness tensor --- //

    double stn[9];
    double Edisp[24];
    double Ehom[9];
    for(i=0;i<9;i++){homMat->mat[i] = 0.0;} // reset

    // For All Elements
    for(m=0;m<elemY;m++)
    {
        for(n=0;n<elemX;n++)
        {
            num = Number[n][m].n; // Element number
            atemp = alpha[num]; // Element area ratio

            if(atemp > 0.0)
            {
                // Read in grid node dof of current element
                tdof[0] = NUM_DOF * Number[n][m].a;
                tdof[1] = tdof[0] + 1;
                tdof[2] = NUM_DOF * Number[n][m].b;
                tdof[3] = tdof[2] + 1;
                tdof[4] = NUM_DOF * Number[n][m].c;
                tdof[5] = tdof[4] + 1;
                tdof[6] = NUM_DOF * Number[n][m].d;
                tdof[7] = tdof[6] + 1;

                // get element displacement arrays
                for(j=0;j<3;j++) // each case
                {
                    temp = j*8;
                    temp2 = j*numDof;
                    for(i=0;i<8;i++) // each dof
                    {
                        Edisp[temp + i] =  hom_disp[temp2 + tdof[i]];
                    }
                }

                // get element (centre) strains for each load case: stn = B x disp
                cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, 3, 3, 8, 0.5, BHom, 8, Edisp, 8, 0.0, stn, 3);
                // 1/2 factor for proper B matrix

                // multiply material property matrix with strains
                cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, 3, 3, 3, 1.0, mat1->mat, 3, stn, 3, 0.0, Ehom, 3);

                // subtract from elem mat matrix & multiply by area ratio and sum into homogenized values
                // do not include some entires as shear test case BCs mess them up
                //for(i=0;i<5;i++){j=Esymm_map[i]; homMat->mat[j] += atemp * (mat1->mat[j] - Ehom[j]);}
                for(i=0;i<9;i++){homMat->mat[i] += atemp * (mat1->mat[i] - Ehom[i]);}
            }
        }
    }

    // also compute other homogenized properties
    for(i=0;i<9;i++){homMat->mat[i] /= area;}
    double E1, E2, v12, v21, G, B;
    E1 = homMat->mat[0] - ((homMat->mat[1] * homMat->mat[1]) / homMat->mat[4]);
    E2 = homMat->mat[4] - ((homMat->mat[1] * homMat->mat[1]) / homMat->mat[0]);
    v12 = homMat->mat[1] / homMat->mat[0];
    v21 = homMat->mat[1] / homMat->mat[4];
    G = homMat->mat[8];
    B = 0.5*(0.5*(homMat->mat[0] + homMat->mat[4]) + homMat->mat[1]);

    // save data
    homMat->g = G;
    homMat->k = B;
    homMat->e = E1;
    homMat->v = v12;

    if(pinfo==3)
    {
        printf("\nE[]:");
        for(i=0;i<9;i+=3)
        {
            printf("\n%12.4e\t%12.4e\t%12.4e",homMat->mat[i],homMat->mat[i+1],homMat->mat[i+2]);
        }

        printf("\nE1=%12.4e\nE2=%12.4e\nv12=%12.4e\nv21=%12.4e\nG=%12.4e\nB=%12.4e",E1, E2, v12, v21, G, B);
    }

    // process displacement fields for sensitivity computation
    int NumNodes = inMesh->NumNodes;
    Coord * NodeCoord = inMesh->NodeCoord;

    // unit strain in x, y and shear
    temp = 2*numDof;
    for(i=0;i<NumNodes;i++)
    {
        j = i*2;
        disp[j] = NodeCoord[i].x;
        disp[j+1] = 0.0;
        disp[numDof + j] = 0.0;
        disp[numDof + j + 1] = NodeCoord[i].y;
        disp[temp + j] = NodeCoord[i].y;
        disp[temp + j + 1] = NodeCoord[i].x;
    }

    // difference in disp fields
    temp = 3 * numDof;
    for(i=0;i<temp;i++) { disp[i] -= hom_disp[i]; }
    free(hom_disp);
}

// function to compute homogenized BCs - direct method
void CHakHomogenize::hom_bc(mesh *inMesh, int **fixed, int *freeDof)
{
    // read in data
    double tol = inMesh->tol;
    int NumNodes = inMesh->NumNodes;
    double maxX = inMesh->maxX - tol;
    double maxY = inMesh->maxY - tol;
    Coord * NodeCoord = inMesh->NodeCoord;

    int i;
    int ycnt=0;
    int xcnt=0;
    Coord loc;
    int *ynodes = malloc(NumNodes*sizeof(int));
    int *xnodes = malloc(NumNodes*sizeof(int));

    for(i=0;i<NumNodes;i++)
    {
        loc.x = NodeCoord[i].x;
        loc.y = NodeCoord[i].y;

        // find top and bottom nodes (y=0 or maxY)
        if(loc.y < tol || loc.y > maxY){ynodes[ycnt++] = i;}

        // find left and right set of nodes (x=0 or maxX)
        if(loc.x < tol || loc.x > maxX){xnodes[xcnt++] = i;}
    }

    // First BC
    bool *fixdof_temp = calloc(2*NumNodes, sizeof(bool));
    for(i=0;i<xcnt;i++){fixdof_temp[xnodes[i]*2] = true;} // fix x dof
    for(i=0;i<ycnt;i++){fixdof_temp[ynodes[i]*2 + 1] = true;} // fix y dof

    // use fixdof_temp to create a dof map (saves time later)
    int allDof = NUM_DOF * inMesh->NumNodes; // total dof
    int *map1 = malloc(allDof * sizeof(int));
    int nd = 0; // count free dof
    for(i=0;i<allDof;i++)
    {
        if(fixdof_temp[i]) { map1[i] = -1; } // fixed dof
        else { map1[i] = nd;  // free dof
               nd++; } // move to next dof
    }
    free(fixdof_temp);

    *freeDof = nd; // total free dof
    fixed[0] = map1;

    // Second BC
    fixdof_temp = calloc(2*NumNodes, sizeof(bool));
    for(i=0;i<xcnt;i++){fixdof_temp[xnodes[i]*2 + 1] = true;} // fix y dof
    for(i=0;i<ycnt;i++){fixdof_temp[ynodes[i]*2] = true;} // fix x dof

    // use fixdof_temp to create a dof map (saves time later)
    int *map2 = malloc(allDof * sizeof(int));
    nd = 0; // count free dof
    for(i=0;i<allDof;i++)
    {
        if(fixdof_temp[i]) { map2[i] = -1; } // fixed dof
        else { map2[i] = nd;  // free dof
            nd++; } // move to next dof
    }
    free(fixdof_temp);

    fixed[1] = map2;

    free(xnodes);
    free(ynodes);
}

// function to compute homogenized BCs - Lagrange multiplier method
void CHakHomogenize::hom_bc2(mesh *inMesh, sp_mat *HomBC)
{
    // read data
    int NodeX = inMesh->elemX + 1;
    int NodeY = inMesh->elemY + 1;
    int numDof = NUM_DOF*NodeX*NodeY; // total dof
    int **Nodes2 = inMesh->Nodes2;

    // total number of displacement BCs
    int tot_bc = 2 * (NodeX + NodeY);
    HomBC->ne = 2 * tot_bc; set_sp_mat(HomBC);

    // now search for coincident nodes on opposite boundaries - build constraint equations
    int i;
    int cnt_ent=0; // count number of entries
    int cnt = numDof; // count number of constraint equations (starting from original num dof)

    // left and right edge displacements
    for(i=1;i<=NodeY;i++)
    {
        HomBC->jcn[cnt_ent] = cnt;
        HomBC->irn[cnt_ent] = Nodes2[1][i] * 2; // x-dof on left
        HomBC->A[cnt_ent++] = 1.0;

        HomBC->jcn[cnt_ent] = cnt++;
        HomBC->irn[cnt_ent] = Nodes2[NodeX][i] * 2; // x-dof on right
        HomBC->A[cnt_ent++] = -1.0;

        HomBC->jcn[cnt_ent] = cnt;
        HomBC->irn[cnt_ent] = Nodes2[1][i] * 2 + 1; // y-dof on left
        HomBC->A[cnt_ent++] = 1.0;

        HomBC->jcn[cnt_ent] = cnt++;
        HomBC->irn[cnt_ent] = Nodes2[NodeX][i] * 2 + 1; // y-dof on right
        HomBC->A[cnt_ent++] = -1.0;
    }

    // botton and top edge displacements
    for(i=1;i<=NodeX;i++)
    {
        HomBC->jcn[cnt_ent] = cnt;
        HomBC->irn[cnt_ent] = Nodes2[i][1] * 2; // x-dof on bottom
        HomBC->A[cnt_ent++] = 1.0;

        HomBC->jcn[cnt_ent] = cnt++;
        HomBC->irn[cnt_ent] = Nodes2[i][NodeY] * 2; // x-dof on right
        HomBC->A[cnt_ent++] = -1.0;

        HomBC->jcn[cnt_ent] = cnt;
        HomBC->irn[cnt_ent] = Nodes2[i][1] * 2 + 1; // y-dof on bottom
        HomBC->A[cnt_ent++] = 1.0;

        HomBC->jcn[cnt_ent] = cnt++;
        HomBC->irn[cnt_ent] = Nodes2[i][NodeY] * 2 + 1; // y-dof on right
        HomBC->A[cnt_ent++] = -1.0;
    }

    /*for(i=0;i<HomBC->ne;i++)
    {
        printf("\n%i, %i, %f",HomBC->irn[i],HomBC->jcn[i],HomBC->A[i]);
    }*/
}
