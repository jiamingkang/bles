/*
 * CHakStress.cpp
 *
 *  Created on: 22 Feb 2015
 *      Author: knai20-admin
 */

#include "CHakStress.h"

CHakStress::CHakStress() {
	// TODO Auto-generated constructor stub

}

CHakStress::~CHakStress() {
	// TODO Auto-generated destructor stub
}

//
//  Stress.c
//  BLES5_Hstress
//
//  Created by Christopher Brampton on 23/10/2014.
//  Copyright (c) 2014 Christopher Brampton. All rights reserved.
//

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "ABFG.h"
#include "EMatrix.h"
#include "Solve.h"
#include "Sens.h"
#include "Stress.h"

/*Function to caculate the Pnorm VonMises Stress in the structure and the Maximum stress estimate*/
double CHakStress::PnormStressCalc(mesh *inMesh, isoMat *inMat, double PN, double relax, double *alpha, double *disp, int itt, double *EStress, double *lsf, bool *fixed, double *Cstress)
{
    /*Intialise function varibles*/
    FILE *outfile;                                              /*File varible for output files*/
    char plotname[40];                                          /*variable to change names of plotting output files*/
    int i, j, k, n, m, x, y, count;                             /*loop counts*/
    int num, temp;                                              /*Element Number*/
    int *tnodes;
    double dtemp, ftemp, etemp, PNtemp, StressMaxTrue;			/*Temporary Varribles*/
    double *U;                                                  /*Elemental Displacement Vector*/
    double nx, ny, gax, gay;                                    /*Guass point Co-ordinates*/
    double n1, n2, n3, n4;                                      /*Shape Functions matrix inputs*/
    double dnx1, dnx2, dnx3, dnx4;                              /*Shape Functions Differentiated WRT to X*/
    double dny1, dny2, dny3, dny4;                              /*Shape Functions Differentiated WRT to Y*/
    double StressX, StressY, StressXY, StressVM;                /*Stress vector, axial stress in x,y and the shear stress and the VonMises Stress*/
    double h = inMesh->h;
    double t = inMesh->t;
    double h_fact;
    double cx, cy;
    int NumNodes = inMesh->NumNodes;
    int elemX = inMesh->elemX;
    int elemY = inMesh->elemY;
    int NumElem = inMesh->NumElem;
    Elem **Number = inMesh->Number;
    Coord *NodeCoord = inMesh->NodeCoord;
    double h2 = h*0.5;
    double Arelax,A;
    double e11 = inMat->mat[0];
    double v12 = inMat->mat[1];
    double g33 = inMat->mat[8];
    double lsfGuass;
    int fixedTemp, print;
    double *EStressX = malloc(elemX*elemY*(sizeof(double)));
    double *EStressY = malloc(elemX*elemY*(sizeof(double)));
    double *EStressXY = malloc(elemX*elemY*(sizeof(double)));

    double ga = 0.5773502692;
    double FF = 4.0;
    Pos Po[4] = {{-1,-1},{1,-1},{1,1},{-1,1}};/**/
    /*double ga = 1.0;
    double FF = 1.0;
    Pos Po[1] = {{0,0}};/**/

    tnodes = malloc(4*sizeof(int));
    U = malloc(8*sizeof(double));

    StressMaxTrue = 0;
    PNtemp = 0;

    for(n=0;n<elemX;n++)
    {
        for(m=0;m<elemY;m++)
        {
            /*Get the node and element numbers*/
            num = Number[n][m].n;
            EStress[num]=0;
            EStressX[num] = 0;
            EStressY[num] = 0;
            EStressXY[num] = 0;
            fixedTemp = 0;

            //if((num == 4735)||(num == 4736)||(num == 4737)||(num == 4835)||(num == 4836)||(num == 4935)||(num == 4936)){print = 1;}
            //if((num == 4038)||(num == 4039)||(num == 4138)||(num == 4139)||(num == 4238)||(num == 4239)||(num == 4338)||(num == 4339)||(num == 4438)||(num == 4439)){print = 1;}
            //else{print = 0;}

            tnodes[0] = Number[n][m].a;
            tnodes[1] = Number[n][m].b;
            tnodes[2] = Number[n][m].c;
            tnodes[3] = Number[n][m].d;

            for(i=0;i<4;i++){fixedTemp += fixed[tnodes[i]];}

            if((alpha[num]>0.000001))//&&(fixedTemp<4))
            {
                /*if(print==1)
                {
                    //printf("\n\n num: %i, ( %i, %i )", num, n, m);
                    printf("\n num: %i", num);
                    //printf("\t Area = %f \nU: ", alpha[num]);
                }/**/
                /*Get the elemental displacement vector, 8 enteries for the 8 degrees of freedom*/
                for(i=0;i<4;i++)
                {
                    j = i*2;
                    temp = tnodes[i] * 2;	/*Get the X dof number for this node*/

                    U[j] = disp[temp];		/*Put x displacement into elemental displacement vector*/
                    U[j+1] = disp[temp+1];	/*Put y displacement into elemental displacement vector*/
                    if(print==1)
                    {
                        printf("\n%i\t%f\t%f", num, (NodeCoord[tnodes[i]].x + U[j]), -1*(NodeCoord[tnodes[i]].y + U[j+1]));
                        //printf("\n%i\t%f\t%f", num, (NodeCoord[tnodes[i]].x), (NodeCoord[tnodes[i]].y));
                    }/**/
                }

                if(print==1)
                {
                    printf("\n%i\t%f\t%f\n", num, (NodeCoord[tnodes[0]].x + U[0]), -1*(NodeCoord[tnodes[0]].y + U[1]));
                    //printf("\n%i\t%f\t%f\n", num, (NodeCoord[tnodes[0]].x), (NodeCoord[tnodes[0]].y));
                }/**/

                for(i=0;i<FF;i++)
                {
                    /*local guass point co-ordinates*/
                    nx = ga * Po[i].x;
                    ny = ga * Po[i].y;

                    A = alpha[num]/(h*h);
                    //printf("\nA = %f", A);
                    /*Get the relaxed area value*/
                    Arelax = pow(A, relax-1);/**/
                    Arelax = 1.0;

                    /*Get the shape functions*/
                    n1 = (1-nx)*(1-ny)*0.25;
                    n2 = (1+nx)*(1-ny)*0.25;
                    n3 = (1+nx)*(1+ny)*0.25;
                    n4 = (1-nx)*(1+ny)*0.25;

                    lsfGuass = 0.0;

                    lsfGuass = n1*lsf[tnodes[0]];
                    lsfGuass += n2*lsf[tnodes[1]];
                    lsfGuass += n3*lsf[tnodes[2]];
                    lsfGuass += n4*lsf[tnodes[3]];
                    //if(print==1){printf("\nlsfGuass: %f", lsfGuass);}

                    /*work out global co-ordinates of element center, from co-ords of node 1 and element edge length*/
                    cx = NodeCoord[tnodes[0]].x + h2;
                    cy = NodeCoord[tnodes[0]].y + h2;

                    if(lsfGuass>0.0) /*If the guass point is inside the structure the lsf value will be greater then 0*/
                    {

                        /*Get the shape functions differentiated wrt to x and y*/
                        dnx1 = (-0.5)*(1-ny)*h;
                        dnx2 = (0.5)*(1-ny)*h;
                        dnx3 = (0.5)*(1+ny)*h;
                        dnx4 = (-0.5)*(1+ny)*h;

                        dny1 = (-0.5)*(1-nx)*h;
                        dny2 = (-0.5)*(1+nx)*h;
                        dny3 = (0.5)*(1+nx)*h;
                        dny4 = (0.5)*(1-nx)*h;

                        /*Now we can get the stress vector at the guass point, Stress = EBU*/
                        StressX = A*e11*(dnx1*U[0] + dnx2*U[2] + dnx3*U[4] + dnx4*U[6]) + A*v12*(dny1*U[1] + dny2*U[3] + dny3*U[5] + dny4*U[7]);
                        StressY = A*e11*(dny1*U[1] + dny2*U[3] + dny3*U[5] + dny4*U[7]) + A*v12*(dnx1*U[0] + dnx2*U[2] + dnx3*U[4] + dnx4*U[6]);
                        StressXY = A*g33*(dny1*U[0] + dnx1*U[1] + dny2*U[2] + dnx2*U[3] + dny3*U[4] + dnx3*U[5] + dny4*U[6] + dnx4*U[7]);

                        //if(print==1){printf("\n StressX = %f, StressY = %f, StressXY = %f", StressX, StressY, StressXY);}

                        /*Next we get the Von Mises Stress at the guass points*/
                        dtemp = (StressX*StressX + StressY*StressY - StressX*StressY + 3*StressXY*StressXY);
                        StressVM = sqrt(dtemp);

                        /*Now we relax the stress relative to the element area to stop low density elements having infinate stress and screwing up the model*/

                        EStress[num] += (1.0/FF)*StressVM*Arelax;

                        EStressX[num] += ((1.0/FF)*(StressX));
                        EStressY[num] += ((1.0/FF)*(StressY));
                        EStressXY[num] += ((1.0/FF)*(StressXY));



                        //if(print==1){printf("\n EStress[%i] += 0.25 * %f * %f = %f \tEStress[%i] = %f", num, StressVM, Arelax, 0.25*StressVM*Arelax, num, EStress[num]);}
                    }
                }

                //if(print == 1){printf("\n%i %f %f %f %f", num, EStressX[num], EStressY[num], EStressXY[num],EStress[num]);}

                //printf("\nFixtemp[%i][%i] = %i", n, m, fixedTemp);
                if(fixedTemp!=4)
                {
                    PNtemp += pow(EStress[num], PN);

                    if(EStress[num]>StressMaxTrue)
                    {
                        StressMaxTrue = EStress[num];
                    }
                }
                //else{EStress[num] = 0;}
            }
        }
    }

    /*sprintf(plotname,"Estress%i.txt",itt);
    outfile = fopen(plotname, "w");
    if(outfile == NULL){
        printf("\nFailed to open Node Strain Energy writefile\n");
    }
    else{
        for(m=0;m<elemY;m++)
        {
            for(n=0;n<elemX;n++)
            {
                num = Number[n][m].n;
                //printf("\n (%i, %i): EStress[%i] = %f", n, m, num, EStress[num]);
                fprintf(outfile,"%f\t", EStress[num]);
            }
            fprintf(outfile,"\n");
        }
    }
    fclose(outfile);/**/

    /*now calcualte the Pnorm Stress*/
    dtemp = PNtemp;
    ftemp = 1.0/PN;
    PNtemp = pow(dtemp, ftemp);

    printf("\n PNtemp = %f", PNtemp);


    //Now normalise the pnorm stress
    if(itt==0){Cstress[4] = 1.0;}
    else
    {
        Cstress[3] = Cstress[4];                                                    //Set Cstress[3] to be last iterations corection factor
        Cstress[4] = (Cstress[0]*(Cstress[1]/Cstress[2])) + ((1-Cstress[0])*Cstress[3]);   //Calculate the correction factor in Cstress[4]
        Cstress[1] = StressMaxTrue;                                                 //Set Cstress1 to this Maximum Stress for next iteration
        Cstress[2] = PNtemp;                                                        //Set Cstress2 to this Pnorm Stress for next iteration
    }

    PNtemp = Cstress[4]*PNtemp;         //Calculate normalised Stress
    Cstress[5] = StressMaxTrue;
    printf("\n C = %f", Cstress[4]);
    printf("\n Corrected PNtemp = %f", PNtemp);
    printf("\n True Max Stress = %f", StressMaxTrue);

    char filename[40];
    sprintf(filename,"EstressX");
    OutPLotStresssVTK2(inMesh, lsf, EStressX, 3, itt, filename);
    sprintf(filename,"EstressY");
    OutPLotStresssVTK2(inMesh, lsf, EStressY, 3, itt, filename);
    sprintf(filename,"EstressXY");
    OutPLotStresssVTK2(inMesh, lsf, EStressXY, 3, itt, filename);

    free(EStressX);
    free(EStressY);
    free(EStressXY);
    free(U);
    free(tnodes);
    return(PNtemp);
}

// calculate the loading for the adjoint problem
void CHakStress::StressAdjoint(mesh *inMesh, isoMat *inMat, double *adjoint, double PN, double relax, double *alpha, double *disp, double NumDof, int itt, double *lsf, int *fixDof, double *EStress, bool *fixed)
{
    FILE *outfile;                                              /*File varible for output files*/
    char plotname[40];                                          /*variable to change names of plotting output files*/
    int i, j, k, n, m, x, y, count;                             /*loop counts*/
    int num, temp;                                              /*Element Number*/
    int *tnodes;
    double dtemp, ftemp, etemp, PNtemp, StressMaxTrue;			/*Temporary Varribles*/
    double *U, *Tadjoint;                                       /*Elemental Displacement Vector*/
    double nx, ny, gax, gay;                                    /*Guass point Co-ordinates*/
    double n1, n2, n3, n4;                                      /*Shape Functions matrix inputs*/
    double dnx1, dnx2, dnx3, dnx4;                              /*Shape Functions Differentiated WRT to X*/
    double dny1, dny2, dny3, dny4;                              /*Shape Functions Differentiated WRT to Y*/
    double StressX, StressY, StressXY, StressVM;                /*Stress vector, axial stress in x,y and the shear stress and the VonMises Stress*/
    double h = inMesh->h;
    double cx, cy;
    int NumNodes = inMesh->NumNodes;
    int elemX = inMesh->elemX;
    int elemY = inMesh->elemY;
    int NumElem = inMesh->NumElem;
    Elem **Number = inMesh->Number;
    Coord *NodeCoord = inMesh->NodeCoord;
    double h2 = h*0.5;
    double Arelax,A;
    double e11 = inMat->mat[0];
    double v12 = inMat->mat[1];
    double g33 = inMat->mat[8];
    int fixedTemp;
    double lsfGuass;
    double ga = 0.5773502692;
    double FF = 4;
    Pos Po[4] = {{-1,-1},{1,-1},{1,1},{-1,1}};/**/
    /*double ga = 0.0;
    double FF = 1.0;
    Pos Po[1] = {{0,0}};/**/
    tnodes = malloc(4*sizeof(int));
    U = malloc(8*sizeof(double));
    Tadjoint = malloc(8*sizeof(double));

    //printf("/n\t\t1111AAAAA11111/n fixDof[0] = %i", fixDof[0]);

    /*Right the formulation for each entry in the adjoint load vector is p*(OVM^(p-1))*(DOVM/D)EBU */
    for(i=0;i<NumDof;i++){adjoint[i]=0;}

    for(n=0;n<elemX;n++)
    {
        for(m=0;m<elemY;m++)
        {
            num = Number[n][m].n;
            fixedTemp = 0;
            /*Get the element Nodes*/
            tnodes[0] = Number[n][m].a;
            tnodes[1] = Number[n][m].b;
            tnodes[2] = Number[n][m].c;
            tnodes[3] = Number[n][m].d;
            for(i=0;i<4;i++){fixedTemp += fixed[tnodes[i]];}

            if((alpha[num]>0.000001)&&(EStress[num]>0.00000001))//&&(fixedTemp<4))
            {

                /*Get the elemental displacement vector, 8 enteries for the 8 degrees of freedom*/
                for(i=0;i<4;i++)
                {
                    j = i*2;
                    temp = tnodes[i] * 2;	/*Get the X dof number for this node*/

                    U[j] = disp[temp];		/*Put x displacement into elemental displacement vector*/
                    U[j+1] = disp[temp+1];	/*Put y displacement into elemental displacement vector*/
                    Tadjoint[j] = 0;		/*Initialise this to zero*/
                    Tadjoint[j+1] = 0;
                }

                for(i=0;i<FF;i++)
                {
                    /*local guass point co-ordinates*/
                    nx = ga * Po[i].x*h2;
                    ny = ga * Po[i].y*h2;

                    A = alpha[num]/(h*h);
                    /*Get the relaxed area value*/
                    Arelax = pow(A, relax-1);/**/
                    Arelax = 1;

                    /*Get the shape functions*/
                    n1 = (1-nx)*(1-ny)*0.25;
                    n2 = (1+nx)*(1-ny)*0.25;
                    n3 = (1+nx)*(1+ny)*0.25;
                    n4 = (1-nx)*(1+ny)*0.25;

                    lsfGuass = 0.0;

                    lsfGuass = n1*lsf[tnodes[0]];
                    lsfGuass += n2*lsf[tnodes[1]];
                    lsfGuass += n3*lsf[tnodes[2]];
                    lsfGuass += n4*lsf[tnodes[3]];

                    if(lsfGuass>0.0) /*If the guass point is inside the structure the lsf value will be greater then 0*/
                    {
                        /*Get the shape functions differentiated wrt to x and y*/
                        dnx1 = (-0.5)*(1-ny)*h;
                        dnx2 = (0.5)*(1-ny)*h;
                        dnx3 = (0.5)*(1+ny)*h;
                        dnx4 = (-0.5)*(1+ny)*h;

                        dny1 = (-0.5)*(1-nx)*h;
                        dny2 = (-0.5)*(1+nx)*h;
                        dny3 = (0.5)*(1+nx)*h;
                        dny4 = (0.5)*(1-nx)*h;

                        /*Calcuate the directional stresses for each gauss point*/
                        StressX = A*e11*(dnx1*U[0] + dnx2*U[2] + dnx3*U[4] + dnx4*U[6]) + A*v12*(dny1*U[1] + dny2*U[3] + dny3*U[5] + dny4*U[7]);
                        StressY = A*e11*(dny1*U[1] + dny2*U[3] + dny3*U[5] + dny4*U[7]) + A*v12*(dnx1*U[0] + dnx2*U[2] + dnx3*U[4] + dnx4*U[6]);
                        StressXY = A*g33*(dny1*U[0] + dnx1*U[1] + dny2*U[2] + dnx2*U[3] + dny3*U[4] + dnx3*U[5] + dny4*U[6] + dnx4*U[7]);/**/

                        dtemp = (StressX*StressX + StressY*StressY - StressX*StressY + 3*StressXY*StressXY);
                        StressVM = sqrt(dtemp);
                        /*etemp = (1.0/(2*StressVM));/**/

                        /*Calculate (OVM)^(p-1)*/
                        dtemp = pow((StressVM*Arelax) , PN);
                        //dtemp = pow((StressVM) , PN);
                        ftemp = pow(StressVM , 2);/**/
                        PNtemp = PN*dtemp/ftemp;

                        /*Get this guass point contribution to the adjoint load vector for this element*/
                        Tadjoint[0] += PNtemp*((StressX - 0.5*StressY)*A*e11*dnx1 + (-0.5*StressX + StressY)*A*v12*dnx1 + 3*StressXY*A*g33*dny1);
                        Tadjoint[1] += PNtemp*((StressX - 0.5*StressY)*A*v12*dny1 + (-0.5*StressX + StressY)*A*e11*dny1 + 3*StressXY*A*g33*dnx1);
                        Tadjoint[2] += PNtemp*((StressX - 0.5*StressY)*A*e11*dnx2 + (-0.5*StressX + StressY)*A*v12*dnx2 + 3*StressXY*A*g33*dny2);
                        Tadjoint[3] += PNtemp*((StressX - 0.5*StressY)*A*v12*dny2 + (-0.5*StressX + StressY)*A*e11*dny2 + 3*StressXY*A*g33*dnx2);
                        Tadjoint[4] += PNtemp*((StressX - 0.5*StressY)*A*e11*dnx3 + (-0.5*StressX + StressY)*A*v12*dnx3 + 3*StressXY*A*g33*dny3);
                        Tadjoint[5] += PNtemp*((StressX - 0.5*StressY)*A*v12*dny3 + (-0.5*StressX + StressY)*A*e11*dny3 + 3*StressXY*A*g33*dnx3);
                        Tadjoint[6] += PNtemp*((StressX - 0.5*StressY)*A*e11*dnx4 + (-0.5*StressX + StressY)*A*v12*dnx4 + 3*StressXY*A*g33*dny4);
                        Tadjoint[7] += PNtemp*((StressX - 0.5*StressY)*A*v12*dny4 + (-0.5*StressX + StressY)*A*e11*dny4 + 3*StressXY*A*g33*dnx4);
                    }
                }
                /*Now add this element's adjoint load vector to the global load vector*/
                for(i=0;i<4;i++)
                {
                    j = i*2;
                    temp = tnodes[i] * 2;	/*Get the X dof number for this node*/
                    //printf("\nfixdof[%i] = %i\t fixdof[%i] = %i",temp, fixDof[temp], temp +1, fixDof[temp+1]);

                    if(fixDof[temp] > -1)
                    {
                        adjoint[fixDof[temp]] += Tadjoint[j]*(1.0/(FF));		/*Put x Elemental Load into adjoint vector*/
                    }
                    if(fixDof[temp+1] > -1)
                    {
                        adjoint[fixDof[temp+1]] += Tadjoint[j+1]*(1.0/(FF));	/*Put y Elemental Load into adjoint vector*/
                    }
                }
            }
        }
    }

    free(U);
    free(tnodes);
    free(Tadjoint);
}

// calculate the stress sensitivies using least squares of integration points for AFG method
void CHakStress::AFG_Stress_Sens(mesh *inMesh, boundary *bound_in, double *alpha, isoMat *inMat,  double *Nsens,
                     double *disp, double *StAdj, int StressNumDual, int numCase, double *wgt, Coord *gCoord, double aMin, int mode, double Gspn, double PN, double RN, double *lsf, double *EStress, int *fixDof, bool *fixed, int itt)
{
    FILE *outfile;                                              /*File varible for output files*/
    char plotname[40];
    // read in mesh data
    double h = inMesh->h;
    int NumNodes = inMesh->NumNodes;
    int elemX = inMesh->elemX;
    int elemY = inMesh->elemY;
    int NumElem = inMesh->NumElem;
    Elem **Number = inMesh->Number;
    Coord *NodeCoord = inMesh->NodeCoord;
    int numDual = 1; //Stress Sensitvity interpolated in its own way so is by itself here. StressNumDual points to which sensitvity stress is for Nsens placement.
    double rad = 3.0 * h; // hard coded for two elements around a node
    rad *=rad; // input squared dist to Lsens
    Coord *AuxNodes = bound_in->AuxNodes;
    int numBound = bound_in->NumBound;
    Bseg *Boundary = bound_in->Bound;
    int Ntot = NumNodes + bound_in->NumAux;
    double alMin = (aMin < 1.0e-6) ? 1.0e-6 : aMin; // numerical tollerance for sens calc
    if(mode==1){ alMin = (aMin < 1.0e-2) ? 1.0e-2 : aMin; } // increase for eigenvalue sensitivities
    else if(mode==2){alMin = 0.0;} // for compliant mechanisms
    double *lsfGuass;
    lsfGuass = malloc(4*NumElem*sizeof(double));
    double ga = 0.5773502692;
    Pos Po[4] = {{-1,-1},{1,-1},{1,1},{-1,1}};/**/
    double nx, ny, n1, n2, n3, n4;
    double *ESens;
    ESens = malloc(NumElem*sizeof(double));
    double *ENSens;
    ENSens = malloc(NumNodes*sizeof(double));
    double pi2 = 1.5707963267949;
    double *Trust;
    Trust = malloc(NumElem*sizeof(double));         //Varible to check for elements of jagged edges, defining how much we trust the stress/sensitvity values at the gp
    double Tcount;

    int i,j,n,m,o,num,node,ex,ey,xMin,xMax,yMin,yMax,num2;
    double atemp,w1,w2,cx,cy,dtemp,wtemp,Ttemp;
    if(mode==2){w1 = -1.0/wgt[0]; w2 = wgt[1]*w1*w1;} // extra weights for complinat mech disp constraint
    int Gcount = NumElem * 4;		// varible to track number of points evaluated
    int tnodes[4];
    int gpoints = NumElem * 4; // number of sensitivity values (4 per element)
    double *gSens = calloc(gpoints*numDual,sizeof(double));	// define memory for gauss point sensitivities
    int nCount = 0;
    isoMat *mptr; // pointer for material


    // Step 1. Compute gauss point stress senstivities

    StressSensitivity(inMesh, inMat, PN, RN, alpha, disp, StAdj, itt, EStress, ESens, gCoord, gSens, lsfGuass, lsf, fixDof, fixed, Gspn);


    // Step 2. Compute node sensitivities using least squares method

    double *sens_temp = malloc(numDual * sizeof(double)); // array for sens of each dual state

    // calculate smoothed sensitivities for all boundary (or non-OUT) nodes
    bool *done = calloc(Ntot, sizeof(bool));
    //bool *bc = levelset->bound; // (not used)
    double ftemp;

    double *x1, *y1, *x2, *y2, *Theta, *BL;
    double xtemp1, ytemp1, xtemp2, ytemp2,dx,dy;
    int node1, node2;
    x1 = malloc(Ntot*sizeof(double));
    x2 = malloc(Ntot*sizeof(double));
    y1 = malloc(Ntot*sizeof(double));
    y2 = malloc(Ntot*sizeof(double));
    Theta = malloc(Ntot*sizeof(double));
    BL = malloc(Ntot*sizeof(double));

    for(n=0;n<Ntot;n++)
    {
        x1[n] = -1;
        y1[n] = -1;
        x2[n] = -1;
        y2[n] = -1;
        BL[n] = 0;
    }

    //NOW WE WILL GET THE NORMAL ORIENTATION AT THE BOUNDARY
    // first search through each boundary segment and get the orientation of the boundary
    for(n=0;n<numBound;n++)
    {
        node1=Boundary[n].n1; // 1st node number
        node2=Boundary[n].n2; // 2nd node number
        //For First Boundary Node
        if(node2 < NumNodes) // If attached node is a grid node
        {
            xtemp1 = NodeCoord[node2].x; //Store co-ordinates of attached grid node
            ytemp1 = NodeCoord[node2].y;
        }
        else // Otherwise attached node is an auxillary node
        {
            m = node2 - NumNodes; // position in AuxNodes array
            xtemp1 = AuxNodes[m].x; //Store co-ordinates of attached grid node
            ytemp1 = AuxNodes[m].y;
        }
        //For Second Boundary Node
        if(node1 < NumNodes) // If attached node is a grid node
        {
            xtemp2 = NodeCoord[node1].x; //Store co-ordinates of attached grid node
            ytemp2 = NodeCoord[node1].y;
        }
        else // Otherwise attached node is an auxillary node
        {
            m = node1 - NumNodes; // position in AuxNodes array
            xtemp2 = AuxNodes[m].x; //Store co-ordinates of attached grid node
            ytemp2 = AuxNodes[m].y;
        }
        //Now store these results in the correct array location
        if(x1[node1] == -1){ x1[node1] = xtemp1;    y1[node1] = ytemp1;}
        else if(x2[node1] == -1){ x2[node1] = xtemp1;    y2[node1] = ytemp1;}
        else{printf("\n\n!!!!1__ERROR_______NODE %i Attached to more then two boundary segments________\nx1[%i] = %f\tx2[%i] = %f\ty1[%i] = %f\ty2[%i] = %f\n\n", node1, node1, x1[node1], node1, x2[node1], node1, y1[node1], node1, y2[node1]);}

        if(x1[node2] == -1){ x1[node2] = xtemp2;    y1[node2] = ytemp2;}
        else if(x2[node2] == -1){ x2[node2] = xtemp2;    y2[node2] = ytemp2;}
        else{printf("\n\n!!!!2__ERROR_______NODE %i Attached to more then two boundary segments________\nx1[%i] = %f\tx2[%i] = %f\ty1[%i] = %f\ty2[%i] = %f\n\n", node2, node2, x1[node2], node2, x2[node2], node2, y1[node2], node2, y2[node2]);}
    }

    //Now go through all nodes and calcualte normal.
    for(n=0;n<Ntot;n++)
    {
        if((x1[n]!=-1)&&(x2[n]!=-1))
        {
            dx = (x2[n]<x1[n]) ? (x1[n]-x2[n]):(x2[n]-x1[n]);
            dy = (x2[n]<x1[n]) ? (y1[n]-y2[n]):(y2[n]-y1[n]);
            ftemp = atan(dy/dx);
            if(dy==0){ftemp = 0;}
            else if(dx==0){ftemp = pi2;}
            else{ftemp = atan(dy/dx);}
            //printf("\ndx: %f\t dy: %f\t ftemp: %f", dx, dy, ftemp);

            //Get Boundary Lenght
            if(n < NumNodes) // If attached node is a grid node
            {
                cx = NodeCoord[n].x;
                cy = NodeCoord[n].y;
            }
            else // Otherwise attached node is an auxillary node
            {
                m = n - NumNodes;
                cx = AuxNodes[m].x;
                cy = AuxNodes[m].y;
            }
            xtemp1 = (x1[n]-cx)*(x1[n]-cx);
            xtemp2 = (x2[n]-cx)*(x2[n]-cx);
            ytemp1 = (y1[n]-cy)*(y1[n]-cy);
            ytemp2 = (y2[n]-cy)*(y2[n]-cy);

            BL[n] = sqrt(xtemp1+ytemp1) + sqrt(xtemp2+ytemp2);
            BL[n] = (BL[n]<0.01) ? 0.01:BL[n];
        }
        else if((x1[n]==-1)&&(x2[n]!=-1))
        {
            printf("\n\n!!!!ERROR_______NODE %i Attached to only one boundary segment________\n\n", n);
            if(n < NumNodes) // If attached node is a grid node
            {
                dx = (NodeCoord[n].x<x1[n]) ? (x1[n]-NodeCoord[n].x):(NodeCoord[n].x-x1[n]);
                dy = (NodeCoord[n].x<x1[n]) ? (y1[n]-NodeCoord[n].y):(NodeCoord[n].y-y1[n]);
                if(dy==0){ftemp = 0;}
                else if(dx==0){ftemp = pi2;}
                else{ftemp = atan(dy/dx);}
                //printf("\ndx: %f\t dy: %f\t ftemp: %f", dx, dy, ftemp);
                cx = NodeCoord[n].x;
                cy = NodeCoord[n].y;

                xtemp1 = (x1[n]-cx)*(x1[n]-cx);
                xtemp2 = (x2[n]-cx)*(x2[n]-cx);
                BL[n] = sqrt(xtemp1+ytemp1);
                BL[n] = (BL[n]<0.01) ? 0.01:BL[n];
            }
            else // Otherwise attached node is an auxillary node
            {
                m = n - NumNodes; // position in AuxNodes array
                dx = (AuxNodes[m].x<x1[n]) ? (x1[n]-AuxNodes[m].x):(AuxNodes[m].x-x1[n]);
                dy = (AuxNodes[m].x<x1[n]) ? (y1[n]-AuxNodes[m].y):(AuxNodes[m].y-y1[n]);
                if(dy==0){ftemp = 0;}
                else if(dx==0){ftemp = pi2;}
                else{ftemp = atan(dy/dx);}
                //printf("\ndx: %f\t dy: %f\t ftemp: %f", dx, dy, ftemp);
                cx = AuxNodes[m].x;
                cy = AuxNodes[m].y;

                xtemp1 = (x1[n]-cx)*(x1[n]-cx);
                xtemp2 = (x2[n]-cx)*(x2[n]-cx);
                BL[n] = sqrt(xtemp1+ytemp1);
                BL[n] = (BL[n]<0.01) ? 0.01:BL[n];
            }
        }
        Theta[n] = ftemp-pi2;
    }

    //Establish the trust region
    for(i=0;i<elemX;i++)
    {
        for(j=0;j<elemY;j++)
        {
            num = Number[i][j].n;
            if(alpha[num]>aMin)
            {
                cx = NodeCoord[Number[i][j].a].x + 0.5*h;
                cy = NodeCoord[Number[i][j].a].y + 0.5*h;
                wtemp = 0;
                Ttemp = 0;
                for(n=0;n<elemX;n++)
                {
                    for(m=0;m<elemY;m++)
                    {
                        num2 = Number[n][m].n;
                        dx = cx-(NodeCoord[Number[n][m].a].x + 0.5*h);
                        dy = cy-(NodeCoord[Number[n][m].a].y + 0.5*h);

                        dx *= dx;
                        dy *= dy;

                        dtemp = sqrt(dx+dy);
                        if(dtemp<2)
                        {
                            Ttemp += (2-dtemp)*alpha[num2];
                            wtemp += (2-dtemp);
                        }
                    }
                }
                Trust[num] = Ttemp/wtemp;

            }
            else{Trust[num] = 0;}

        }
    }

    double theta_temp;
    //sprintf(plotname,"Nsens.txt");
    //outfile = fopen(plotname, "w");
    // first search through each boundary segment
    for(n=0;n<numBound;n++)
    {
        ey = Boundary[n].e / elemX;
        ex = Boundary[n].e - elemX*ey; // element indices
        xMin = ex-10; xMax=ex+10;
        yMin = ey-10; yMax=ey+10; // search limits
        xMin = (xMin < 0) ? 0 : xMin;
        yMin = (yMin < 0) ? 0 : yMin;
        xMax = (xMax > elemX) ? elemX : xMax;
        yMax = (yMax > elemY) ? elemY : yMax; // adjust search limits

        node=Boundary[n].n1; // 1st node number
        theta_temp = Theta[node];
        if(!done[node]) // if sens not already computed
        {
            ftemp = 1.0;
            if(node < NumNodes) // If node is a grid node
            {
                num = LsensStress(&NodeCoord[node], xMax, xMin, yMax, yMin, alMin, alpha, rad, gCoord, gSens, Number, 4, numDual, sens_temp, lsfGuass, fixed, EStress, theta_temp,Trust, h);
                //if(bc[node]){ftemp = -1.0;}
            }
            else // Otherwise node is an auxillary node
            {
                m = node - NumNodes; // position in AuxNodes array
                num = LsensStress(&AuxNodes[m], xMax, xMin, yMax, yMin, alMin, alpha, rad, gCoord, gSens, Number, 4, numDual, sens_temp, lsfGuass, fixed, EStress, theta_temp,Trust, h);
            }
            //for(i=0;i<numDual;i++)
            i = StressNumDual;
            {
                m = Ntot*i + node; // point to correct place in Nsens
                Nsens[m] = ftemp * sens_temp[i];// * BL[node]; // multiply smoothed sensitivities by weight
                //fprintf(outfile,"%i\t%0.7f\n", m, Nsens[m]);
            }
            done[node] = true;
            nCount++;
        }
        node=Boundary[n].n2; // 2nd node number
        theta_temp = Theta[node];
        if(!done[node]) // if sens not already computed
        {
            ftemp = 1.0;
            if(node < NumNodes) // If node is a grid node
            {
                num = LsensStress(&NodeCoord[node], xMax, xMin, yMax, yMin, alMin, alpha, rad, gCoord, gSens, Number, 4, numDual, sens_temp, lsfGuass, fixed, EStress, theta_temp,Trust, h);
                //if(bc[node]){ftemp = -1.0;}
            }
            else // Otherwise node is an auxillary node
            {
                m = node - NumNodes; // position in AuxNodes array
                num = LsensStress(&AuxNodes[m], xMax, xMin, yMax, yMin, alMin, alpha, rad, gCoord, gSens, Number, 4, numDual, sens_temp, lsfGuass, fixed, EStress, theta_temp,Trust, h);
            }
            //for(i=0;i<numDual;i++)
            i = StressNumDual;
            {
                m = Ntot*i + node; // point to correct place in Nsens
                Nsens[m] = ftemp * sens_temp[i];// * BL[node]; // multiply smoothed sensitivities by weight
                //fprintf(outfile,"%i\t%0.7f\n", m, Nsens[m]);
            }
            done[node] = true;
            nCount++;
        }
    }
    //fclose(outfile);

    free(Theta);




    int **Nodes2 = inMesh->Nodes2;
    int NodeX = inMesh->NodeX-2;
    int NodeY = inMesh->NodeY-2;


    int ind;

     sprintf(plotname,"Sensitivity%i.vtk",itt);
     outfile = fopen(plotname, "w");
     if(outfile == NULL){
     printf("\nFailed to open Gsens writefile\n");
     }
     fprintf(outfile,"# vtk DataFile Version 3.0\n");
     fprintf(outfile,"Para0\n");
     fprintf(outfile,"ASCII\n");
     fprintf(outfile,"DATASET RECTILINEAR_GRID\n");
     fprintf(outfile,"DIMENSIONS %i %i %i\n", (2*NodeX)-1, (2*NodeY)-1, 1);
     fprintf(outfile,"X_COORDINATES %i double\n", (2*NodeX)-1);
     for(i=0;i<(2*NodeX)-1;i++){fprintf(outfile,"%0.1f ", i/2.0);}
     fprintf(outfile,"\nY_COORDINATES %i double\n", (2*NodeY)-1);
     for(i=0;i<(2*NodeY)-1;i++){fprintf(outfile,"%0.1f ",i/2.0);}
     fprintf(outfile,"\nZ_COORDINATES %i double\n", 1);
     for(i=0;i<1;i++){fprintf(outfile,"%i ",i);}

     fprintf(outfile,"\n\nCELL_DATA %i\n", 4*NumElem);
     fprintf(outfile,"SCALARS Trust double 1\n");
     fprintf(outfile,"LOOKUP_TABLE default\n");

     /*for(n=0;n<NumElem;n++)
     {
     fprintf(outfile,"%f\n", ESens[n]);
     }*/

    for(m=0;m<elemY;m++)
    {
        for(i=0;i<2;i++)
        {
            for(n=0;n<elemX;n++)
            {
                ind = 4*n+4*elemX*m;
                //printf("\n%i", ind);
                if(i==0)
                {
                    fprintf(outfile,"%f\n", gSens[ind]);
                    fprintf(outfile,"%f\n", gSens[ind+1]);
                }
                else if(i==1)
                {
                    fprintf(outfile,"%f\n", gSens[ind+3]);
                    fprintf(outfile,"%f\n", gSens[ind+2]);
                }
            }
        }
    }

    fclose(outfile);/**/

    // clear memory
    free(x1);
    free(x2);
    free(y1);
    free(y2);
    free(BL);
    free(done);
    free(gSens);
    free(sens_temp);
    free(lsfGuass);
    free(ESens);
    free(ENSens);
}

//Function to Calculate the Stress Sensitivities in all the guass points
void CHakStress::StressSensitivity(mesh *inMesh, isoMat *inMat, double PN, double relax, double *alpha, double *disp, double *StAdj, int itt, double *EStress, double *ESens, Coord *gCoord, double *gSens, double *lsfGuass, double *lsf, int *fixdof, bool *fixed, double Gspn)
{
    FILE *outfile;                                              /*File varible for output files*/
    char plotname[40];                                          /*variable to change names of plotting output files*/
    int i, j, k, n, m, x, y, count;                             /*loop counts*/
    int num, temp;                                              /*Element Number*/
    int *tnodes;
    double dtemp, ftemp, etemp, PNtemp, StressMaxTrue;			/*Temporary Varribles*/
    double *U, *Tadj;                                       /*Elemental Displacement Vector*/
    double nx, ny, gax, gay;                                    /*Guass point Co-ordinates*/
    double n1, n2, n3, n4;                                      /*Shape Functions matrix inputs*/
    double dnx1, dnx2, dnx3, dnx4;                              /*Shape Functions Differentiated WRT to X*/
    double dny1, dny2, dny3, dny4;                              /*Shape Functions Differentiated WRT to Y*/
    double StressX, StressY, StressXY, StressVM;                /*Stress vector, axial stress in x,y and the shear stress and the VonMises Stress*/
    double h = inMesh->h;
    double cx, cy;
    int NumNodes = inMesh->NumNodes;
    int elemX = inMesh->elemX;
    int elemY = inMesh->elemY;
    int NumElem = inMesh->NumElem;
    Elem **Number = inMesh->Number;
    Coord *NodeCoord = inMesh->NodeCoord;
    double h2 = h*0.5;
    double Arelax, A, StressSens;
    double e11 = inMat->mat[0];
    double v12 = inMat->mat[1];
    double g33 = inMat->mat[8];
    int fixedTemp;
    double ga = 0.5773502692;
    double FF = 4.0;
    Pos Po[4] = {{-1,-1},{1,-1},{1,1},{-1,1}};/**/
    /*double ga = 1.0;
    double FF = 1.0;
    Pos Po[1] = {{0,0}};/**/

    tnodes = malloc(4*sizeof(int));
    U = malloc(8*sizeof(double));
    Tadj = malloc(8*sizeof(double));

    /*Right the formulation for each entry in the adjoint load vector is p*(OVM^(p-1))*(DOVM/D)EBU */
    for(i=0;i<(FF*NumElem);i++){gSens[i]=0;}

    int Gcount = 0;
    double MinG = 10000000000;
    double MaxG = -10000000000;

    for(m=0;m<elemY;m++)
    {
        for(n=0;n<elemX;n++)
        {
            /*Get the node and element numbers*/
            num = Number[n][m].n;
            ESens[num]=0;

            for(i=0;i<FF;i++)
            {
                gSens[Gcount+i] = 0;
                lsfGuass[Gcount+i] = 0;
            }

            fixedTemp = 0;

            tnodes[0] = Number[n][m].a;
            tnodes[1] = Number[n][m].b;
            tnodes[2] = Number[n][m].c;
            tnodes[3] = Number[n][m].d;

            //for(i=0;i<4;i++){fixedTemp += fixed[tnodes[i]];}

            //if(((n==158)||(n==159)) && ((m==39)||(m==40)||(m==38)||(m==41))){fixedTemp=0;}
            //else{fixedTemp=0;}
            //printf("\nn = %i\tm = %i\tfixedTemp = %i",n,m,fixedTemp);

            if((alpha[num]>0.000001)&&(EStress[num]>0.00000001))//&&(fixedTemp<4))
            {
                /*Get the elemental displacement vector, 8 enteries for the 8 degrees of freedom*/
                for(i=0;i<4;i++)
                {
                    j = i*2;
                    temp = tnodes[i] * 2;	/*Get the X dof number for this node*/

                    U[j] = disp[temp];		/*Put x displacement into elemental displacement vector*/
                    U[j+1] = disp[temp+1];	/*Put y displacement into elemental displacement vector*/

                    Tadj[j] = StAdj[temp];
                    Tadj[j+1] = StAdj[temp+1];
                }

                A = alpha[num];/**/

                for(i=0;i<FF;i++)
                {
                    /*local guass point co-ordinates*/
                    nx = ga * Po[i].x*h2;
                    ny = ga * Po[i].y*h2;

                    A = alpha[num]/(h*h);
                    /*Get the relaxed area value*/
                    Arelax = pow(A, relax-1);/**/
                    Arelax = 1;

                    /*Get the shape functions*/
                    n1 = (1-nx)*(1-ny)*0.25;
                    n2 = (1+nx)*(1-ny)*0.25;
                    n3 = (1+nx)*(1+ny)*0.25;
                    n4 = (1-nx)*(1+ny)*0.25;

                    lsfGuass[Gcount] = 0.0;

                    lsfGuass[Gcount] = n1*lsf[tnodes[0]];
                    lsfGuass[Gcount] += n2*lsf[tnodes[1]];
                    lsfGuass[Gcount] += n3*lsf[tnodes[2]];
                    lsfGuass[Gcount] += n4*lsf[tnodes[3]];
                    /*printf("\nlsfGuass: %f", lsfGuass);/**/

                    /*work out global co-ordinates of element center, from co-ords of node 1 and element edge length*/
                    cx = NodeCoord[tnodes[0]].x + h2;
                    cy = NodeCoord[tnodes[0]].y + h2;

                    if(lsfGuass[Gcount]>0.0) /*If the guass point is inside the structure the lsf value will be greater then 0*/
                    {

                        /*Get the shape functions differentiated wrt to x and y*/
                        dnx1 = (-0.5)*(1-ny)*h;
                        dnx2 = (0.5)*(1-ny)*h;
                        dnx3 = (0.5)*(1+ny)*h;
                        dnx4 = (-0.5)*(1+ny)*h;

                        dny1 = (-0.5)*(1-nx)*h;
                        dny2 = (-0.5)*(1+nx)*h;
                        dny3 = (0.5)*(1+nx)*h;
                        dny4 = (0.5)*(1-nx)*h;

                        /*Now we can get the stress vector at the guass point, Stress = EBU*/
                        StressX = A*e11*(dnx1*U[0] + dnx2*U[2] + dnx3*U[4] + dnx4*U[6]) + A*v12*(dny1*U[1] + dny2*U[3] + dny3*U[5] + dny4*U[7]);
                        StressY = A*e11*(dny1*U[1] + dny2*U[3] + dny3*U[5] + dny4*U[7]) + A*v12*(dnx1*U[0] + dnx2*U[2] + dnx3*U[4] + dnx4*U[6]);
                        StressXY = A*g33*(dny1*U[0] + dnx1*U[1] + dny2*U[2] + dnx2*U[3] + dny3*U[4] + dnx3*U[5] + dny4*U[6] + dnx4*U[7]);

                        /*Next we get the Von Mises Stress at the guass points*/
                        dtemp = (StressX*StressX + StressY*StressY - StressX*StressY + 3*StressXY*StressXY);
                        StressVM = Arelax*sqrt(dtemp);
                        dtemp = pow(fabs(StressVM),PN);

                        StressSens = 0;
                        /*Now calculate u*k/adj*/
                        StressSens  = (((U[0]*dnx1 + U[2]*dnx2 + U[4]*dnx3 + U[6]*dnx4)*e11 + (U[1]*dny1 + U[3]*dny2 + U[5]*dny3 + U[7]*dny4)*v12)*dnx1 + (U[0]*dny1 + U[1]*dnx1 + U[2]*dny2 + U[3]*dnx2 + U[4]*dny3 + U[5]*dnx3 + U[6]*dny4 + U[7]*dnx4)*g33*dny1) * Tadj[0];
                        StressSens += (((U[0]*dnx1 + U[2]*dnx2 + U[4]*dnx3 + U[6]*dnx4)*v12 + (U[1]*dny1 + U[3]*dny2 + U[5]*dny3 + U[7]*dny4)*e11)*dny1 + (U[0]*dny1 + U[1]*dnx1 + U[2]*dny2 + U[3]*dnx2 + U[4]*dny3 + U[5]*dnx3 + U[6]*dny4 + U[7]*dnx4)*g33*dnx1) * Tadj[1];
                        StressSens += (((U[0]*dnx1 + U[2]*dnx2 + U[4]*dnx3 + U[6]*dnx4)*e11 + (U[1]*dny1 + U[3]*dny2 + U[5]*dny3 + U[7]*dny4)*v12)*dnx2 + (U[0]*dny1 + U[1]*dnx1 + U[2]*dny2 + U[3]*dnx2 + U[4]*dny3 + U[5]*dnx3 + U[6]*dny4 + U[7]*dnx4)*g33*dny2) * Tadj[2];
                        StressSens += (((U[0]*dnx1 + U[2]*dnx2 + U[4]*dnx3 + U[6]*dnx4)*v12 + (U[1]*dny1 + U[3]*dny2 + U[5]*dny3 + U[7]*dny4)*e11)*dny2 + (U[0]*dny1 + U[1]*dnx1 + U[2]*dny2 + U[3]*dnx2 + U[4]*dny3 + U[5]*dnx3 + U[6]*dny4 + U[7]*dnx4)*g33*dnx2) * Tadj[3];
                        StressSens += (((U[0]*dnx1 + U[2]*dnx2 + U[4]*dnx3 + U[6]*dnx4)*e11 + (U[1]*dny1 + U[3]*dny2 + U[5]*dny3 + U[7]*dny4)*v12)*dnx3 + (U[0]*dny1 + U[1]*dnx1 + U[2]*dny2 + U[3]*dnx2 + U[4]*dny3 + U[5]*dnx3 + U[6]*dny4 + U[7]*dnx4)*g33*dny3) * Tadj[4];
                        StressSens += (((U[0]*dnx1 + U[2]*dnx2 + U[4]*dnx3 + U[6]*dnx4)*v12 + (U[1]*dny1 + U[3]*dny2 + U[5]*dny3 + U[7]*dny4)*e11)*dny3 + (U[0]*dny1 + U[1]*dnx1 + U[2]*dny2 + U[3]*dnx2 + U[4]*dny3 + U[5]*dnx3 + U[6]*dny4 + U[7]*dnx4)*g33*dnx3) * Tadj[5];
                        StressSens += (((U[0]*dnx1 + U[2]*dnx2 + U[4]*dnx3 + U[6]*dnx4)*e11 + (U[1]*dny1 + U[3]*dny2 + U[5]*dny3 + U[7]*dny4)*v12)*dnx4 + (U[0]*dny1 + U[1]*dnx1 + U[2]*dny2 + U[3]*dnx2 + U[4]*dny3 + U[5]*dnx3 + U[6]*dny4 + U[7]*dnx4)*g33*dny4) * Tadj[6];
                        StressSens += (((U[0]*dnx1 + U[2]*dnx2 + U[4]*dnx3 + U[6]*dnx4)*v12 + (U[1]*dny1 + U[3]*dny2 + U[5]*dny3 + U[7]*dny4)*e11)*dny4 + (U[0]*dny1 + U[1]*dnx1 + U[2]*dny2 + U[3]*dnx2 + U[4]*dny3 + U[5]*dnx3 + U[6]*dny4 + U[7]*dnx4)*g33*dnx4) * Tadj[7];

                        dtemp = pow(fabs(StressVM),PN);

                        //if(itt==56){printf("\n%f\t%f",dtemp, StressSens);}

                        gSens[Gcount] = -0.5*Gspn*((1/A)*dtemp -  StressSens);

                        //if(itt==100){gSens[Gcount] = 0.5*Gspn*(dtemp);}
                        //if(itt==100){gSens[Gcount] = 0.5*Gspn*(StressSens);}
                        ESens[num] += (1/4.0)*gSens[Gcount];

                        MinG = (gSens[Gcount]<MinG) ? gSens[Gcount]:MinG;

                        Gcount++;
                    }
                    else
                    {

                        gSens[Gcount] = 0;
                        Gcount++;
                    }
                }


            }
            else
            {
                for(i=0;i<FF;i++)
                {
                    gSens[Gcount] = 0;
                    Gcount++;
                }
                /*gSens[Gcount] = 0;
                Gcount++;
                gSens[Gcount] = 0;
                Gcount++;
                gSens[Gcount] = 0;
                Gcount++;/**/
            }

        }
    }

    //for(i=0;i<Gcount;i++){gSens[i] -= MinG;}

    free(tnodes);
    free(U);
    free(Tadj);
}

// Function that calculates the stress sensitivity of a node by a least squares (2nd order) filter of near-by gauss points
int CHakStress::LsensStress(Coord *pt, int xMax, int xMin, int yMax, int yMin, double aMin, double *alpha, double r2,
                Coord *gCoord, double *gSens, Elem **Number, int wFlag, int numDual, double *out, double *lsfGuass, bool *fixed, double *EStress, double theta, double *Trust, double h)
{
    int i,j,n,m,num,ind;  // Incrementor
    double atemp,ftemp,x2,y2,xb,yb,xa,ya;  // variables for squared co-ordinates
    double wgt;	 // varible for weight function
    int pts[500]; // array of global gPoint numbers
    double dist[500]; // array of distances
    double lsf[500]; // array of distances
    double aPts[500]; // array of area ratios (of the points)
    double r2isle;
    ftemp = sqrt(r2);
    atemp = ftemp/2.0;
    //r2isle = atemp*atemp;
    double Ra = 4.0*h;
    double Rb = 2.0*h;
    double Ra2 = Ra*Ra;
    double Rb2 = Rb*Rb;
    double da2, db2, dx, da, db, dy;
    int FF = 4;
    r2isle = 2.0*h;
    int fixedTemp;
    int *tnodes;
    tnodes = malloc(4*sizeof(int));

    int count = 0;	 // initialize count of number of points used to zero
    int IslandCount = 0;

    xa = pt->x;
    ya = pt->y;

    int print = 0;
    //if((xa<40.01)&&(xa>36.99)&&(ya<48.01)&&(ya>42.99)){print = 1;}

    /*yMin = 0;
     yMax = 80;
     xMin = 0;
     xMax = 160;/**/

    if(print==1)
    {
        //printf("\n\nTheta = %f", theta);
        printf("\nNode\t %f\t%f", xa, ya);
    }

    // only look at nearby elements
    for(m=yMin;m<yMax;m++)
    {
        //printf("\nm: %i", m);
        for(n=xMin;n<xMax;n++)
        {
            //printf("\tn: %i\n", n);
            num = Number[n][m].n;	// Element number
            atemp = alpha[num];		// area ratio

            fixedTemp = 0;

            tnodes[0] = Number[n][m].a;
            tnodes[1] = Number[n][m].b;
            tnodes[2] = Number[n][m].c;
            tnodes[3] = Number[n][m].d;

            for(i=0;i<4;i++){fixedTemp += fixed[tnodes[i]];}

            // check element is not OUT
            if((atemp > aMin)&&(EStress[num]>0.00000001)&&(fixedTemp<4))
                //if((atemp > 0.5)&&(EStress[num]>0.00000001)&&(fixedTemp<4))
            {
                ind = FF*num; // point to correct place in gCoord
                for(i=0;i<FF;i++) // for all gauss points
                {
                    ftemp = pt->x - gCoord[ind].x;
                    x2 = ftemp * ftemp;	 // squared x-distance from current gauss point to node of interest
                    ftemp = pt->y - gCoord[ind].y;
                    y2 = ftemp * ftemp;	 // squared y-distance from current gauss point to node of interest
                    dx = sqrt(x2);
                    dy = sqrt(y2);

                    da = dx*cos(theta) + dy*sin(theta);             //Transfer co-ordinate into terms of the the primary ellipse axis
                    db = -1*dx*sin(theta) + dy*cos(theta);          //That is normal to the structural boundary.

                    da2 = da*da;
                    db2 = db*db;

                    ftemp = (da2/Ra2) + (db2/Rb2);  // The radius formulation of the ellipise

                    if((ftemp < 1)&&(lsfGuass[ind]>0.0)) // if point is inside ellipse, and guass point inside the structure then add data to arrays
                    {

                        ftemp = x2+y2;
                        dist[count] = (Ra-sqrt(ftemp))/Ra; // save distance
                        //dist[count] = sqrt(ftemp);
                        aPts[count] = atemp; // save area ratio
                        //aPts[count] = Trust[num]; // save area ratio
                        lsf[count] = lsfGuass[ind]<1.5 ? 0.666666*lsfGuass[ind]:1.0;
                        pts[count++] = ind;  // save gPoint number
                        if(ftemp < r2isle)
                        {
                            IslandCount++;
                        }
                    }
                    ind++; // next gPoint
                }
            }
        }
    }
    //printf("\nA\nA");
    /* if((xa == 0)&&((ya==1)||(ya==79)))
     {/**/
     //printf("\n%f\t%f", xa, ya);
     /*printf("\ncount: %i\t IslandCount: %i)", count, IslandCount);
     printf("\nx: Min/Max (%i, %i)", xMin, xMax);
     printf("\ny: Min/Max (%i, %i)", yMin, yMax);
     }/**/

    if((count < 8)||(IslandCount < 0)) // if not enough points within radius then point must be on an island!
    {
        //printf("\nWarning island node found: count= %i \t Islecount = %i \t ", count, IslandCount);
        for(j=0;j<numDual;j++)
        {
            ftemp = -10.0;
            // compute average
            /*for(i=0;i<count;i++)
             {
             m = numDual*pts[i]; // point to correct place in gSens
             ftemp += gSens[m + j] * aPts[i]; // still weight ny area ratio
             }*/
            out[j] = ftemp;  // return average
            //printf("\t%12.4e ",ftemp);
        }
        return(1);
    }

    double *A = malloc( 6 * count * sizeof(double)); // array to store gpoint values
    double *B = malloc(numDual * count * sizeof(double)); // array to store rhs' of equation

    double wtemp = 0;
    double dtemp = 0;

    double MaxSens = -10000000000;
    double MinSens = 10000000000;

    for(i=0;i<count;i++)
    {
        ind = pts[i]; // gPoint number
        //if(dist[i]>0.5)
        {
            MaxSens = (gSens[ind]>MaxSens) ? gSens[ind]:MaxSens;
            MinSens = (gSens[ind]<MinSens) ? gSens[ind]:MinSens;
        }
        // calculate the weight function
        switch (wFlag)
        {
            case 1:
                wgt = 1.0; // no weighting
                break;
            case 2:
                wgt = 1.0 / sqrt(dist[i]); // weighted by inverse distance
                break;
            case 3:
                wgt = sqrt(aPts[i]); // weighted by area
                break;
            case 4:
                //wgt = sqrt(aPts[i] / dist[i]); // weighted by area & inverse distance
                //wgt = pow(aPts[i], (1.0/6.0)) / pow(dist[i], (1.0/6.0)); /*weighted by area & inverse distance*/
                wgt = pow(aPts[i], (1.0/6.0)) * pow(dist[i], (1.0/6.0)); /*weighted by area & inverse distance*/
                //wgt = aPts[i] * dist[i] * lsf[i]; /*weighted by area & inverse distance*/

                //wgt = pow(aPts[i], (1.0/2.0)) / pow(dist[i], (1.0)); /*weighted by area & inverse distance*/

                //wgt = (pow(aPts[i], (1/6.0))*pow(lsf[i], (1/2.0))) / pow(dist[i], (1.0/2.0)); /*weighted by area & inverse distance*/
                if(print==1)
                {
                    printf("\n%f\t%f\t%f", gCoord[ind].x, gCoord[ind].y, gSens[ind]);
                    //printf("\n%f\t%f\t%f\t%f\t%f\t%f\t%f", gCoord[ind].x, gCoord[ind].y, lsfGuass[ind], gSens[ind], aPts[i], dist[i], wgt);
                }
                break;
            default:
                printf("\nERROR!! wFlag out of range!");
        }

        xb = gCoord[ind].x - pt->x; // relative x-coord
        yb = gCoord[ind].y - pt->y; // relative y-coord

        A[i] = wgt;
        A[count + i] = xb * wgt;
        A[(2 * count) + i] = yb * wgt;
        A[(3 * count) + i] = xb * yb * wgt;
        A[(4 * count) + i] = xb * xb * wgt;
        A[(5 * count) + i] = yb * yb * wgt;

        //n = numDual*i; // point to correct place in B
        m = numDual*pts[i]; // point to correct place in gSens
        for(j=0;j<numDual;j++)
        {
            n = j * count;
            B[n + i] = gSens[m + j] * wgt; // save senstivity for each dual case
            dtemp += gSens[m + j] * wgt;
            wtemp += wgt;
        }
    }

    // solve least squares problem using LAPACK routine
    d_lsLPK(6, count, numDual, A, B);

    // Finally evaluate sensitivity at the node using the co-efficients
    for(j=0;j<numDual;j++)
    {
        //out[j] = dtemp/wtemp;
        out[j] = B[count*j];
        if(out[j]<MinSens){out[j]=MinSens;}
        else if(out[j]>MaxSens){out[j]=MaxSens;}

        //if(out[j] > 1.0e3)
        //{ printf("sens %i = %12.4e",j+1,out[j]); }
        //if(print==1){printf("\t%12.4e",out[j]);}
    }

    free(A);
    free(B);
    return(0);
}



// calculate the stress sensitivies using least squares of integration points for AFG method
void CHakStress::AFG_Stress_Sens_hole(mesh *inMesh, boundary *bound_in, double *alpha, isoMat *inMat,  double *Nsens, double *disp, double *StAdj, int StressNumDual, int numCase, double *wgt, Coord *gCoord, double aMin, int mode, double Gspn, double PN, double RN, double *lsf, double *EStress, int *fixDof, bool *fixed, int itt, int h_count, int *h_index, int *h_EmapX, int *h_EmapY, int *h_posN, int *h_posE, double *h_Esens, double *h_Nsens, int ND)
{
    int FF = 4;
    FILE *outfile;                                              /*File varible for output files*/
    char plotname[40];
    // read in mesh data
    double h = inMesh->h;
    int NumNodes = inMesh->NumNodes;
    int elemX = inMesh->elemX;
    int elemY = inMesh->elemY;
    int NumElem = inMesh->NumElem;
    Elem **Number = inMesh->Number;
    Coord *NodeCoord = inMesh->NodeCoord;
    int numDual = 1; //Stress Sensitvity interpolated in its own way so is by itself here. StressNumDual points to which sensitvity stress is for Nsens placement.
    double rad = 3.0 * h; // hard coded for two elements around a node
    rad *=rad; // input squared dist to Lsens
    Coord *AuxNodes = bound_in->AuxNodes;
    int numBound = bound_in->NumBound;
    Bseg *Boundary = bound_in->Bound;
    int Ntot = NumNodes + bound_in->NumAux;
    double alMin = (aMin < 1.0e-6) ? 1.0e-6 : aMin; // numerical tollerance for sens calc
    if(mode==1){ alMin = (aMin < 1.0e-2) ? 1.0e-2 : aMin; } // increase for eigenvalue sensitivities
    else if(mode==2){alMin = 0.0;} // for compliant mechanisms
    double *lsfGuass;
    lsfGuass = malloc(4*NumElem*sizeof(double));
    double ga = 0.5773502692;
    Pos Po[4] = {{-1,-1},{1,-1},{1,1},{-1,1}};/**/
    double nx, ny, n1, n2, n3, n4;
    double *ESens;
    ESens = malloc(NumElem*sizeof(double));
    double *ENSens;
    ENSens = malloc(NumNodes*sizeof(double));
    double pi2 = 1.5707963267949;
    double *Trust;
    Trust = malloc(NumElem*sizeof(double));         //Varible to check for elements of jagged edges, defining how much we trust the stress/sensitvity values at the gp
    double Tcount;

    int i,j,n,m,o,num,node,ex,ey,xMin,xMax,yMin,yMax,num2;
    double atemp,w1,w2,cx,cy,dtemp,wtemp,Ttemp;
    if(mode==2){w1 = -1.0/wgt[0]; w2 = wgt[1]*w1*w1;} // extra weights for complinat mech disp constraint
    int Gcount = NumElem * 4;		// varible to track number of points evaluated
    int tnodes[4];
    int gpoints = NumElem * 4; // number of sensitivity values (4 per element)
    double *gSens = calloc(gpoints*numDual,sizeof(double));	// define memory for gauss point sensitivities
    int nCount = 0;
    isoMat *mptr; // pointer for material


    // Step 1. Compute gauss point stress senstivities
    StressSensitivity(inMesh, inMat, PN, RN, alpha, disp, StAdj, itt, EStress, ESens, gCoord, gSens, lsfGuass, lsf, fixDof, fixed, Gspn);


    // Step 2. Compute node sensitivities using least squares method

    double *sens_temp = malloc(numDual * sizeof(double)); // array for sens of each dual state

    // calculate smoothed sensitivities for all boundary (or non-OUT) nodes
    bool *done = calloc(Ntot, sizeof(bool));
    //bool *bc = levelset->bound; // (not used)
    double ftemp;

    double *x1, *y1, *x2, *y2, *Theta, *BL;
    double xtemp1, ytemp1, xtemp2, ytemp2,dx,dy;
    int node1, node2;
    x1 = malloc(Ntot*sizeof(double));
    x2 = malloc(Ntot*sizeof(double));
    y1 = malloc(Ntot*sizeof(double));
    y2 = malloc(Ntot*sizeof(double));
    Theta = malloc(Ntot*sizeof(double));
    BL = malloc(Ntot*sizeof(double));

    for(n=0;n<Ntot;n++)
    {
        x1[n] = -1;
        y1[n] = -1;
        x2[n] = -1;
        y2[n] = -1;
        BL[n] = 0;
    }

    //NOW WE WILL GET THE NORMAL ORIENTATION AT THE BOUNDARY
    // first search through each boundary segment and get the orientation of the boundary
    for(n=0;n<numBound;n++)
    {
        node1=Boundary[n].n1; // 1st node number
        node2=Boundary[n].n2; // 2nd node number
        //For First Boundary Node
        if(node2 < NumNodes) // If attached node is a grid node
        {
            xtemp1 = NodeCoord[node2].x; //Store co-ordinates of attached grid node
            ytemp1 = NodeCoord[node2].y;
        }
        else // Otherwise attached node is an auxillary node
        {
            m = node2 - NumNodes; // position in AuxNodes array
            xtemp1 = AuxNodes[m].x; //Store co-ordinates of attached grid node
            ytemp1 = AuxNodes[m].y;
        }
        //For Second Boundary Node
        if(node1 < NumNodes) // If attached node is a grid node
        {
            xtemp2 = NodeCoord[node1].x; //Store co-ordinates of attached grid node
            ytemp2 = NodeCoord[node1].y;
        }
        else // Otherwise attached node is an auxillary node
        {
            m = node1 - NumNodes; // position in AuxNodes array
            xtemp2 = AuxNodes[m].x; //Store co-ordinates of attached grid node
            ytemp2 = AuxNodes[m].y;
        }
        //Now store these results in the correct array location
        if(x1[node1] == -1){ x1[node1] = xtemp1;    y1[node1] = ytemp1;}
        else if(x2[node1] == -1){ x2[node1] = xtemp1;    y2[node1] = ytemp1;}
        else{printf("\n\n!!!!1__ERROR_______NODE %i Attached to more then two boundary segments________\nx1[%i] = %f\tx2[%i] = %f\ty1[%i] = %f\ty2[%i] = %f\n\n", node1, node1, x1[node1], node1, x2[node1], node1, y1[node1], node1, y2[node1]);}

        if(x1[node2] == -1){ x1[node2] = xtemp2;    y1[node2] = ytemp2;}
        else if(x2[node2] == -1){ x2[node2] = xtemp2;    y2[node2] = ytemp2;}
        else{printf("\n\n!!!!2__ERROR_______NODE %i Attached to more then two boundary segments________\nx1[%i] = %f\tx2[%i] = %f\ty1[%i] = %f\ty2[%i] = %f\n\n", node2, node2, x1[node2], node2, x2[node2], node2, y1[node2], node2, y2[node2]);}
    }

    //Now go through all nodes and calcualte normal.
    for(n=0;n<Ntot;n++)
    {
        if((x1[n]!=-1)&&(x2[n]!=-1))
        {
            dx = (x2[n]<x1[n]) ? (x1[n]-x2[n]):(x2[n]-x1[n]);
            dy = (x2[n]<x1[n]) ? (y1[n]-y2[n]):(y2[n]-y1[n]);
            ftemp = atan(dy/dx);
            if(dy==0){ftemp = 0;}
            else if(dx==0){ftemp = pi2;}
            else{ftemp = atan(dy/dx);}
            //printf("\ndx: %f\t dy: %f\t ftemp: %f", dx, dy, ftemp);

            //Get Boundary Lenght
            if(n < NumNodes) // If attached node is a grid node
            {
                cx = NodeCoord[n].x;
                cy = NodeCoord[n].y;
            }
            else // Otherwise attached node is an auxillary node
            {
                m = n - NumNodes;
                cx = AuxNodes[m].x;
                cy = AuxNodes[m].y;
            }
            xtemp1 = (x1[n]-cx)*(x1[n]-cx);
            xtemp2 = (x2[n]-cx)*(x2[n]-cx);
            ytemp1 = (y1[n]-cy)*(y1[n]-cy);
            ytemp2 = (y2[n]-cy)*(y2[n]-cy);

            BL[n] = sqrt(xtemp1+ytemp1) + sqrt(xtemp2+ytemp2);
            BL[n] = (BL[n]<0.01) ? 0.01:BL[n];
        }
        else if((x1[n]==-1)&&(x2[n]!=-1))
        {
            printf("\n\n!!!!ERROR_______NODE %i Attached to only one boundary segment________\n\n", n);
            if(n < NumNodes) // If attached node is a grid node
            {
                dx = (NodeCoord[n].x<x1[n]) ? (x1[n]-NodeCoord[n].x):(NodeCoord[n].x-x1[n]);
                dy = (NodeCoord[n].x<x1[n]) ? (y1[n]-NodeCoord[n].y):(NodeCoord[n].y-y1[n]);
                if(dy==0){ftemp = 0;}
                else if(dx==0){ftemp = pi2;}
                else{ftemp = atan(dy/dx);}
                //printf("\ndx: %f\t dy: %f\t ftemp: %f", dx, dy, ftemp);
                cx = NodeCoord[n].x;
                cy = NodeCoord[n].y;

                xtemp1 = (x1[n]-cx)*(x1[n]-cx);
                xtemp2 = (x2[n]-cx)*(x2[n]-cx);
                BL[n] = sqrt(xtemp1+ytemp1);
                BL[n] = (BL[n]<0.01) ? 0.01:BL[n];
            }
            else // Otherwise attached node is an auxillary node
            {
                m = n - NumNodes; // position in AuxNodes array
                dx = (AuxNodes[m].x<x1[n]) ? (x1[n]-AuxNodes[m].x):(AuxNodes[m].x-x1[n]);
                dy = (AuxNodes[m].x<x1[n]) ? (y1[n]-AuxNodes[m].y):(AuxNodes[m].y-y1[n]);
                if(dy==0){ftemp = 0;}
                else if(dx==0){ftemp = pi2;}
                else{ftemp = atan(dy/dx);}
                //printf("\ndx: %f\t dy: %f\t ftemp: %f", dx, dy, ftemp);
                cx = AuxNodes[m].x;
                cy = AuxNodes[m].y;

                xtemp1 = (x1[n]-cx)*(x1[n]-cx);
                xtemp2 = (x2[n]-cx)*(x2[n]-cx);
                BL[n] = sqrt(xtemp1+ytemp1);
                BL[n] = (BL[n]<0.01) ? 0.01:BL[n];
            }
        }
        Theta[n] = ftemp-pi2;
    }

    //Establish the trust region
    for(i=0;i<elemX;i++)
    {
        for(j=0;j<elemY;j++)
        {
            num = Number[i][j].n;
            if(alpha[num]>aMin)
            {
                cx = NodeCoord[Number[i][j].a].x + 0.5*h;
                cy = NodeCoord[Number[i][j].a].y + 0.5*h;
                wtemp = 0;
                Ttemp = 0;
                for(n=0;n<elemX;n++)
                {
                    for(m=0;m<elemY;m++)
                    {
                        num2 = Number[n][m].n;
                        dx = cx-(NodeCoord[Number[n][m].a].x + 0.5*h);
                        dy = cy-(NodeCoord[Number[n][m].a].y + 0.5*h);

                        dx *= dx;
                        dy *= dy;

                        dtemp = sqrt(dx+dy);
                        if(dtemp<2)
                        {
                            Ttemp += (2-dtemp)*alpha[num2];
                            wtemp += (2-dtemp);
                        }
                    }
                }
                Trust[num] = Ttemp/wtemp;

            }
            else{Trust[num] = 0;}

        }
    }

    double theta_temp;
    //sprintf(plotname,"Nsens.txt");
    //outfile = fopen(plotname, "w");
    // first search through each boundary segment
    for(n=0;n<numBound;n++)
    {
        ey = Boundary[n].e / elemX;
        ex = Boundary[n].e - elemX*ey; // element indices
        xMin = ex-10; xMax=ex+10;
        yMin = ey-10; yMax=ey+10; // search limits
        xMin = (xMin < 0) ? 0 : xMin;
        yMin = (yMin < 0) ? 0 : yMin;
        xMax = (xMax > elemX) ? elemX : xMax;
        yMax = (yMax > elemY) ? elemY : yMax; // adjust search limits

        node=Boundary[n].n1; // 1st node number
        theta_temp = Theta[node];
        if(!done[node]) // if sens not already computed
        {
            ftemp = 1.0;
            if(node < NumNodes) // If node is a grid node
            {
                num = LsensStress(&NodeCoord[node], xMax, xMin, yMax, yMin, alMin, alpha, rad, gCoord, gSens, Number, 4, numDual, sens_temp, lsfGuass, fixed, EStress, theta_temp,Trust, h);
                //if(bc[node]){ftemp = -1.0;}
            }
            else // Otherwise node is an auxillary node
            {
                m = node - NumNodes; // position in AuxNodes array
                num = LsensStress(&AuxNodes[m], xMax, xMin, yMax, yMin, alMin, alpha, rad, gCoord, gSens, Number, 4, numDual, sens_temp, lsfGuass, fixed, EStress, theta_temp,Trust, h);
            }
            //for(i=0;i<numDual;i++)
            i = StressNumDual;
            {
                m = Ntot*i + node; // point to correct place in Nsens
                Nsens[m] = ftemp * sens_temp[i] ;//* BL[node]; // multiply smoothed sensitivities by weight
                //fprintf(outfile,"%i\t%0.7f\n", m, Nsens[m]);
            }
            done[node] = true;
            nCount++;
        }
        node=Boundary[n].n2; // 2nd node number
        theta_temp = Theta[node];
        if(!done[node]) // if sens not already computed
        {
            ftemp = 1.0;
            if(node < NumNodes) // If node is a grid node
            {
                num = LsensStress(&NodeCoord[node], xMax, xMin, yMax, yMin, alMin, alpha, rad, gCoord, gSens, Number, 4, numDual, sens_temp, lsfGuass, fixed, EStress, theta_temp,Trust, h);
                //if(bc[node]){ftemp = -1.0;}
            }
            else // Otherwise node is an auxillary node
            {
                m = node - NumNodes; // position in AuxNodes array
                num = LsensStress(&AuxNodes[m], xMax, xMin, yMax, yMin, alMin, alpha, rad, gCoord, gSens, Number, 4, numDual, sens_temp, lsfGuass, fixed, EStress, theta_temp,Trust, h);
            }
            //for(i=0;i<numDual;i++)
            i = StressNumDual;
            {
                m = Ntot*i + node; // point to correct place in Nsens
                Nsens[m] = ftemp * sens_temp[i] ;//* BL[node]; // multiply smoothed sensitivities by weight
                //fprintf(outfile,"%i\t%0.7f\n", m, Nsens[m]);
            }
            done[node] = true;
            nCount++;
        }
    }
    //fclose(outfile);

    free(Theta);

   /* sprintf(plotname,"Nsens.txt");
    outfile = fopen(plotname, "r");
    if(outfile == NULL){
        printf("\nFailed to open Gsens writefile\n");
    }
    else{
        for(i=0;i<nCount;i++){fscanf(outfile, "%i", &m);
            fscanf(outfile, "%lf", &Nsens[m]);}
    }
    fclose(outfile);

    sprintf(plotname,"ESens%i.txt",itt);
    outfile = fopen(plotname, "w");
    if(outfile == NULL){
        printf("\nFailed to open Node Strain Energy writefile\n");
    }
    else{
        for(m=0;m<elemY;m++)
        {
            for(n=0;n<elemX;n++)
            {
                ///Get the node and element numbers
                num = Number[n][m].n;
                //printf("\n (%i, %i): EStress[%i] = %f", n, m, num, EStress[num]);
                fprintf(outfile,"%0.20f\t", ESens[num]);
            }
            fprintf(outfile,"\n");
        }
    }
    fclose(outfile);/**/

    //Get the Hole node sensitvities
    int A,B,C,D,l,k;
    int ind1,ind2,ind3,ind4;
    double *max_sens = malloc(numDual*sizeof(double));
    double *min_sens = malloc(numDual*sizeof(double));

    l = h_posN[ND];
    //for(i=0;i<NumNodes;i++){h_Nsens[l+i] = 0.0;}

    for(i=0;i<h_count;i++)
    {
        n = h_EmapX[i];
        m = h_EmapY[i];

        //printf("\n%i: n: %i, m: %i", i, n, m);

        num = Number[n][m].n;
        A  = Number[n][m].a;
        B  = Number[n][m].b;
        C  = Number[n][m].c;
        D  = Number[n][m].d;

        for(j=0;j<numDual;j++)
        {
            if(FF==4)
            {
                ind1 = (4*num)*numDual+j; // point to place in gSens
                ind2 = (4*num+1)*numDual+j;
                ind3 = (4*num+2)*numDual+j;
                ind4 = (4*num+3)*numDual+j;
                k = h_posE[j];
                l = h_posN[j];
                h_Esens[k+i] = 0.25*(gSens[ind1] + gSens[ind2] + gSens[ind3] + gSens[ind4]);
                if(h_index[A]==1){h_Nsens[l+A] += 0.25*h_Esens[k+i];}
                if(h_index[B]==1){h_Nsens[l+B] += 0.25*h_Esens[k+i];}
                if(h_index[C]==1){h_Nsens[l+C] += 0.25*h_Esens[k+i];}
                if(h_index[D]==1){h_Nsens[l+D] += 0.25*h_Esens[k+i];}
            }
            else if(FF==1)
            {
                 ind1 = (1*num)*numDual+j; // point to place in gSens
                k = h_posE[j];
                l = h_posN[j];
                h_Esens[k+i] = gSens[ind1];
                if(h_index[A]==1){h_Nsens[l+A] += 0.25*h_Esens[k+i];}
                if(h_index[B]==1){h_Nsens[l+B] += 0.25*h_Esens[k+i];}
                if(h_index[C]==1){h_Nsens[l+C] += 0.25*h_Esens[k+i];}
                if(h_index[D]==1){h_Nsens[l+D] += 0.25*h_Esens[k+i];}
            }
        }
    }

    for(j=0;j<numDual;j++)
    {
        k = h_posN[ND];
        max_sens[j] = 0;
        min_sens[j] = 0;

        for(i=0;i<NumNodes;i++)
        {
            if(h_index[i]==1)
            {
                max_sens[j] = ((h_Nsens[k+i])>max_sens[j]) ? (h_Nsens[k+i]):max_sens[j];
                min_sens[j] = ((h_Nsens[k+i])<min_sens[j]) ? (h_Nsens[k+i]):min_sens[j];
            }
        }
        max_sens[j] = max_sens[j]-min_sens[j];
        //printf("\n Maxsens[%i] = %f", j, max_sens[j]);
    }

    for(j=0;j<numDual;j++)
    {
        k = h_posN[j];
        for(i=0;i<NumNodes;i++)
        {
            if(h_index[i]==1)
            {
                h_Nsens[k+i] = h_Nsens[k+i]/max_sens[j];
            }
            else{h_Nsens[k+i] =1;}
        }
    }

    int **Nodes2 = inMesh->Nodes2;
    int NodeX = inMesh->NodeX-2;
    int NodeY = inMesh->NodeY-2;

    //NodeX++;
    //NodeY++;

    int ind;

    sprintf(plotname,"Sensitivity%i.vtk",itt);
    outfile = fopen(plotname, "w");
    if(outfile == NULL){
        printf("\nFailed to open Gsens writefile\n");
    }
    fprintf(outfile,"# vtk DataFile Version 3.0\n");
    fprintf(outfile,"Para0\n");
    fprintf(outfile,"ASCII\n");
    fprintf(outfile,"DATASET RECTILINEAR_GRID\n");
    fprintf(outfile,"DIMENSIONS %i %i %i\n", (2*NodeX)-1, (2*NodeY)-1, 1);
    fprintf(outfile,"X_COORDINATES %i double\n", (2*NodeX)-1);
    for(i=0;i<(2*NodeX)-1;i++){fprintf(outfile,"%0.1f ", i/2.0);}
    fprintf(outfile,"\nY_COORDINATES %i double\n", (2*NodeY)-1);
    for(i=0;i<(2*NodeY)-1;i++){fprintf(outfile,"%0.1f ",i/2.0);}
    fprintf(outfile,"\nZ_COORDINATES %i double\n", 1);
    for(i=0;i<1;i++){fprintf(outfile,"%i ",i);}

    fprintf(outfile,"\n\nCELL_DATA %i\n", 4*NumElem);
    fprintf(outfile,"SCALARS Trust double 1\n");
    fprintf(outfile,"LOOKUP_TABLE default\n");

    /*for(n=0;n<NumElem;n++)
     {
     fprintf(outfile,"%f\n", ESens[n]);
     }*/

    for(m=0;m<elemY;m++)
    {
        for(i=0;i<2;i++)
        {
            for(n=0;n<elemX;n++)
            {
                ind = 4*n+4*elemX*m;
                //printf("\n%i", ind);
                if(i==0)
                {
                    fprintf(outfile,"%f\n", gSens[ind]);
                    fprintf(outfile,"%f\n", gSens[ind+1]);
                }
                else if(i==1)
                {
                    fprintf(outfile,"%f\n", gSens[ind+3]);
                    fprintf(outfile,"%f\n", gSens[ind+2]);
                }
            }
        }
    }

    fclose(outfile);/**/

    // clear memory
    free(max_sens);
    free(min_sens);
    free(x1);
    free(x2);
    free(y1);
    free(y2);
    free(BL);
    free(done);
    free(gSens);
    free(sens_temp);
    free(lsfGuass);
    free(ESens);
    free(ENSens);
}
