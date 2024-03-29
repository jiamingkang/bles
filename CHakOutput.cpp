/*
	CHakOutput.cpp

  	Created on: Nov 24, 2014
	Author: Khalid Ismail, Peter Dunning, JeeHang Lee

 	-- This file is part of Topology Optimisation Opensource Project,
 	owned by MSO (Multidisciplinary and Structural Optimisation) Research Group
 	(http://people.bath.ac.uk/ens/MSORG/index.html) at University of Bath.

 	The project is led by Dr. Hyunsun Alicia Kim
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
#include "CHakOutput.h"

//
// Constructor / Destructor
//
CHakOutput::CHakOutput() {
	// TODO Auto-generated constructor stub

}

CHakOutput::~CHakOutput() {
	// TODO Auto-generated destructor stub
}

//
// Implementation - Interfaces
//
// knai20@bath.ac.uk : OOD version

// function to write element and node numbering info
void CHakOutput::OutNumber(CHakMesh& m_mesh, char *datafile)
{
    // read data
    int m_elemX = m_mesh.m_elemX;
    int elemY = m_mesh.m_elemY;
    Elem **Number = m_mesh.m_pNumber;
    int i,j; // incrementors
    FILE *outfile;	// File varible for output files
    char plotname[120];	// variable to change names of plotting output files
    sprintf(plotname,"%s_Numbering.txt",datafile); // set name for current output file
    outfile = fopen(plotname, "w");
    if(outfile == NULL){
        printf("\nFailed to open Numbering writefile\n");
    }
    else{
        fprintf(outfile,"n\ta\tb\tc\td\n"); // column headings
        for(j=0;j<elemY;j++)
        {
            for(i=0;i<m_elemX;i++)
            {
                fprintf(outfile,"%i\t",Number[i][j].n);
                fprintf(outfile,"%i\t",Number[i][j].a);
                fprintf(outfile,"%i\t",Number[i][j].b);
                fprintf(outfile,"%i\t",Number[i][j].c);
                fprintf(outfile,"%i\t",Number[i][j].d);
                fprintf(outfile,"\n");
            }
        }
    }
    fclose(outfile);
    printf("Numbering info file written\n");
}

// function to write node co-ordinate information file
void CHakOutput::OutNodeCoord(CHakMesh& m_mesh, char *datafile)
{
    // read data
    int NumNodes = m_mesh.m_numNodes;
    Coord *NodeCoord = m_mesh.m_pNodeCoord;
    
    int i; // incrementors
    FILE *outfile;	// File varible for output files
    char plotname[120];	// variable to change names of plotting output files
    
    sprintf(plotname,"%s_NCrd.txt",datafile); // set name for current output file
    outfile = fopen(plotname, "w");
    if(outfile == NULL){
        printf("\nFailed to open Node_Coord writefile\n");
    }
    else{
        for(i=0;i<NumNodes;i++)
        {
            fprintf(outfile,"%i\t",i);
            fprintf(outfile,"%lf\t",NodeCoord[i].x);
            fprintf(outfile,"%lf\t",NodeCoord[i].y);
            fprintf(outfile,"\n");
        }
    }
    fclose(outfile);
    printf("Node Co-ordinate File written");
}

// fucntion to output lsf (& maybe alpha) in vtk format
void CHakOutput::OutPLotShapeVTK2(CHakMesh& m_mesh, double *lsf, double *alpha, int pinfo, int itt, char *datafile)
{
    // read data
    int NodeX = m_mesh.m_nodeX-2;
    int NodeY = m_mesh.m_nodeY-2;
    int NumNodes = m_mesh.m_numNodes;
    int **Nodes2 = m_mesh.m_pNodes2D;
    
    int i,n,m,num;
    
    FILE *outfile;	// File varible for output files
    char plotname[120];	// variable to change names of plotting output files
    
    // Write initial signed distance value information file
    sprintf(plotname,"%s_Lsf_%i.vtk",datafile,itt); // set name for current output file
    outfile = fopen(plotname, "w");
    
    if(outfile == NULL){
        printf("\nFailed to open signed distance paraview writefile\n");
    }
    else
    {
        fprintf(outfile,"# vtk DataFile Version 3.0\n");
        fprintf(outfile,"Para0\n");
        fprintf(outfile,"ASCII\n");
        fprintf(outfile,"DATASET RECTILINEAR_GRID\n");
        fprintf(outfile,"DIMENSIONS %i %i %i\n", NodeX, NodeY, 1);
        fprintf(outfile,"X_COORDINATES %i int\n", NodeX);
        for(i=0;i<NodeX;i++){fprintf(outfile,"%i ", i);}
        fprintf(outfile,"\nY_COORDINATES %i int\n", NodeY);
        for(i=0;i<NodeY;i++){fprintf(outfile,"%i ",i);}
        fprintf(outfile,"\nZ_COORDINATES %i int\n", 1);
        for(i=0;i<1;i++){fprintf(outfile,"%i ",i);}
        
        fprintf(outfile,"\n\nPOINT_DATA %i\n", NumNodes);
        fprintf(outfile,"SCALARS lsf double 1\n");
        fprintf(outfile,"LOOKUP_TABLE default\n");
        
        // Now to print out the node values
        NodeX++;
        NodeY++;
        
        for(m=1;m<NodeY;m++)
        {
            for(n=1;n<NodeX;n++)
            {
                num = Nodes2[n][m];
                fprintf(outfile,"%lf\n", lsf[num]);
            }
        }
        
        // also output element area ratios (if required)
        if(pinfo==3)
        {
            int NumElem = m_mesh.m_numElem;
            fprintf(outfile,"\n\nCELL_DATA %i\n", NumElem);
            fprintf(outfile,"SCALARS alpha double 1\n");
            fprintf(outfile,"LOOKUP_TABLE default\n");
            
            for(n=0;n<NumElem;n++)
            {
                fprintf(outfile,"%lf\n", alpha[n]);
            }
        }
    }
    fclose(outfile);
    printf("\nSigned Distance info file written (Paraview)\n");
}

// function to output boundary as a mesh for Paraview (with shape sensitivities)
void CHakOutput::OutBoundVTK(CHakMesh& m_mesh, CHakBoundary *bound_in, int num_sens, double **Sens, int itt, char *datafile)
{
    // read data
    int NumNodes = m_mesh.m_numNodes;
    Coord *nc = m_mesh.m_pNodeCoord;
    Coord *AuxNodes = bound_in->m_pAuxNodes;
    int NumBound = bound_in->m_numBound;
    Bseg *Boundary = bound_in->Bound;
    int Ntot = NumNodes + bound_in->m_numAux;
    
    int i,n;
    
    FILE *outfile;	// File varible for output files
    char plotname[120];	// variable to change names of plotting output files
    
    // Write initial signed distance value information file
    sprintf(plotname,"%s_Bound-Mesh_%i.vtk",datafile,itt); // set name for current output file
    outfile = fopen(plotname, "w");
    
    if(outfile == NULL){
        printf("\nFailed to open results paraview input writefile\n");
    }
    else
    {
        fprintf(outfile,"# vtk DataFile Version 3.0\n");
        fprintf(outfile,"Para0\n");
        fprintf(outfile,"ASCII\n");
        fprintf(outfile,"DATASET UNSTRUCTURED_GRID\n");
        fprintf(outfile,"POINTS %i float\n",Ntot);
        for(i=0;i<NumNodes;i++)
        {
            fprintf(outfile,"%lf %lf 0.0\n",nc[i].x,nc[i].y);
        }
        n = Ntot-NumNodes;
        for(i=0;i<n;i++)
        {
            fprintf(outfile,"%lf %lf 0.0\n",AuxNodes[i].x,AuxNodes[i].y);
        }
        
        fprintf(outfile,"\nCELLS %i %i\n",NumBound,(NumBound*3));
        for(i=0;i<NumBound;i++)
        {
            fprintf(outfile,"2 %i %i\n",Boundary[i].n1,Boundary[i].n2);
        }
        
        fprintf(outfile,"\nCELL_TYPES %i\n",NumBound);
        for(i=0;i<NumBound;i++)
        {
            fprintf(outfile,"3\n");
        }
        
        fprintf(outfile,"\nPOINT_DATA %i\n",Ntot);
        for(n=0;n<num_sens;n++)
        {
            fprintf(outfile,"SCALARS Sens%i float 1\n",n+1);
            fprintf(outfile,"LOOKUP_TABLE default\n");
            for(i=0;i<Ntot;i++)
            {
                fprintf(outfile,"%12.4e\n",Sens[n][i]);
            }
        }
    }
    
    fclose(outfile);
    printf("\nBoudnary Mesh Info File written (Paraview)\n");
}

// function to output boundary integration data
void CHakOutput::OutBoundInt(int numFunc, int numLbound, int *Lbound_nums, double *Lbound, int itt, char *datafile)
{
    int i,j,p2;
    FILE *outfile;	// File varible for output files
    char plotname[120];	// variable to change names of plotting output files
    
    sprintf(plotname,"%s_Bint_%i.txt",datafile,itt); // set name
    outfile = fopen(plotname, "w");
    if(outfile == NULL){
        printf("\nFailed to open Boundary integration data writefile\n");
    }
    
    else{
        fprintf(outfile,"point\tobj"); // column headings
        for(i=1;i<numFunc;i++)
        {
            fprintf(outfile, "\tcnst%i",i);
        }
        fprintf(outfile, "\n");
        
        for(i=0;i<numLbound;i++)
        {
            fprintf(outfile,"%i",Lbound_nums[i]); // node num & length
            for(j=0;j<numFunc;j++)
            {
                p2 = j*numLbound + i; // point to correct place in Lbound
                fprintf(outfile,"\t%12.4e",Lbound[p2]);
            }
            fprintf(outfile, "\n");
        }
    }
    fclose(outfile);
    printf("\nBoundary integration data File written\n");
}

// function to output Vext & Grad in vtk format for Paraview
void CHakOutput::OutHJVTK(CHakMesh& m_mesh, double *Vnorm, double *Grad, int itt, char *datafile)
{
    // read data
    int NodeX = m_mesh.m_nodeX-2;
    int NodeY = m_mesh.m_nodeY-2;
    int NumNodes = m_mesh.m_numNodes;
    int **Nodes2 = m_mesh.m_pNodes2D;
    
    int i,k,n,m;
    
    FILE *outfile;	// File varible for output files
    char plotname[120];	// variable to change names of plotting output files
    
    // Write initial signed distance value information file
    sprintf(plotname,"%s_HJ_%i.vtk",datafile,itt); // set name for current output file
    outfile = fopen(plotname, "w");
    
    if(outfile == NULL){
        printf("\nFailed to open results paraview input writefile\n");
    }
    else
    {
        fprintf(outfile,"# vtk DataFile Version 3.0\n");
        fprintf(outfile,"Para0\n");
        fprintf(outfile,"ASCII\n");
        fprintf(outfile,"DATASET RECTILINEAR_GRID\n");
        fprintf(outfile,"DIMENSIONS %i %i %i\n", NodeX, NodeY, 1);
        fprintf(outfile,"X_COORDINATES %i int\n", NodeX);
        for(i=0;i<NodeX;i++){fprintf(outfile,"%i ", i);}
        fprintf(outfile,"\nY_COORDINATES %i int\n", NodeY);
        for(i=0;i<NodeY;i++){fprintf(outfile,"%i ",i);}
        fprintf(outfile,"\nZ_COORDINATES %i int\n", 1);
        for(i=0;i<1;i++){fprintf(outfile,"%i ",i);}
        
        NodeX++;
        NodeY++; // need to increase
        
        fprintf(outfile,"\n\nPOINT_DATA %i\n", NumNodes);
        
        // Now to print out nodal Vext
        fprintf(outfile,"SCALARS Vext float\n");
        fprintf(outfile,"LOOKUP_TABLE default\n");
        for(m=1;m<NodeY;m++)
        {
            for(n=1;n<NodeX;n++)
            {
                k = Nodes2[n][m]; // node number
                fprintf(outfile,"%12.4e\n",Vnorm[k]); // Vext value
            }
        }
        
        // Now to print out nodal grad
        fprintf(outfile,"SCALARS Grad float\n");
        fprintf(outfile,"LOOKUP_TABLE default\n");
        for(m=1;m<NodeY;m++)
        {
            for(n=1;n<NodeX;n++)
            {
                k = Nodes2[n][m]; // node number
                fprintf(outfile,"%12.4e\n",Grad[k]); // Grad value
            }
        }
    }
    fclose(outfile);
    printf("\nVelocity & Gradient File written (Paraview)\n");
}

// function to output displacements (and mode shapes) in vtk format for Paraview
void CHakOutput::OutDispVTK(CHakMesh& m_mesh, int numCase, double *disp, int num_eig, double *vec, int itt, char *datafile)
{
    // read data
    int NodeX = m_mesh.m_nodeX-2;
    int NodeY = m_mesh.m_nodeY-2;
    int NumNodes = m_mesh.m_numNodes;
    int **Nodes2 = m_mesh.m_pNodes2D;
    
    int i,k,n,m,c,ind;
    
    FILE *outfile;	// File varible for output files
    char plotname[120];	// variable to change names of plotting output files
    
    // Write initial signed distance value information file
    sprintf(plotname,"%s_Disp_%i.vtk",datafile,itt); // set name for current output file
    outfile = fopen(plotname, "w");
    
    if(outfile == NULL){
        printf("\nFailed to open results paraview input writefile\n");
    }
    else
    {
        fprintf(outfile,"# vtk DataFile Version 3.0\n");
        fprintf(outfile,"Para0\n");
        fprintf(outfile,"ASCII\n");
        fprintf(outfile,"DATASET RECTILINEAR_GRID\n");
        fprintf(outfile,"DIMENSIONS %i %i %i\n", NodeX, NodeY, 1);
        fprintf(outfile,"X_COORDINATES %i int\n", NodeX);
        for(i=0;i<NodeX;i++){fprintf(outfile,"%i ", i);}
        fprintf(outfile,"\nY_COORDINATES %i int\n", NodeY);
        for(i=0;i<NodeY;i++){fprintf(outfile,"%i ",i);}
        fprintf(outfile,"\nZ_COORDINATES %i int\n", 1);
        for(i=0;i<1;i++){fprintf(outfile,"%i ",i);}
        
        NodeX++;
        NodeY++; // need to incease
        
        fprintf(outfile,"\n\nPOINT_DATA %i\n", NumNodes);
        for(c=0;c<numCase;c++)
        {
            fprintf(outfile,"VECTORS disp%i float\n",c+1);
            
            ind = c*(NUM_DOF*NumNodes); // start of current disp vector
            // Now to print out nodal displacements (for this case)
            for(m=1;m<NodeY;m++)
            {
                for(n=1;n<NodeX;n++)
                {
                    k = (Nodes2[n][m]*NUM_DOF) + ind; // indicate correct place in disp
                    fprintf(outfile,"%12.4e %12.4e 0.0\n",disp[k],disp[k+1]); // disp vector
                }
            }
        }
        
        for(c=0;c<num_eig;c++)
        {
            fprintf(outfile,"VECTORS mode%i float\n",c+1);
            
            ind = c*(NUM_DOF*NumNodes); // start of current eigenvector
            // Now to print out nodal displacements (for this eigenvector)
            for(m=1;m<NodeY;m++)
            {
                for(n=1;n<NodeX;n++)
                {
                    k = (Nodes2[n][m]*NUM_DOF) + ind; // indicate correct place in vec
                    fprintf(outfile,"%12.4e %12.4e 0.0\n",vec[k],vec[k+1]); // eigenvector
                }
            }
        }
        
    }
    fclose(outfile);
    printf("\nDisplacement File written (Paraview)\n");
}

// function to output object & constraint convergence data
void CHakOutput::OutConv(int itt, CHakOptProblem *lsprob, double *Obj, double *constr, char *datafile)
{
    int numCon = lsprob->m_numConstraint;
    int i,j,ind; // incrementor
    FILE *outfile;	// File varible for output files
    char plotname[120];	// variable to change names of plotting output files
    
    sprintf(plotname,"%s_Convergence.txt",datafile); // set name for output file
    outfile = fopen(plotname, "w");
    if(outfile == NULL){
        printf("\nFailed to open Convergence writefile\n");
    }
    else{
        // column headings
        fprintf(outfile,"Iteration\t"); // Iteration
        
        // objective
        switch (lsprob->m_idObjType)
        {
            case 1:
                fprintf(outfile,"Compliance\t");
                break;
            case 2:
                fprintf(outfile,"Volume\t");
                break;
            case 3:
                fprintf(outfile,"1st Frequency\t");
                break;
            case 5:
                fprintf(outfile,"Mechanial Advantage\t");
                break;
            default:
                fprintf(outfile, "Objective\t");
        }
        
        // constraints
        for(i=0;i<numCon;i++)
        {
            switch (lsprob->m_pConstraint[i].type)
            {
                case 1:
                    fprintf(outfile,"Volume\t");
                    break;
                case 2:
                    fprintf(outfile,"Mass\t");
                    break;
                case 3:
                    fprintf(outfile,"Compliance\t");
                    break;
                case 5:
                    fprintf(outfile,"Displacement diff\t");
                    break;
                case 6:
                    fprintf(outfile,"Eigenvalue ratio\t");
                    break;
                case 7:
                    fprintf(outfile,"Input displacement\t");
                    break;
                case 10:
                    fprintf(outfile,"Support cost\t");
                    break;
                default:
                    fprintf(outfile, "Constraint %i\t",i+1);
            }
        }
        fprintf(outfile, "\n");
        
        // DATA
        for(i=0;i<=itt;i++)
        {
            fprintf(outfile,"%i\t",i); // Iteration
            fprintf(outfile,"%12.4e",Obj[i]); // Objective
            
            for(j=0;j<numCon;j++)
            {
                ind = (numCon*i)+j;
                fprintf(outfile,"\t%12.4e",constr[ind]); // Constraint
            }
            fprintf(outfile, "\n");
        }
    }
    fclose(outfile);
    printf("\nConvergence History file written");
}

// function to output covergence history of frequencies
void CHakOutput::OutFreq(int itt, int num_eig, double *freq, char *datafile)
{
    int i,j,ind; // incrementor
    int tot_eig = 2*num_eig;
    FILE *outfile;	// File varible for output files
    char plotname[120];	// variable to change names of plotting output files
    
    sprintf(plotname,"%s_Frequency.txt",datafile); // set name for output file
    outfile = fopen(plotname, "w");
    if(outfile == NULL){
        printf("\nFailed to open Frequency writefile\n");
    }
    else{
        // column headings
        fprintf(outfile,"Iteration");
        for(j=0;j<tot_eig;j++){ fprintf(outfile,"\tFreq %i",j+1); }
        fprintf(outfile, "\n");
        
        // DATA
        ind = 0;
        for(i=0;i<=itt;i++)
        {
            fprintf(outfile,"%i",i); // Iteration
            for(j=0;j<tot_eig;j++)
            {
                fprintf(outfile,"\t%12.4e",freq[ind++]); // frequency
            }
            fprintf(outfile, "\n");
        }
    }
    fclose(outfile);
    printf("\nFrequency History file written");
}

// function to output bar areas
void CHakOutput::OutBars(CHakMesh& m_mesh, int numFunc, double *sens, int pinfo, int itt, char *datafile)
{
    // read data
    int NumBars = m_mesh.m_numBars;
    int NumNodes = m_mesh.m_numNodes;
    Elem **Number = m_mesh.m_pNumber;
    Coord *nc = m_mesh.m_pNodeCoord;
    int elemX = m_mesh.m_elemX;
    int elemY = m_mesh.m_elemY;
    int xend = elemX-1;
    int yend = elemY-1;
    
    int i,j,n,m;
    
    FILE *outfile;	// File varible for output files
    char plotname[120];	// variable to change names of plotting output files
    
    // Write initial signed distance value information file
    sprintf(plotname,"%s_Bars_%i.vtk",datafile,itt); // set name for current output file
    outfile = fopen(plotname, "w");
    
    if(outfile == NULL){
        printf("\nFailed to open bars paraview input writefile\n");
    }
    else
    {
        fprintf(outfile,"# vtk DataFile Version 3.0\n");
        fprintf(outfile,"Para0\n");
        fprintf(outfile,"ASCII\n");
        fprintf(outfile,"DATASET UNSTRUCTURED_GRID\n");
        fprintf(outfile,"POINTS %i float\n",NumNodes);
        for(i=0;i<NumNodes;i++)
        {
            fprintf(outfile,"%lf %lf 0.0\n",nc[i].x,nc[i].y);
        }
        
        // define bar cells
        fprintf(outfile,"\n\nCELLS %i %i\n", NumBars, (3*NumBars));
        
        // loop through all elements
        for(j=0;j<elemY;j++)
        {
            for(i=0;i<elemX;i++)
            {
                // use bottom edge as x-bars
                fprintf(outfile,"2 %i %i\n", Number[i][j].a, Number[i][j].b);
                
                // also use top edge for top row of elements
                if(j == yend)
                {
                    fprintf(outfile,"2 %i %i\n", Number[i][j].d, Number[i][j].c);
                }
            }
        }
        
        // loop through all elements
        for(i=0;i<elemX;i++)
        {
            for(j=0;j<elemY;j++)
            {
                // use left edge as y-bars
                fprintf(outfile,"2 %i %i\n", Number[i][j].a, Number[i][j].d);
                
                // also use right edge for last column of elements
                if(i == xend)
                {
                    fprintf(outfile,"2 %i %i\n", Number[i][j].b, Number[i][j].c);
                }
            }
        }
        
        // define bar cell types (line)
        fprintf(outfile,"\n\nCELL_TYPES %i\n", NumBars);
        for(i=0;i<NumBars;i++){ fprintf(outfile,"3\n"); }
        
        // write out bar areas
        fprintf(outfile,"\nCELL_DATA %i\nSCALARS areas double 1\nLOOKUP_TABLE default\n",NumBars);
        for(i=0;i<NumBars;i++)
        {
            fprintf(outfile,"%lf\n",m_mesh.m_pBarAreas[i]);
        }
        
        // also output bar element sensitivities (if required)
        if(pinfo==3)
        {
            for(n=0;n<numFunc;n++)
            {
                fprintf(outfile,"\n\nSCALARS sens%i double 1\nLOOKUP_TABLE default\n",n);
                m=n*NumBars;
                for(i=0;i<NumBars;i++)
                {
                    fprintf(outfile,"%lf\n",sens[m++]);
                }
            }
        }
    }
    fclose(outfile);
    printf("\nBar Info File written (Paraview)\n");
}

// function to output designable bc varibles
void CHakOutput::OutDesBC(CHakMesh& m_mesh, double *sens, int pinfo, int itt, char *datafile)
{
    // read data
    int *BC_nums = m_mesh.m_pBcNums;
    double *K_bc = m_mesh.m_K_bc;
    int NodeX = m_mesh.m_nodeX-2;
    int NodeY = m_mesh.m_nodeY-2;
    int NumElem = m_mesh.m_numElem;
    
    int i,temp;
    
    FILE *outfile;	// File varible for output files
    char plotname[120];	// variable to change names of plotting output files
    
    sprintf(plotname,"%s_DesBC_%i.vtk",datafile,itt); // set name for current output file
    outfile = fopen(plotname, "w");
    
    if(outfile == NULL){
        printf("\nFailed to open designable bc paraview writefile\n");
    }
    else
    {
        fprintf(outfile,"# vtk DataFile Version 3.0\n");
        fprintf(outfile,"Para0\n");
        fprintf(outfile,"ASCII\n");
        fprintf(outfile,"DATASET RECTILINEAR_GRID\n");
        fprintf(outfile,"DIMENSIONS %i %i %i\n", NodeX, NodeY, 1);
        fprintf(outfile,"X_COORDINATES %i int\n", NodeX);
        for(i=0;i<NodeX;i++){fprintf(outfile,"%i ", i);}
        fprintf(outfile,"\nY_COORDINATES %i int\n", NodeY);
        for(i=0;i<NodeY;i++){fprintf(outfile,"%i ",i);}
        fprintf(outfile,"\nZ_COORDINATES %i int\n", 1);
        for(i=0;i<1;i++){fprintf(outfile,"%i ",i);}
        
        fprintf(outfile,"\n\nCELL_DATA %i\n", NumElem);
        fprintf(outfile,"SCALARS desBC double 1\n");
        fprintf(outfile,"LOOKUP_TABLE default\n");
        
        temp = 0;
        for(i=0;i<NumElem;i++)
        {
            if(i==BC_nums[temp]) {
                fprintf(outfile,"%lf\n", K_bc[temp++]);
            }
            else {
                fprintf(outfile,"0.0\n");
            }
        }
        
        // also output sensitivities (if required)
        if(pinfo==3)
        {
            fprintf(outfile,"\n\nSCALARS comp_sens double 1\nLOOKUP_TABLE default\n");
            temp = 0;
            for(i=0;i<NumElem;i++)
            {
                if(i==BC_nums[temp]) {
                    fprintf(outfile,"%lf\n", sens[temp++]);
                }
                else {
                    fprintf(outfile,"0.0\n");
                }
            }
        }
    }
    
    fclose(outfile);
    printf("\nDesignable BC variable File written (Paraview)\n");
}

// function to output designable material varibles
void CHakOutput::OutDesMat(CHakMesh& m_mesh, double *alpha, double aMin, int num_sens, double *sens, int pinfo, int itt, char *datafile)
{
    // read data
    int NumDesMat = m_mesh.m_numDesignableMaterialt;
    int *mat_elems = m_mesh.m_pMaterialElems;
    double *mat_vars = m_mesh.m_pMaterialVars;
    int NodeX = m_mesh.m_nodeX-2;
    int NodeY = m_mesh.m_nodeY-2;
    int NumElem = m_mesh.m_numElem;
    
    int i,j,temp;
    
    FILE *outfile;	// File varible for output files
    char plotname[120];	// variable to change names of plotting output files
    
    sprintf(plotname,"%s_DesMat_%i.vtk",datafile,itt); // set name for current output file
    outfile = fopen(plotname, "w");
    
    if(outfile == NULL){
        printf("\nFailed to open designable material paraview writefile\n");
    }
    else
    {
        fprintf(outfile,"# vtk DataFile Version 3.0\n");
        fprintf(outfile,"Para0\n");
        fprintf(outfile,"ASCII\n");
        fprintf(outfile,"DATASET RECTILINEAR_GRID\n");
        fprintf(outfile,"DIMENSIONS %i %i %i\n", NodeX, NodeY, 1);
        fprintf(outfile,"X_COORDINATES %i int\n", NodeX);
        for(i=0;i<NodeX;i++){fprintf(outfile,"%i ", i);}
        fprintf(outfile,"\nY_COORDINATES %i int\n", NodeY);
        for(i=0;i<NodeY;i++){fprintf(outfile,"%i ",i);}
        fprintf(outfile,"\nZ_COORDINATES %i int\n", 1);
        for(i=0;i<1;i++){fprintf(outfile,"%i ",i);}
        
        fprintf(outfile,"\n\nCELL_DATA %i\n", NumElem);
        fprintf(outfile,"SCALARS desMat double 1\n");
        fprintf(outfile,"LOOKUP_TABLE default\n");
        
        temp = 0;
        for(i=0;i<NumElem;i++)
        {
            if(i==mat_elems[temp] && alpha[i] > aMin) {
                fprintf(outfile,"%lf\n", mat_vars[temp]);
            }
            else {
                fprintf(outfile,"-1.0\n");
            }
            if(i==mat_elems[temp]){temp++;}
        }
        
        // also output sensitivities (if required)
        if(pinfo==3)
        {
            for(j=0;j<num_sens;j++)
            {
                temp = 0;
                fprintf(outfile,"\n\nSCALARS Sens%i double 1\nLOOKUP_TABLE default\n",j+1);
                for(i=0;i<NumElem;i++)
                {
                    if(i==mat_elems[temp]) {
                        fprintf(outfile,"%lf\n", sens[(j*NumDesMat)+temp]);
                        temp++;
                    }
                    else {
                        fprintf(outfile,"0.0\n");
                    }
                }
            }
        }
    }
    
    fclose(outfile);
    printf("\nDesignable Material variable File written (Paraview)\n");
}

// function to print solution to screen
void CHakOutput::report(int n,  int m, int nu, double **v)
{
    int i;
    int lim = (nu < 10) ? nu : 10;
    int lim2 = (n < 10) ? n : 10;
    
    printf("\nPrinting solution ...");
    printf("\n\nX values\n--------");
    for(i=0;i<lim2;i++)
    {
        printf("\nx[%i] = %12.6e",i,v[0][i]);
    }
    
    printf("\n\nZ values\n--------");
    for(i=0;i<lim2;i++)
    {
        printf("\nz[%i] = %12.6e",i,v[1][i]);
    }
    
    if(nu > 0)
    {
        printf("\n\nS values\n--------");
        for(i=0;i<lim;i++)
        {
            printf("\ns[%i] = %12.6e",i,v[2][i]);
        }
        
        printf("\n\nW values\n--------");
        for(i=0;i<lim;i++)
        {
            printf("\nw[%i] = %12.6e",i,v[3][i]);
        }
    }
    
    printf("\n\nY values\n--------");
    for(i=0;i<m;i++)
    {
        printf("\ny[%i] = %12.6e",i,v[4][i]);
    }
    
    printf("\n\n-------------------");
}


//
// Implementation - Interfaces
//
// knai20@bath.ac.uk: Non OOD functions

// function to write element and node numbering info
void CHakOutput::OutNumber(mesh *inMesh, char *datafile)
{
	// read data
	int m_elemX = inMesh->elemX;
	int elemY = inMesh->elemY;
	Elem **Number = inMesh->Number;
	int i,j; // incrementors
	FILE *outfile;	// File varible for output files
	char plotname[120];	// variable to change names of plotting output files

	sprintf(plotname,"%s_Numbering.txt",datafile); // set name for current output file
	outfile = fopen(plotname, "w");
	if(outfile == NULL){
		printf("\nFailed to open Numbering writefile\n");
	}
	else{
		fprintf(outfile,"n\ta\tb\tc\td\n"); // column headings
		for(j=0;j<elemY;j++)
		{
			for(i=0;i<m_elemX;i++)
			{
				fprintf(outfile,"%i\t",Number[i][j].n);
				fprintf(outfile,"%i\t",Number[i][j].a);
				fprintf(outfile,"%i\t",Number[i][j].b);
				fprintf(outfile,"%i\t",Number[i][j].c);
				fprintf(outfile,"%i\t",Number[i][j].d);
				fprintf(outfile,"\n");
			}
		}
	}
	fclose(outfile);
	printf("Numbering info file written\n");
}

// function to write node co-ordinate information file
void CHakOutput::OutNodeCoord(mesh *inMesh, char *datafile)
{
	// read data
	int NumNodes = inMesh->NumNodes;
	Coord *NodeCoord = inMesh->NodeCoord;

	int i; // incrementors
	FILE *outfile;	// File varible for output files
	char plotname[120];	// variable to change names of plotting output files

	sprintf(plotname,"%s_NCrd.txt",datafile); // set name for current output file
	outfile = fopen(plotname, "w");
	if(outfile == NULL){
		printf("\nFailed to open Node_Coord writefile\n");
	}
	else{
		for(i=0;i<NumNodes;i++)
		{
			fprintf(outfile,"%i\t",i);
			fprintf(outfile,"%lf\t",NodeCoord[i].x);
			fprintf(outfile,"%lf\t",NodeCoord[i].y);
			fprintf(outfile,"\n");
		}
	}
	fclose(outfile);
	printf("Node Co-ordinate File written");
}

// fucntion to output lsf (& maybe alpha) in vtk format
void CHakOutput::OutPLotShapeVTK2(mesh *inMesh, double *lsf, double *alpha, int pinfo, int itt, char *datafile)
{
	// read data
	int NodeX = inMesh->NodeX-2;
	int NodeY = inMesh->NodeY-2;
	int NumNodes = inMesh->NumNodes;
	int **Nodes2 = inMesh->Nodes2;

	int i,n,m,num;

	FILE *outfile;	// File varible for output files
	char plotname[120];	// variable to change names of plotting output files

	// Write initial signed distance value information file
	sprintf(plotname,"%s_Lsf_%i.vtk",datafile,itt); // set name for current output file
	outfile = fopen(plotname, "w");

	if(outfile == NULL){
		printf("\nFailed to open signed distance paraview writefile\n");
	}
	else
	{
		fprintf(outfile,"# vtk DataFile Version 3.0\n");
		fprintf(outfile,"Para0\n");
		fprintf(outfile,"ASCII\n");
		fprintf(outfile,"DATASET RECTILINEAR_GRID\n");
		fprintf(outfile,"DIMENSIONS %i %i %i\n", NodeX, NodeY, 1);
		fprintf(outfile,"X_COORDINATES %i int\n", NodeX);
		for(i=0;i<NodeX;i++){fprintf(outfile,"%i ", i);}
		fprintf(outfile,"\nY_COORDINATES %i int\n", NodeY);
		for(i=0;i<NodeY;i++){fprintf(outfile,"%i ",i);}
		fprintf(outfile,"\nZ_COORDINATES %i int\n", 1);
		for(i=0;i<1;i++){fprintf(outfile,"%i ",i);}

		fprintf(outfile,"\n\nPOINT_DATA %i\n", NumNodes);
		fprintf(outfile,"SCALARS lsf double 1\n");
		fprintf(outfile,"LOOKUP_TABLE default\n");

		// Now to print out the node values
		NodeX++;
		NodeY++;

        for(m=1;m<NodeY;m++)
        {
            for(n=1;n<NodeX;n++)
            {
				num = Nodes2[n][m];
				fprintf(outfile,"%lf\n", lsf[num]);
			}
		}

		// also output element area ratios (if required)
		if(pinfo==3)
		{
			int NumElem = inMesh->NumElem;
			fprintf(outfile,"\n\nCELL_DATA %i\n", NumElem);
			fprintf(outfile,"SCALARS alpha double 1\n");
			fprintf(outfile,"LOOKUP_TABLE default\n");

			for(n=0;n<NumElem;n++)
			{
				fprintf(outfile,"%lf\n", alpha[n]);
			}
		}
	}
	fclose(outfile);
	printf("\nSigned Distance info file written (Paraview)\n");
}

// function to output boundary as a mesh for Paraview (with shape sensitivities)
void CHakOutput::OutBoundVTK(mesh *inMesh, boundary *bound_in, int num_sens, double **Sens, int itt, char *datafile)
{
	// read data
	int NumNodes = inMesh->NumNodes;
	Coord *nc = inMesh->NodeCoord;
	Coord *AuxNodes = bound_in->AuxNodes;
	int NumBound = bound_in->NumBound;
	Bseg *Boundary = bound_in->Bound;
	int Ntot = NumNodes + bound_in->NumAux;

	int i,n;

	FILE *outfile;	// File varible for output files
	char plotname[120];	// variable to change names of plotting output files

	// Write initial signed distance value information file
	sprintf(plotname,"%s_Bound-Mesh_%i.vtk",datafile,itt); // set name for current output file
	outfile = fopen(plotname, "w");

	if(outfile == NULL){
		printf("\nFailed to open results paraview input writefile\n");
	}
	else
	{
		fprintf(outfile,"# vtk DataFile Version 3.0\n");
		fprintf(outfile,"Para0\n");
		fprintf(outfile,"ASCII\n");
		fprintf(outfile,"DATASET UNSTRUCTURED_GRID\n");
		fprintf(outfile,"POINTS %i float\n",Ntot);
		for(i=0;i<NumNodes;i++)
		{
			fprintf(outfile,"%lf %lf 0.0\n",nc[i].x,nc[i].y);
		}
		n = Ntot-NumNodes;
		for(i=0;i<n;i++)
		{
			fprintf(outfile,"%lf %lf 0.0\n",AuxNodes[i].x,AuxNodes[i].y);
		}

		fprintf(outfile,"\nCELLS %i %i\n",NumBound,(NumBound*3));
		for(i=0;i<NumBound;i++)
		{
			fprintf(outfile,"2 %i %i\n",Boundary[i].n1,Boundary[i].n2);
		}

		fprintf(outfile,"\nCELL_TYPES %i\n",NumBound);
		for(i=0;i<NumBound;i++)
		{
			fprintf(outfile,"3\n");
		}

		fprintf(outfile,"\nPOINT_DATA %i\n",Ntot);
		for(n=0;n<num_sens;n++)
		{
			fprintf(outfile,"SCALARS Sens%i float 1\n",n+1);
			fprintf(outfile,"LOOKUP_TABLE default\n");
			for(i=0;i<Ntot;i++)
			{
				fprintf(outfile,"%12.4e\n",Sens[n][i]);
			}
		}
	}

	fclose(outfile);
	printf("\nBoudnary Mesh Info File written (Paraview)\n");
}

// function to output boundary integration data
void CHakOutput::OutBoundInt(int numFunc, int numLbound, int *Lbound_nums, double *Lbound, int itt, char *datafile)
{
	int i,j,p2;
	FILE *outfile;	// File varible for output files
	char plotname[120];	// variable to change names of plotting output files

	sprintf(plotname,"%s_Bint_%i.txt",datafile,itt); // set name
	outfile = fopen(plotname, "w");
	if(outfile == NULL){
		printf("\nFailed to open Boundary integration data writefile\n");
	}

	else{
		fprintf(outfile,"point\tobj"); // column headings
		for(i=1;i<numFunc;i++)
		{
			fprintf(outfile, "\tcnst%i",i);
		}
		fprintf(outfile, "\n");

		for(i=0;i<numLbound;i++)
		{
			fprintf(outfile,"%i",Lbound_nums[i]); // node num & length
			for(j=0;j<numFunc;j++)
			{
				p2 = j*numLbound + i; // point to correct place in Lbound
				fprintf(outfile,"\t%12.4e",Lbound[p2]);
			}
			fprintf(outfile, "\n");
		}
	}
	fclose(outfile);
	printf("\nBoundary integration data File written\n");
}

// function to output Vext & Grad in vtk format for Paraview
void CHakOutput::OutHJVTK(mesh *inMesh, double *Vnorm, double *Grad, int itt, char *datafile)
{
	// read data
	int NodeX = inMesh->NodeX-2;
	int NodeY = inMesh->NodeY-2;
	int NumNodes = inMesh->NumNodes;
	int **Nodes2 = inMesh->Nodes2;

	int i,k,n,m;

	FILE *outfile;	// File varible for output files
	char plotname[120];	// variable to change names of plotting output files

	// Write initial signed distance value information file
	sprintf(plotname,"%s_HJ_%i.vtk",datafile,itt); // set name for current output file
	outfile = fopen(plotname, "w");

	if(outfile == NULL){
		printf("\nFailed to open results paraview input writefile\n");
	}
	else
	{
		fprintf(outfile,"# vtk DataFile Version 3.0\n");
		fprintf(outfile,"Para0\n");
		fprintf(outfile,"ASCII\n");
		fprintf(outfile,"DATASET RECTILINEAR_GRID\n");
		fprintf(outfile,"DIMENSIONS %i %i %i\n", NodeX, NodeY, 1);
		fprintf(outfile,"X_COORDINATES %i int\n", NodeX);
		for(i=0;i<NodeX;i++){fprintf(outfile,"%i ", i);}
		fprintf(outfile,"\nY_COORDINATES %i int\n", NodeY);
		for(i=0;i<NodeY;i++){fprintf(outfile,"%i ",i);}
		fprintf(outfile,"\nZ_COORDINATES %i int\n", 1);
		for(i=0;i<1;i++){fprintf(outfile,"%i ",i);}

		NodeX++;
		NodeY++; // need to increase

		fprintf(outfile,"\n\nPOINT_DATA %i\n", NumNodes);

		// Now to print out nodal Vext
		fprintf(outfile,"SCALARS Vext float\n");
		fprintf(outfile,"LOOKUP_TABLE default\n");
		for(m=1;m<NodeY;m++)
		{
			for(n=1;n<NodeX;n++)
			{
				k = Nodes2[n][m]; // node number
				fprintf(outfile,"%12.4e\n",Vnorm[k]); // Vext value
			}
		}

		// Now to print out nodal grad
		fprintf(outfile,"SCALARS Grad float\n");
		fprintf(outfile,"LOOKUP_TABLE default\n");
		for(m=1;m<NodeY;m++)
		{
			for(n=1;n<NodeX;n++)
			{
				k = Nodes2[n][m]; // node number
				fprintf(outfile,"%12.4e\n",Grad[k]); // Grad value
			}
		}
	}
	fclose(outfile);
	printf("\nVelocity & Gradient File written (Paraview)\n");
}

// function to output displacements (and mode shapes) in vtk format for Paraview
void CHakOutput::OutDispVTK(mesh *inMesh, int numCase, double *disp, int num_eig, double *vec, int itt, char *datafile)
{
	// read data
	int NodeX = inMesh->NodeX-2;
	int NodeY = inMesh->NodeY-2;
	int NumNodes = inMesh->NumNodes;
	int **Nodes2 = inMesh->Nodes2;

	int i,k,n,m,c,ind;

	FILE *outfile;	// File varible for output files
	char plotname[120];	// variable to change names of plotting output files

	// Write initial signed distance value information file
	sprintf(plotname,"%s_Disp_%i.vtk",datafile,itt); // set name for current output file
	outfile = fopen(plotname, "w");

	if(outfile == NULL){
		printf("\nFailed to open results paraview input writefile\n");
	}
	else
	{
		fprintf(outfile,"# vtk DataFile Version 3.0\n");
		fprintf(outfile,"Para0\n");
		fprintf(outfile,"ASCII\n");
		fprintf(outfile,"DATASET RECTILINEAR_GRID\n");
		fprintf(outfile,"DIMENSIONS %i %i %i\n", NodeX, NodeY, 1);
		fprintf(outfile,"X_COORDINATES %i int\n", NodeX);
		for(i=0;i<NodeX;i++){fprintf(outfile,"%i ", i);}
		fprintf(outfile,"\nY_COORDINATES %i int\n", NodeY);
		for(i=0;i<NodeY;i++){fprintf(outfile,"%i ",i);}
		fprintf(outfile,"\nZ_COORDINATES %i int\n", 1);
		for(i=0;i<1;i++){fprintf(outfile,"%i ",i);}

		NodeX++;
		NodeY++; // need to incease

		fprintf(outfile,"\n\nPOINT_DATA %i\n", NumNodes);
		for(c=0;c<numCase;c++)
		{
			fprintf(outfile,"VECTORS disp%i float\n",c+1);

			ind = c*(NUM_DOF*NumNodes); // start of current disp vector
			// Now to print out nodal displacements (for this case)
			for(m=1;m<NodeY;m++)
			{
				for(n=1;n<NodeX;n++)
				{
					k = (Nodes2[n][m]*NUM_DOF) + ind; // indicate correct place in disp
					fprintf(outfile,"%12.4e %12.4e 0.0\n",disp[k],disp[k+1]); // disp vector
				}
			}
		}

		for(c=0;c<num_eig;c++)
		{
			fprintf(outfile,"VECTORS mode%i float\n",c+1);

			ind = c*(NUM_DOF*NumNodes); // start of current eigenvector
			// Now to print out nodal displacements (for this eigenvector)
			for(m=1;m<NodeY;m++)
			{
				for(n=1;n<NodeX;n++)
				{
					k = (Nodes2[n][m]*NUM_DOF) + ind; // indicate correct place in vec
					fprintf(outfile,"%12.4e %12.4e 0.0\n",vec[k],vec[k+1]); // eigenvector
				}
			}
		}

	}
	fclose(outfile);
	printf("\nDisplacement File written (Paraview)\n");
}

// function to output object & constraint convergence data
void CHakOutput::OutConv(int itt, prob *lsprob, double *Obj, double *constr, char *datafile)
{
	int numCon = lsprob->num;
	int i,j,ind; // incrementor
	FILE *outfile;	// File varible for output files
	char plotname[120];	// variable to change names of plotting output files

	sprintf(plotname,"%s_Convergence.txt",datafile); // set name for output file
	outfile = fopen(plotname, "w");
	if(outfile == NULL){
		printf("\nFailed to open Convergence writefile\n");
	}
	else{
		// column headings
		fprintf(outfile,"Iteration\t"); // Iteration

		// objective
		switch (lsprob->obj)
		{
			case 1:
				fprintf(outfile,"Compliance\t");
				break;
			case 2:
				fprintf(outfile,"Volume\t");
				break;
			case 3:
				fprintf(outfile,"1st Frequency\t");
				break;
			case 5:
				fprintf(outfile,"Mechanial Advantage\t");
				break;
			default:
				fprintf(outfile, "Objective\t");
		}

		// constraints
		for(i=0;i<numCon;i++)
		{
			switch (lsprob->con[i].type)
			{
				case 1:
					fprintf(outfile,"Volume\t");
					break;
				case 2:
					fprintf(outfile,"Mass\t");
					break;
				case 3:
					fprintf(outfile,"Compliance\t");
					break;
				case 5:
					fprintf(outfile,"Displacement diff\t");
					break;
				case 6:
					fprintf(outfile,"Eigenvalue ratio\t");
					break;
				case 7:
					fprintf(outfile,"Input displacement\t");
					break;
				case 10:
					fprintf(outfile,"Support cost\t");
					break;
				default:
					fprintf(outfile, "Constraint %i\t",i+1);
			}
		}
		fprintf(outfile, "\n");

		// DATA
		for(i=0;i<=itt;i++)
		{
			fprintf(outfile,"%i\t",i); // Iteration
			fprintf(outfile,"%12.4e",Obj[i]); // Objective

			for(j=0;j<numCon;j++)
			{
				ind = (numCon*i)+j;
				fprintf(outfile,"\t%12.4e",constr[ind]); // Constraint
			}
			fprintf(outfile, "\n");
		}
	}
	fclose(outfile);
	printf("\nConvergence History file written");
}

// function to output covergence history of frequencies
void CHakOutput::OutFreq(int itt, int num_eig, double *freq, char *datafile)
{
	int i,j,ind; // incrementor
	int tot_eig = 2*num_eig;
	FILE *outfile;	// File varible for output files
	char plotname[120];	// variable to change names of plotting output files

	sprintf(plotname,"%s_Frequency.txt",datafile); // set name for output file
	outfile = fopen(plotname, "w");
	if(outfile == NULL){
		printf("\nFailed to open Frequency writefile\n");
	}
	else{
		// column headings
		fprintf(outfile,"Iteration");
		for(j=0;j<tot_eig;j++){ fprintf(outfile,"\tFreq %i",j+1); }
		fprintf(outfile, "\n");

		// DATA
		ind = 0;
		for(i=0;i<=itt;i++)
		{
			fprintf(outfile,"%i",i); // Iteration
			for(j=0;j<tot_eig;j++)
			{
				fprintf(outfile,"\t%12.4e",freq[ind++]); // frequency
			}
			fprintf(outfile, "\n");
		}
	}
	fclose(outfile);
	printf("\nFrequency History file written");
}

// function to output bar areas
void CHakOutput::OutBars(mesh *inMesh, int numFunc, double *sens, int pinfo, int itt, char *datafile)
{
	// read data
	int NumBars = inMesh->NumBars;
	int NumNodes = inMesh->NumNodes;
	Elem **Number = inMesh->Number;
	Coord *nc = inMesh->NodeCoord;
	int m_elemX = inMesh->elemX;
	int elemY = inMesh->elemY;
	int xend = m_elemX-1;
	int yend = elemY-1;

	int i,j,n,m;

	FILE *outfile;	// File varible for output files
	char plotname[120];	// variable to change names of plotting output files

	// Write initial signed distance value information file
	sprintf(plotname,"%s_Bars_%i.vtk",datafile,itt); // set name for current output file
	outfile = fopen(plotname, "w");

	if(outfile == NULL){
		printf("\nFailed to open bars paraview input writefile\n");
	}
	else
	{
		fprintf(outfile,"# vtk DataFile Version 3.0\n");
		fprintf(outfile,"Para0\n");
		fprintf(outfile,"ASCII\n");
		fprintf(outfile,"DATASET UNSTRUCTURED_GRID\n");
		fprintf(outfile,"POINTS %i float\n",NumNodes);
		for(i=0;i<NumNodes;i++)
		{
			fprintf(outfile,"%lf %lf 0.0\n",nc[i].x,nc[i].y);
		}

		// define bar cells
		fprintf(outfile,"\n\nCELLS %i %i\n", NumBars, (3*NumBars));

		// loop through all elements
		for(j=0;j<elemY;j++)
		{
			for(i=0;i<m_elemX;i++)
			{
				// use bottom edge as x-bars
				fprintf(outfile,"2 %i %i\n", Number[i][j].a, Number[i][j].b);

				// also use top edge for top row of elements
				if(j == yend)
				{
					fprintf(outfile,"2 %i %i\n", Number[i][j].d, Number[i][j].c);
				}
			}
		}

		// loop through all elements
		for(i=0;i<m_elemX;i++)
		{
			for(j=0;j<elemY;j++)
			{
				// use left edge as y-bars
				fprintf(outfile,"2 %i %i\n", Number[i][j].a, Number[i][j].d);

				// also use right edge for last column of elements
				if(i == xend)
				{
					fprintf(outfile,"2 %i %i\n", Number[i][j].b, Number[i][j].c);
				}
			}
		}

		// define bar cell types (line)
		fprintf(outfile,"\n\nCELL_TYPES %i\n", NumBars);
		for(i=0;i<NumBars;i++){ fprintf(outfile,"3\n"); }

		// write out bar areas
		fprintf(outfile,"\nCELL_DATA %i\nSCALARS areas double 1\nLOOKUP_TABLE default\n",NumBars);
		for(i=0;i<NumBars;i++)
		{
			fprintf(outfile,"%lf\n",inMesh->bar_areas[i]);
		}

		// also output bar element sensitivities (if required)
		if(pinfo==3)
		{
			for(n=0;n<numFunc;n++)
			{
				fprintf(outfile,"\n\nSCALARS sens%i double 1\nLOOKUP_TABLE default\n",n);
				m=n*NumBars;
				for(i=0;i<NumBars;i++)
				{
					fprintf(outfile,"%lf\n",sens[m++]);
				}
			}
		}
	}
	fclose(outfile);
	printf("\nBar Info File written (Paraview)\n");
}

// function to output designable bc varibles
void CHakOutput::OutDesBC(mesh *inMesh, double *sens, int pinfo, int itt, char *datafile)
{
	// read data
	int *BC_nums = inMesh->BC_nums;
	double *K_bc = inMesh->K_bc;
	int NodeX = inMesh->NodeX-2;
	int NodeY = inMesh->NodeY-2;
	int NumElem = inMesh->NumElem;

	int i,temp;

	FILE *outfile;	// File varible for output files
	char plotname[120];	// variable to change names of plotting output files

	sprintf(plotname,"%s_DesBC_%i.vtk",datafile,itt); // set name for current output file
	outfile = fopen(plotname, "w");

	if(outfile == NULL){
		printf("\nFailed to open designable bc paraview writefile\n");
	}
	else
	{
		fprintf(outfile,"# vtk DataFile Version 3.0\n");
		fprintf(outfile,"Para0\n");
		fprintf(outfile,"ASCII\n");
		fprintf(outfile,"DATASET RECTILINEAR_GRID\n");
		fprintf(outfile,"DIMENSIONS %i %i %i\n", NodeX, NodeY, 1);
		fprintf(outfile,"X_COORDINATES %i int\n", NodeX);
		for(i=0;i<NodeX;i++){fprintf(outfile,"%i ", i);}
		fprintf(outfile,"\nY_COORDINATES %i int\n", NodeY);
		for(i=0;i<NodeY;i++){fprintf(outfile,"%i ",i);}
		fprintf(outfile,"\nZ_COORDINATES %i int\n", 1);
		for(i=0;i<1;i++){fprintf(outfile,"%i ",i);}

		fprintf(outfile,"\n\nCELL_DATA %i\n", NumElem);
		fprintf(outfile,"SCALARS desBC double 1\n");
		fprintf(outfile,"LOOKUP_TABLE default\n");

		temp = 0;
		for(i=0;i<NumElem;i++)
		{
			if(i==BC_nums[temp]) {
				fprintf(outfile,"%lf\n", K_bc[temp++]);
			}
			else {
				fprintf(outfile,"0.0\n");
			}
		}

		// also output sensitivities (if required)
		if(pinfo==3)
		{
			fprintf(outfile,"\n\nSCALARS comp_sens double 1\nLOOKUP_TABLE default\n");
			temp = 0;
			for(i=0;i<NumElem;i++)
			{
				if(i==BC_nums[temp]) {
					fprintf(outfile,"%lf\n", sens[temp++]);
				}
				else {
					fprintf(outfile,"0.0\n");
				}
			}
		}
	}

	fclose(outfile);
	printf("\nDesignable BC variable File written (Paraview)\n");
}

// function to output designable material varibles
void CHakOutput::OutDesMat(mesh *inMesh, double *alpha, double aMin, int num_sens, double *sens, int pinfo, int itt, char *datafile)
{
	// read data
    int NumDesMat = inMesh->NumDesMat;
	int *mat_elems = inMesh->mat_elems;
	double *mat_vars = inMesh->mat_vars;
	int NodeX = inMesh->NodeX-2;
	int NodeY = inMesh->NodeY-2;
	int NumElem = inMesh->NumElem;

	int i,j,temp;

	FILE *outfile;	// File varible for output files
	char plotname[120];	// variable to change names of plotting output files

	sprintf(plotname,"%s_DesMat_%i.vtk",datafile,itt); // set name for current output file
	outfile = fopen(plotname, "w");

	if(outfile == NULL){
		printf("\nFailed to open designable material paraview writefile\n");
	}
	else
	{
		fprintf(outfile,"# vtk DataFile Version 3.0\n");
		fprintf(outfile,"Para0\n");
		fprintf(outfile,"ASCII\n");
		fprintf(outfile,"DATASET RECTILINEAR_GRID\n");
		fprintf(outfile,"DIMENSIONS %i %i %i\n", NodeX, NodeY, 1);
		fprintf(outfile,"X_COORDINATES %i int\n", NodeX);
		for(i=0;i<NodeX;i++){fprintf(outfile,"%i ", i);}
		fprintf(outfile,"\nY_COORDINATES %i int\n", NodeY);
		for(i=0;i<NodeY;i++){fprintf(outfile,"%i ",i);}
		fprintf(outfile,"\nZ_COORDINATES %i int\n", 1);
		for(i=0;i<1;i++){fprintf(outfile,"%i ",i);}

		fprintf(outfile,"\n\nCELL_DATA %i\n", NumElem);
		fprintf(outfile,"SCALARS desMat double 1\n");
		fprintf(outfile,"LOOKUP_TABLE default\n");

		temp = 0;
		for(i=0;i<NumElem;i++)
		{
			if(i==mat_elems[temp] && alpha[i] > aMin) {
				fprintf(outfile,"%lf\n", mat_vars[temp]);
			}
			else {
				fprintf(outfile,"-1.0\n");
			}
            if(i==mat_elems[temp]){temp++;}
		}

		// also output sensitivities (if required)
		if(pinfo==3)
		{
            for(j=0;j<num_sens;j++)
            {
                temp = 0;
                fprintf(outfile,"\n\nSCALARS Sens%i double 1\nLOOKUP_TABLE default\n",j+1);
                for(i=0;i<NumElem;i++)
                {
                    if(i==mat_elems[temp]) {
                        fprintf(outfile,"%lf\n", sens[(j*NumDesMat)+temp]);
                        temp++;
                    }
                    else {
                        fprintf(outfile,"0.0\n");
                    }
                }
            }
		}
	}

	fclose(outfile);
	printf("\nDesignable Material variable File written (Paraview)\n");
}

// function to print solution to screen
void CHakOutput::report(int n,  int m, int nu, double **v)
{
	int i;
	int lim = (nu < 10) ? nu : 10;
	int lim2 = (n < 10) ? n : 10;

	printf("\nPrinting solution ...");
	printf("\n\nX values\n--------");
	for(i=0;i<lim2;i++)
	{
		printf("\nx[%i] = %12.6e",i,v[0][i]);
	}

	printf("\n\nZ values\n--------");
	for(i=0;i<lim2;i++)
	{
		printf("\nz[%i] = %12.6e",i,v[1][i]);
	}

	if(nu > 0)
	{
		printf("\n\nS values\n--------");
		for(i=0;i<lim;i++)
		{
			printf("\ns[%i] = %12.6e",i,v[2][i]);
		}

		printf("\n\nW values\n--------");
		for(i=0;i<lim;i++)
		{
			printf("\nw[%i] = %12.6e",i,v[3][i]);
		}
	}

	printf("\n\nY values\n--------");
	for(i=0;i<m;i++)
	{
		printf("\ny[%i] = %12.6e",i,v[4][i]);
	}

	printf("\n\n-------------------");
}



//Plot hole information
void CHakOutput::OutPLotHoleVTK(mesh *inMesh, double *h_Nsens, int *h_EmapX, int *h_EmapY, double *h_Esens, int *h_index, int h_count, double *h_lsf, double *h_area, int num_sens, int pinfo, int itt, char *datafile)
{
    // read data
    int NodeX = inMesh->NodeX-2;
    int NodeY = inMesh->NodeY-2;
    int NumNodes = inMesh->NumNodes;
    int **Nodes2 = inMesh->Nodes2;
    int elemX = inMesh->elemX;
    int elemY = inMesh->elemY;
    Elem **Number = inMesh->Number;

    int i,j,k,n,m,num;

    FILE *outfile;	// File varible for output files
    char plotname[120];	// variable to change names of plotting output files

    // Write initial signed distance value information file
    sprintf(plotname,"%s_hole_%i.vtk",datafile,itt); // set name for current output file
    outfile = fopen(plotname, "w");

    if(outfile == NULL){
        printf("\nFailed to open signed distance paraview writefile\n");
    }
    else
    {
        fprintf(outfile,"# vtk DataFile Version 3.0\n");
        fprintf(outfile,"Para0\n");
        fprintf(outfile,"ASCII\n");
        fprintf(outfile,"DATASET RECTILINEAR_GRID\n");
        fprintf(outfile,"DIMENSIONS %i %i %i\n", NodeX, NodeY, 1);
        fprintf(outfile,"X_COORDINATES %i int\n", NodeX);
        for(i=0;i<NodeX;i++){fprintf(outfile,"%i ", i);}
        fprintf(outfile,"\nY_COORDINATES %i int\n", NodeY);
        for(i=0;i<NodeY;i++){fprintf(outfile,"%i ",i);}
        fprintf(outfile,"\nZ_COORDINATES %i int\n", 1);
        for(i=0;i<1;i++){fprintf(outfile,"%i ",i);}

        fprintf(outfile,"\n\nPOINT_DATA %i\n", NumNodes);
        fprintf(outfile,"SCALARS index double 1\n");
        fprintf(outfile,"LOOKUP_TABLE default\n");

        // Now to print out the node values
        NodeX++;
        NodeY++;

        for(m=1;m<NodeY;m++)
        {
            for(n=1;n<NodeX;n++)
            {
                num = Nodes2[n][m];
                fprintf(outfile,"%i\n", h_index[num]);
            }
        }

        for(i=0;i<num_sens;i++)
        {
            fprintf(outfile,"SCALARS hole_sens%i double 1\n", i);
            fprintf(outfile,"LOOKUP_TABLE default\n");

            for(m=1;m<NodeY;m++)
            {
                for(n=1;n<NodeX;n++)
                {
                    num = Nodes2[n][m];
                    j = NumNodes*i;
                    fprintf(outfile,"%lf\n", h_Nsens[j+num]);
                }
            }
        }
        fprintf(outfile,"SCALARS hole_lsf double 1\n");
        fprintf(outfile,"LOOKUP_TABLE default\n");

        for(m=1;m<NodeY;m++)
        {
            for(n=1;n<NodeX;n++)
            {
                num = Nodes2[n][m];
                fprintf(outfile,"%lf\n", h_lsf[num]);
            }
        }

        int NumElem = inMesh->NumElem;
        fprintf(outfile,"\n\nCELL_DATA %i\n", NumElem);
        fprintf(outfile,"SCALARS alpha double 1\n");
        fprintf(outfile,"LOOKUP_TABLE default\n");

        for(n=0;n<NumElem;n++)
        {
            fprintf(outfile,"%lf\n", h_area[n]);
        }

    }
    fclose(outfile);
    printf("\nHole sensivtiy info file written (Paraview)\n");

    sprintf(plotname,"%s_hole_%i.txt",datafile,itt); // set name for current output file
    outfile = fopen(plotname, "w");
    /*for(k=0;k<num_sens;k++)
    {
        n = NumNodes*k;
        fprintf(outfile,"%f", h_Nsens[n+Number[0][0].a]);
        for(i=0;i<elemX;i++)
        {
            fprintf(outfile,"\t%f", Ihole_sens[n+Number[i][0].b]);
        }
        for(j=0;j<elemY;j++)
        {
            fprintf(outfile,"\n%f", Ihole_sens[n+Number[0][j].d]);
            for(i=0;i<elemX;i++)
            {
                fprintf(outfile,"\t%f", Ihole_sens[n+Number[i][j].c]);
            }
        }
        fprintf(outfile,"\n\n\n");
    }*/

    fprintf(outfile,"%f", h_lsf[Number[0][0].a]);
    for(i=0;i<elemX;i++)
    {
        fprintf(outfile,"\t%f", h_lsf[Number[i][0].b]);
    }
    for(j=0;j<elemY;j++)
    {
        fprintf(outfile,"\n%f", h_lsf[Number[0][j].d]);
        for(i=0;i<elemX;i++)
        {
            fprintf(outfile,"\t%f", h_lsf[Number[i][j].c]);
        }
    }

    fclose(outfile);/**/

}

void CHakOutput::OutPLotStresssVTK2(mesh *inMesh, double *lsf, double *EStress, int pinfo, int itt, char *datafile)
{
    // read data
    int NodeX = inMesh->NodeX-2;
    int NodeY = inMesh->NodeY-2;
    int NumNodes = inMesh->NumNodes;
    int **Nodes2 = inMesh->Nodes2;

    int i,n,m,num;

    FILE *outfile;	// File varible for output files
    char plotname[120];	// variable to change names of plotting output files

    // Write initial signed distance value information file
    sprintf(plotname,"A_%s_Stress_%i.vtk",datafile,itt); // set name for current output file
    outfile = fopen(plotname, "w");

    if(outfile == NULL){
        printf("\nFailed to open results paraview input writefile\n");
    }
    else
    {
        fprintf(outfile,"# vtk DataFile Version 3.0\n");
        fprintf(outfile,"Para0\n");
        fprintf(outfile,"ASCII\n");
        fprintf(outfile,"DATASET RECTILINEAR_GRID\n");
        fprintf(outfile,"DIMENSIONS %i %i %i\n", NodeX, NodeY, 1);
        fprintf(outfile,"X_COORDINATES %i int\n", NodeX);
        for(i=0;i<NodeX;i++){fprintf(outfile,"%i ", i);}
        fprintf(outfile,"\nY_COORDINATES %i int\n", NodeY);
        for(i=0;i<NodeY;i++){fprintf(outfile,"%i ",i);}
        fprintf(outfile,"\nZ_COORDINATES %i int\n", 1);
        for(i=0;i<1;i++){fprintf(outfile,"%i ",i);}

        fprintf(outfile,"\n\nPOINT_DATA %i\n", NumNodes);
        fprintf(outfile,"SCALARS lsf double 1\n");
        fprintf(outfile,"LOOKUP_TABLE default\n");

        // Now to print out the node values
        NodeX++;
        NodeY++;

        for(m=1;m<NodeY;m++)
        {
            for(n=1;n<NodeX;n++)
            {
                num = Nodes2[n][m];
                fprintf(outfile,"%f\n", lsf[num]);
            }
        }

        // also output element area ratios (if required)
        if(pinfo<5)
        {
            int NumElem = inMesh->NumElem;
            fprintf(outfile,"\n\nCELL_DATA %i\n", NumElem);
            fprintf(outfile,"SCALARS Element_Stress double 1\n");
            fprintf(outfile,"LOOKUP_TABLE default\n");

            for(n=0;n<NumElem;n++)
            {
                fprintf(outfile,"%f\n", EStress[n]);
            }
        }
    }
    fclose(outfile);
    printf("\nSigned Distance Info File written (Paraview)\n");
}

// function to output object & constraint convergence data for stress function
void OutConv(int itt, prob *lsprob, double *Obj, double *constr, char *datafile, double *MaxStress)
{
	int numCon = lsprob->num;
	int i,j,ind; // incrementor
	FILE *outfile;	// File varible for output files
	char plotname[120];	// variable to change names of plotting output files

	sprintf(plotname,"%s_Convergence.txt",datafile); // set name for output file
	outfile = fopen(plotname, "w");
	if(outfile == NULL){
		printf("\nFailed to open Convergence writefile\n");
	}
	else{
		// column headings
		fprintf(outfile,"Iteration\t"); // Iteration

		// objective
		switch (lsprob->obj)
		{
			case 1:
				fprintf(outfile,"Compliance\t");
				break;
            case 27:
                fprintf(outfile,"Pnorm_Stress\t");
                break;
			case 2:
				fprintf(outfile,"Volume\t");
				break;
			case 3:
				fprintf(outfile,"1st Frequency\t");
				break;
			case 5:
				fprintf(outfile,"Mechanial Advantage\t");
				break;
			default:
				fprintf(outfile, "Objective\t");
		}

		// constraints
		for(i=0;i<numCon;i++)
		{
			switch (lsprob->con[i].type)
			{
				case 1:
					fprintf(outfile,"Volume\t");
					break;
				case 2:
					fprintf(outfile,"Mass\t");
					break;
				case 3:
					fprintf(outfile,"Compliance\t");
					break;
				case 5:
					fprintf(outfile,"Displacement diff\t");
					break;
				case 6:
					fprintf(outfile,"Eigenvalue ratio\t");
					break;
				case 7:
					fprintf(outfile,"Input displacement\t");
					break;
				case 10:
					fprintf(outfile,"Support cost\t");
					break;
                case 27:
                    fprintf(outfile,"Pnorm_Stress\t");
                    break;
				default:
					fprintf(outfile, "Constraint %i\t",i+1);
			}
		}
        if(MaxStress[1]>0){fprintf(outfile, "Max_Stress\t");}
		fprintf(outfile, "\n");

		// DATA
		for(i=0;i<=itt;i++)
		{
			fprintf(outfile,"%i\t",i); // Iteration
			fprintf(outfile,"%12.4e",Obj[i]); // Objective

			for(j=0;j<numCon;j++)
			{
				ind = (numCon*i)+j;
				fprintf(outfile,"\t%12.4e",constr[ind]); // Constraint
			}
            if(MaxStress[i]>0){fprintf(outfile, "\t%12.4e",MaxStress[i]);}
			fprintf(outfile, "\n");
		}
	}
	fclose(outfile);
	printf("\nConvergence History file written");
}

