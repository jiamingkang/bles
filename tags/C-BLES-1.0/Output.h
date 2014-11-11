/*
 *  Output.h
 *
 *  Created by Peter Dunning
 *	Version 5.4 @ 23/04/2014
 */

#include "ls_types.h"

// function to write element and node numbering info
void OutNumber(mesh *inMesh, char *datafile);

// function to write node co-ordinate information file
void OutNodeCoord(mesh *inMesh, char *datafile);

// function to output lsf (& maybe alpha) in vtk format
void OutPLotShapeVTK2(mesh *inMesh, double *lsf, double *alpha, int pinfo, int itt, char *datafile);

// function to output boundary as a mesh for Paraview (with shape sensitivities)
void OutBoundVTK(mesh *inMesh, boundary *bound_in, int num_sens, double **Sens, int itt, char *datafile);

// function to output boundary integration data
void OutBoundInt(int numFunc, int numLbound, int *Lbound_nums, double *Lbound, int itt, char *datafile);

// function to output Vext & Grad in vtk format for Paraview
void OutHJVTK(mesh *inMesh, double *Vnorm, double *Grad, int itt, char *datafile);

// function to output displacements (and mode shapes) in vtk format for Paraview
void OutDispVTK(mesh *inMesh, int numCase, double *disp, int num_eig, double *vec, int itt, char *datafile);

// function to output object & constraint convergence data
void OutConv(int itt, prob *lsprob, double *Obj, double *constr, char *datafile);

// function to output covergence history of frequencies
void OutFreq(int itt, int num_eig, double *freq, char *datafile);

// function to output bar areas
void OutBars(mesh *inMesh, int numFunc, double *sens, int pinfo, int itt, char *datafile);

// function to output designable bc varibles
void OutDesBC(mesh *inMesh, double *sens, int pinfo, int itt, char *datafile);

// function to output designable material varibles
void OutDesMat(mesh *inMesh, double *alpha, double aMin, int num_sens, double *sens, int pinfo, int itt, char *datafile);
