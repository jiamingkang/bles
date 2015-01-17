/*
	CHakOutput.h

  	Created on: Nov 24, 2014
	Author: Peter Dunning, JeeHang Lee

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

#ifndef CHakOutput_H_
#define CHakOutput_H_

// include
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "CommonTypes.h"

class CHakOutput {
public:
	CHakOutput();
	virtual ~CHakOutput();


// knai20@bath.ac.uk: OOD version
    
public:
    // function to write element and node numbering info
    void OutNumber(CHakMesh& m_mesh, char *datafile);
    
    // function to write node co-ordinate information file
    void OutNodeCoord(CHakMesh& m_mesh, char *datafile);
    
    // function to output lsf (& maybe alpha) in vtk format
    void OutPLotShapeVTK2(CHakMesh& m_mesh, double *lsf, double *alpha, int pinfo, int itt, char *datafile);
    
    // function to output boundary as a mesh for Paraview (with shape sensitivities)
    void OutBoundVTK(CHakMesh& m_mesh, boundary *bound_in, int num_sens, double **Sens, int itt, char *datafile);
    
    // function to output boundary integration data
    void OutBoundInt(int numFunc, int numLbound, int *Lbound_nums, double *Lbound, int itt, char *datafile);
    
    // function to output Vext & Grad in vtk format for Paraview
    void OutHJVTK(CHakMesh& m_mesh, double *Vnorm, double *Grad, int itt, char *datafile);
    
    // function to output displacements (and mode shapes) in vtk format for Paraview
    void OutDispVTK(CHakMesh& m_mesh, int numCase, double *disp, int num_eig, double *vec, int itt, char *datafile);
    
    // function to output object & constraint convergence data
    void OutConv(int itt, prob *lsprob, double *Obj, double *constr, char *datafile);
    
    // function to output covergence history of frequencies
    void OutFreq(int itt, int num_eig, double *freq, char *datafile);
    
    // function to output bar areas
    void OutBars(CHakMesh& m_mesh, int numFunc, double *sens, int pinfo, int itt, char *datafile);
    
    // function to output designable bc varibles
    void OutDesBC(CHakMesh& m_mesh, double *sens, int pinfo, int itt, char *datafile);
    
    // function to output designable material varibles
    void OutDesMat(CHakMesh& m_mesh, double *alpha, double aMin, int num_sens, double *sens, int pinfo, int itt, char *datafile);
    
    // function to print solution to screen
    void report(int n,  int m, int nu, double **v);

 // knai20@bath.ac.uk: Non OOD version

public:
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

	// function to print solution to screen
	void report(int n,  int m, int nu, double **v);
};

#endif /* CHakOutput_H_ */
