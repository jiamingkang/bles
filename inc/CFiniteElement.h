/*
   CFiniteElement.h

    Created on: 24 Nov 2014
    Author: Peter Dunning, Khalid Ismail

	-- This file is part of Topology Optimisation Opensource Project,
 	owned by MSO (Multidisciplinary and Structural Optimisation) Research Group
 	(http://people.bath.ac.uk/ens/MSORG/index.html) at University of Bath.

 	The project is led by Dr. Alicia Kim
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

#ifndef CFINITEELEMENT_H_
#define CFINITEELEMENT_H_


class CFiniteElement {
public:
	CFiniteElement();
	virtual ~CFiniteElement();

public:
	// Assembles global stiffness (& maybe mass) matrix for AFG method in triplet format (for MA57 solver)
	void AFG_Matrix(int mass, double **KE, double **ME, sp_mat *Kg, sp_mat *Mg, sp_mat *lump_mass, double *alpha,
						mesh *inMesh, isoMat *inMat, double aMin, double mMin);

	// function to compute a cut element area
	double InArea(int eStat, int eNum, int *Lnodes, short *NodeStat, int NumNodes, Coord *NodeCoord,
				  Coord *AuxNodes, int NumBound, Bseg *Boundary);

	// For odd (i + j = odd) Matrix enrites
	double Aodd(int i,int j,double e,double g);

	// For even (i + j = even) Matrix enrites
	double Aeven(int i,int j,double e,double g);

	// Complete square IN element stiffness matix computation
	// only for isotropic material
	void KEMatrix(double *KE, isoMat *mat, double t);

	// consistent mass matrix for a 2D plane element
	void MEMatrix(double *ME, double rho, double area, double t);

	// Function that assemebles the element matricies into the global matrix (in triplet form)
	// Symmetric matrix - upper triangle
	void Assemble2(int *K_begin, int *M_begin, int nNod, int nDof, int *tnodes, double *KE, double *ME,
						sp_mat *Kg, sp_mat *Mg, int mass);

	// function to free memory for a sparse matrix
	void free_sp_mat(sp_mat *m);

	// function to create memory for a sparse matrix
	void set_sp_mat(sp_mat *m);

	// function to remove dof from a sparse matrix
	void rem_sp_mat(sp_mat *m, int *map, int inc);

};

#endif /* CFINITEELEMENT_H_ */
