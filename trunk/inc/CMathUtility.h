/*
	CMathUtility.h

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

#ifndef CMATHUTILITY_H_
#define CMATHUTILITY_H_

//
// Includes
//


//
// Pre-definitions
//
#define TINY 1.0E-5
#define SMALL 1.0E-8
#define BIG 1.0E6

class CMathUtility {
public:
	CMathUtility();
	virtual ~CMathUtility();

public:
	// function to compute gauss point coords
	void Gauss_Coord(mesh *inMesh, Coord *gCoord);

	// Function that caculates the area of any Polygon
	double PolyArea(int N, Coord *point);

	// Function to determine if two lines (p1 -> p2 & p3 -> p4) cross
	short LineCross(Coord *pts);

	// compute determinant of a 2x2 matrix
	double det2(double a, double b, double c, double d)
	{
		// compute determinant of a 2x2 matrix
		return( (a * d) - (b * c) );
	}

	// sub-fucntion for sorting in ascending order
	int dcmpfunc(const void * p1, const void * p2);

	// function to handle divide by zero
	double divZero(double num, double denom);

	// function to find an initial point
	int initialize(int n,  int m, int nu, double **v, double *A, double *b, double *u, double *c);

	int initialize2(int n,  int m, int nu, double **v, double *A, double *b, double *u, double *c);

	// function to compute del(v) - predictor
	int predictor(int n,  int m, int nu, double **v, double **delP, double *A, double *b, double *u, double *c);

	// function to compute delC(v) - corrector
	void corrector(int n, int nu, double mu, double **v, double **delP);

	// function to compute the centering parameter
	double centering(int n, int m, int nu, double lim, double **v, double **delP);

	// fucntion to find the duality gap
	double dual_gap(int n, int nu, double **v);

	// update the variables
	void update(int n, int nu, int m, double **v, double **delP);

	// compute the stopping criterion
	double stopping(int n,  int m, int nu, double **v, double *A, double lenb, double *b,
				 double lenu, double *u, double lenc, double *c);

private:
	CSolver m_cSolver;

};

#endif /* CMATHUTILITY_H_ */
