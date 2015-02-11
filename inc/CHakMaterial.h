/*
   CHakMaterial.h

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

#ifndef CHakMaterial_H_
#define CHakMaterial_H_

#include "CommonTypes.h"

class CHakMaterial {
public:
	CHakMaterial();
	virtual ~CHakMaterial();

public:
// private:
	double m_e;		// modulus?
	double m_v;		// Poisson's ratio?
	double m_rho;	// density?
	double m_k;		// bulk?
	double m_g;		// shear?
	double m_mat[9];	// material property matrix

// private:
	// number of materials 
	int m_cntMat;
};

#endif /* CHakMaterial_H_ */
