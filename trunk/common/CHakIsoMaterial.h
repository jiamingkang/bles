/*
	CHakIsoMaterial.h

 	Created on: Dec 17, 2014
	Author: JeeHang Lee

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

#ifndef CHAKISOMATERIAL_H_
#define CHAKISOMATERIAL_H_

// isotropic material
class CHakIsoMaterial
{
public:
	CHakIsoMaterial() {};
	virtual ~CHakisoMaterial() {};

// properties
public:
	// todo: add get/set properties

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
}

#endif // CHAKISOMATERIAL_H_