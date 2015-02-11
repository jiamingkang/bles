/*
 * CHakMultiMaterial.cpp
 *
 *  Created on: 11 Feb 2015
 *      Author: Khalid Ismail
 */

#include "CHakMultiMaterial.h"

CHakMultiMaterial::CHakMultiMaterial() {
	// TODO Auto-generated constructor stub

}

CHakMultiMaterial::~CHakMultiMaterial() {
	// TODO Auto-generated destructor stub
}

//
// Implementation - Interfaces
//

// function to compute elastic modulus from two materials
double CHakMultiMaterial::HS_mat(double alpha, double hs_int, CHakMaterial *mat1, CHakMaterial *mat2)
{
    double kmax, kmin, gmax, gmin, khs, ghs;
    double fact, ftemp, ftemp2;
    double amin = 1.0-alpha;

    // bulk modulus
    fact = amin*mat1->m_k + alpha*mat2->m_k;
    ftemp = mat2->m_k - mat1->m_k;
    ftemp *= ftemp;
    ftemp *= amin*alpha;

    kmax = fact - (ftemp / (amin*mat2->m_k + alpha*mat1->m_k + mat2->m_g));
    kmin = fact - (ftemp / (amin*mat2->m_k + alpha*mat1->m_k + mat1->m_g));

    khs = hs_int*kmax + (1.0-hs_int)*kmin;

    // shear modulus
    fact = amin*mat1->m_g + alpha*mat2->m_g;
    ftemp = mat2->m_g - mat1->m_g;
    ftemp *= ftemp;
    ftemp *= amin*alpha;

    ftemp2 = (mat2->m_k * mat2->m_g) / (mat2->m_k + 2.0*mat2->m_g);
    gmax = fact - (ftemp / (amin*mat2->m_g + alpha*mat1->m_g + ftemp2));

    ftemp2 = (mat1->m_k * mat1->m_g) / (mat1->m_k + 2.0*mat1->m_g);
    gmin = fact - (ftemp / (amin*mat2->m_g + alpha*mat1->m_g + ftemp2));

    ghs = hs_int*gmax + (1.0-hs_int)*gmin;

    // Elastic modulus
    ftemp = (9.0*khs*ghs)/(3.0*khs + ghs);
    return (ftemp);
}
