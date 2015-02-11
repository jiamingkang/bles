/*
 * CHakMultiMaterial.h
 *
 *  Created on: 11 Feb 2015
 *      Author: Khalid Ismail
 */

#ifndef CHAKMULTIMATERIAL_H_
#define CHAKMULTIMATERIAL_H_

class CHakMultiMaterial {
public:
	CHakMultiMaterial();
	virtual ~CHakMultiMaterial();

public:
	// function to compute elastic modulus from two materials
	double HS_mat(double alpha, double hs_int, CHakMaterial *mat1, CHakMaterial *mat2);
};

#endif /* CHAKMULTIMATERIAL_H_ */
