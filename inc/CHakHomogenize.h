/*
 * CHakHomogenize.h
 *
 *  Created on: 22 Feb 2015
 *      Author: knai20-admin
 */

#ifndef CHAKHOMOGENIZE_H_
#define CHAKHOMOGENIZE_H_

const static double BHom[24] = {-1.0, 0.0, 1.0, 0.0, 1.0, 0.0, -1.0, 0.0,
                                0.0, -1.0, 0.0, -1.0, 0.0, 1.0, 0.0, 1.0,
                                -1.0, -1.0, -1.0, 1.0, 1.0, 1.0, 1.0, -1.0};

const static int Esymm_map[5] = {0, 1, 3, 4, 8};


class CHakHomogenize {
public:
	CHakHomogenize();
	virtual ~CHakHomogenize();

public:
	// function to homogenize the
	void homogenize(mesh *inMesh, isoMat *mat1, sp_mat *K1, double *alpha, int **fixed, int freeDof, isoMat *homMat, double *disp);

	// function to homogenize the material properties of the microstructure - using constraint equn periodic BCs
	void homogenize2(mesh *inMesh, isoMat *mat1, sp_mat *K1, double *alpha, int **fixed, int freeDof, sp_mat *CE, isoMat *homMat, double *disp, int pinfo);

	// function to compute homogenized BCs
	void hom_bc(mesh *inMesh, int **fixed, int *freeDof);

	// function to compute homogenized BCs - Lagrange multiplier method
	void hom_bc2(mesh *inMesh, sp_mat *HomBC);
};

#endif /* CHAKHOMOGENIZE_H_ */
