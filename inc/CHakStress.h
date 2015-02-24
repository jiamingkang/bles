/*
 * CHakStress.h
 *
 *  Created on: 22 Feb 2015
 *      Author: knai20-admin
 */

#ifndef CHAKSTRESS_H_
#define CHAKSTRESS_H_

class CHakStress {
public:
	CHakStress();
	virtual ~CHakStress();

	//Methods
public:
	// calcualte the pnorm stresss
	double PnormStressCalc(mesh *inMesh, isoMat *inMat, double PN, double relax, double *alpha, double *disp, int itt, double *EStress, double *lsf, bool *fixed, double *Cstress);

	// calculate the loading for the adjoint problem
	void StressAdjoint(mesh *inMesh, isoMat *inMat, double *adjoint, double PN, double relax, double *alpha, double *disp, double NumDof, int itt, double *lsf, int *fixdof, double *EStress, bool *fixed);

	// calculate the stress sensitivies using least squares of integration points for AFG method
	void AFG_Stress_Sens(mesh *inMesh, boundary *bound_in, double *alpha, isoMat *inMat,  double *Nsens, double *disp, double *StAdj, int StressNumDual, int numCase, double *wgt, Coord *gCoord, double aMin, int mode, double Gspn, double PN, double RN, double *lsf, double *EStress, int *fixdof, bool *fixed, int itt);

	//Function to Calculate the Stress Sensitivities in all the guass points
	void StressSensitivity(mesh *inMesh, isoMat *inMat, double PN, double RN, double *alpha, double *disp, double *StAdj, int itt, double *EStress, double * ESens, Coord *gCoord, double *gSens, double *lsfGuass, double *lsf, int *fixdof, bool *fixed, double Gspn);

	// Function that calculates the stress sensitivity of a node by a least squares (2nd order) filter of near-by gauss points
	int LsensStress(Coord *pt, int xMax, int xMin, int yMax, int yMin, double aMin, double *alpha, double r2,
	                Coord *gCoord, double *gSens, Elem **Number, int wFlag, int numDual, double *out, double *lsfGuass, bool *fixed, double *EStress, double theta, double *Trust, double h);

	void GaStressSens_Q4(int *tnodes, double **prim, double **dual, double alpha, double h,
	                     isoMat *inMat, int Gcount, double *gSens, int numCase, int numDual, double *wgt, double Gspn, double PN, double relax);

	// calculate the stress sensitivies using least squares of integration points for AFG method
	void AFG_Stress_Sens_hole(mesh *inMesh, boundary *bound_in, double *alpha, isoMat *inMat,double *Nsens,
			double *disp, double *StAdj, int StressNumDual, int numCase, double *wgt, Coord *gCoord,
			double aMin, int mode, double Gspn, double PN, double RN, double *lsf, double *EStress,
			int *fixdof, bool *fixed, int itt, int h_count, int *h_index, int *h_EmapX, int *h_EmapY,
			int *h_posN, int *h_posE, double *h_Esens, double *h_Nsens, int ND);

};

#endif /* CHAKSTRESS_H_ */
