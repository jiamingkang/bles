/*
 * CHakSparseUnsymmetric.h
 *
 *  Created on: 23 Feb 2015
 *      Author: knai20-admin
 */

#ifndef CHAKSPARSEUNSYMMETRIC_H_
#define CHAKSPARSEUNSYMMETRIC_H_

#include "common/CHakSparseMatrix.h"

class CHakSparseUnsymmetric: public CHakSparseMatrix {
public:
	CHakSparseUnsymmetric();
	virtual ~CHakSparseUnsymmetric();

	int row; // number rows
	int col; // number columns
	int len; // length of storage arrays
	int *irn; // row number indicator
	int *jcn; // column number indicator
	double *A; // matrix entries
public:
	matrix_spun(int r, int c, int e ) // constructor (set rows * cols and initial storage)
	{
		row=r;
		col=c;
		len=e;
		A = new double [e];
		irn = new int [e];
		jcn = new int [e];
	}

	~matrix_spun() // destructor (delete memory)
	{
		len = 0;
		row = 0;
		col = 0;
		delete[] irn;
		delete[] jcn;
		delete[] A;
	}

	void mOut(int stop) {
		int i;
		stop = (stop<len) ? stop : len;
		for(i=0;i<stop;i++) {
			std::cout << "\n" << irn[i] << "\t" << jcn[i] << "\t" << A[i];
		}
	}

	double* getA(){return(A);}
	int* get_irn(){return(irn);}
	int* get_jcn(){return(jcn);}
	void setLen(int l){len=l;}
	int getRow(){return(row);}
	int getCol(){return(col);}
	int getLen(){return(len);}

	void transpose() // method to swap row and column numbers
	{
		int *r, *c, t;
		r = jcn;
		c = irn;
		// simply swap the row * col pointers over
		irn = r;
		jcn = c;
		// and row and col lengths
		t = row;
		row = col;
		col = t;
	}

	// function to remove rows from a matrix
	void removeCols(int, int*);
	// function to expand a spsy matrix into a full matrix
	void expand(double*);
	// function to remove zeros from a spsy matrix
	void remZeros();
	// increase row numbers of a spun matrix
	void upIrn(int);
};

#endif /* CHAKSPARSEUNSYMMETRIC_H_ */
