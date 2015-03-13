/*
 * CHakSparseSymmetric.h
 *
 *  Created on: 23 Feb 2015
 *      Author: knai20-admin
 */

#ifndef CHAKSPARSESYMMETRIC_H_
#define CHAKSPARSESYMMETRIC_H_

#include "common/CHakSparseMatrix.h"

// class for a sparse symmetric matrix  //original "matrix_spsy"
class CHakSparseSymmetric: public CHakSparseMatrix {
public:
	CHakSparseSymmetric();
	virtual ~CHakSparseSymmetric();

	int ord; // matrix order
	int len; // length of storage arrays
    int f_ind; // start of indexing, 0 = C, 1 = FORTRAN
	int *irn; // row number indicator
	int *jcn; // column number indicator
	double *A; // matrix entries
public:
	CHakSparseSymmetric(int n, int e) // constructor (set order and initial storage)
	{
		ord=n;
		len=e;
        f_ind=0; // assume zero
		A = new double [e];
		irn = new int [e];
		jcn = new int [e];
	}
	CHakSparseSymmetric(int n, int e, int f) // constructor (set order and initial storage ... and indexing order)
	{
		ord=n;
		len=e;
        f_ind=f;
		A = new double [e];
		irn = new int [e];
		jcn = new int [e];
	}
	CHakSparseSymmetric(const CHakSparseSymmetric &obj) // copy constructor
    {
        ord = obj.ord;
		len=obj.len;
        f_ind=obj.f_ind;
		A = new double [len];
		irn = new int [len];
		jcn = new int [len];

        int i;
        for(i=0;i<len;i++)
        {
            A[i] = obj.A[i];
            irn[i] = obj.irn[i];
            jcn[i] = obj.jcn[i];
        }
    }
	// overloaded constructor to make a new matrix by adding two existing ones
	CHakSparseSymmetric(CHakSparseSymmetric*, CHakSparseSymmetric*);
	CHakSparseSymmetric(CHakSparseSymmetric*, CHakSparseSymmetric*,double);

	~CHakSparseSymmetric()// destructor (delete memory)
	{
		len = 0;
		ord = 0;
		delete[] irn;
		delete[] jcn;
		delete[] A;
	}

	void mOut(int stop) {
		int i;
		stop = (stop<len) ? stop : len;
		for(i=0;i<stop;i++) {
			if(fabs(A[i])>0.0) {
			std::cout << "\n" << irn[i] << "\t" << jcn[i] << "\t" << A[i];
			}
		}
	}
	double* getA(){return(A);}
	int* get_irn(){return(irn);}
	int* get_jcn(){return(jcn);}
	void removeDof(int, int*, int); // function to remove fixed dof
	void setLen(int l){len=l;}
	int getOrd(){return(ord);}
	int getLen(){return(len);}
    int getf_ind(){return(f_ind);}
	void transpose() // method to swap row and column numbers
	{
		int *row, *col;
		row = jcn;
		col = irn;
		// simply swap the row * col pointers over
		irn = row;
		jcn = col;
	}

	// fucntion to compact entries in a square symmetric matrix
	void compact();
	// scale a spare symmetric matrix
	void scale(double);
	// function to expand a spsy matrix into a full matrix
	void expand(double*);
	// function to find max absolute diagonal entry
	double maxDiag();
    // multiply a sparse matrix by a vector
    void mvec(double *vin, double *vout);
};

#endif /* CHAKSPARSESYMMETRIC_H_ */
