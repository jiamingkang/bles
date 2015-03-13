/*
 * CHakSparseSymmetric.cpp
 *
 *  Created on: 23 Feb 2015
 *      Author: knai20-admin
 */

#include "CHakSparseSymmetric.h"

CHakSparseSymmetric::CHakSparseSymmetric() {
	// TODO Auto-generated constructor stub

}

CHakSparseSymmetric::~CHakSparseSymmetric() {
	// TODO Auto-generated destructor stub
}

// ------------------------------------------------------------ //
// Functions for the CHakSparseSymmetric class (sparse symmetic matrix) //
// ------------------------------------------------------------ //

// function to remove fixed dof from a spare symmetric matrix
void CHakSparseSymmetric::removeDof(int num, int *fixed, int plus)
{
	int i,j,k,temp;

    // create a map for the reduced dof
    int *dofMap;
    dofMap = new int [ord];

    for(i=0;i<ord;i++)
    {
        for(j=0;j<num;j++)
		{
            if(i == fixed[j])
			{
				dofMap[i] = -1;
				break;
			}

			else if( (i > fixed[j]) && (i < fixed[j+1]) )
			{
                dofMap[i] = i - j + plus; // plus one for fortran solvers
				break;
			}
		}
    }

    // reduce number of entries using the map
    temp = 0;
	for(i=0;i<len;i++)
	{
		j = dofMap[irn[i]-f_ind];
        k = dofMap[jcn[i]-f_ind];

		if( (j > -1) && (k > -1) ) // if dof not fixed
		{
			A[temp] = A[i];
			irn[temp] = j+f_ind;
			jcn[temp++] = k+f_ind;
		}
	}

	delete [] dofMap; // no longer required

	// Update the number of index2 entries
    f_ind += plus; // update indeing start
	len = temp;
	ord -= (num-2); // reduce order
}

// fucntion to compact entries in a square symmetric matrix
void CHakSparseSymmetric::compact()
{
	int i,m,o,temp;
	double sum;

	// consolodate duplicate entries
	temp=0;
	o=2*ord;
	for(i=0;i<len;i++)
	{
		if(irn[i] < o) // dont consider if already added before
		{
			sum = A[i];
			for(m=i+1;m<1000;m++)
			{
				// if row and col numbers are the same
				if( (irn[i] == irn[m]) && (jcn[i] == jcn[m]) )
				{
					sum += A[m]; // sum entries
					irn[m] = o;
					jcn[m] = o; // ensure duplicate entries only added once
				}
			}
			if(fabs(sum) > 1e-12)
			{
				A[temp] = sum;
				irn[temp] = irn[i];
				jcn[temp++] = jcn[i]; // store summed entries
			}
		}
	}

	// Update the number of index2 entries and point to reduced arrays
	len = temp;
}

// function to add one spsy matrix to another
// implemented as an overloaded constructor for the matrix_spsy class
CHakSparseSymmetric::CHakSparseSymmetric(CHakSparseSymmetric *one, CHakSparseSymmetric *two)
{
	int i;

    int f1 = one->getf_ind();
	int l1 = one->getLen();
	int ord1 = one->getOrd();
	int *r1 = one->get_irn();
	int *c1 = one->get_jcn();
	double *A1 = one->getA();

    int f2 = two->getf_ind();
	int l2 = two->getLen();
	int ord2 = two->getOrd();
	int *r2 = two->get_irn();
	int *c2 = two->get_jcn();
	double *A2 = two->getA();

    f_ind=f1;
    if(f1==f2){f2=0;} // no shift for 2nd input
    else{f2 = f1-f2;} // shift if differnet

	if(ord1 != ord2) // check input matrix orders are the same
	{
		std::cout << "Error constructing a spsy matrix by adding! ";
		std::cout << "Input matices are of different order";
		std::cout << "\n---return empty matrix";

		len = 0;
		ord = 0;
	}

	else
	{
		len = l1+l2;
		ord = ord1;
		A = new double [len];
		irn = new int [len];
		jcn = new int [len];

		for(i=0;i<l1;i++) // first input matrix
		{
			irn[i] = r1[i];
			jcn[i] = c1[i];
			A[i]   = A1[i];
		}

		for(i=0;i<l2;i++) // second input matrix
		{
			irn[i+l1] = r2[i]+f2;
			jcn[i+l1] = c2[i]+f2;
			A[i+l1]   = A2[i];
		}
	}

}

// function to add one spsy matrix to another (times by a factor)
// implemented as an overloaded constructor for the matrix_spsy class
CHakSparseSymmetric::CHakSparseSymmetric(CHakSparseSymmetric *one, CHakSparseSymmetric *two, double fact)
{
	int i;

    int f1 = one->getf_ind();
	int l1 = one->getLen();
	int ord1 = one->getOrd();
	int *r1 = one->get_irn();
	int *c1 = one->get_jcn();
	double *A1 = one->getA();

    int f2 = two->getf_ind();
	int l2 = two->getLen();
	int ord2 = two->getOrd();
	int *r2 = two->get_irn();
	int *c2 = two->get_jcn();
	double *A2 = two->getA();

    f_ind=f1;
    if(f1==f2){f2=0;} // no shift for 2nd input
    else{f2 = f1-f2;} // shift if differnet

	if(ord1 != ord2) // check input matrix orders are the same
	{
		std::cout << "Error constructing a spsy matrix by adding! ";
		std::cout << "Input matices are of different order";
		std::cout << "\n---return empty matrix";

		len = 0;
		ord = 0;
	}

	else
	{
		len = l1+l2;
		ord = ord1;
		A = new double [len];
		irn = new int [len];
		jcn = new int [len];

		for(i=0;i<l1;i++) // first input matrix
		{
			irn[i] = r1[i];
			jcn[i] = c1[i];
			A[i]   = A1[i];
		}

		for(i=0;i<l2;i++) // second input matrix
		{
			irn[i+l1] = r2[i]+f2;
			jcn[i+l1] = c2[i]+f2;
			A[i+l1]   = A2[i]*fact;
		}
	}

}

// function to expand a spsy matrix into a full matrix
void CHakSparseSymmetric::expand(double *out)
{
	int num,i,r,c;

	// initialise to zero
	num = ord*ord;
	for(i=0;i<num;i++)
	{
		out[i] = 0.0;
	}

	// add stored entries to the output matrix
	num = len;
	for(i=0; i<num; i++)
	{
		r = irn[i];
		c = jcn[i];
		out[index2(r,c,ord)] += A[i];
		if(r != c) {
			out[index2(c,r,ord)] += A[i];
		}

	}
}

// function to find max absolute diagonal entry
double CHakSparseSymmetric::maxDiag()
{
	int i;
	double ft1;

	double max = 0.0;

	//find absolutue maximum diagonal value
	for(i=0;i<len;i++)
	{
		if(irn[i] == jcn[i])
		{
			ft1 = fabs(A[i]);
			max = (max > ft1) ? max : ft1;
		}
	}

	return(max);
}

// multiply a sparse matrix by a vector
void CHakSparseSymmetric::mvec(double *vin, double *vout)
{
	int i,rn,cn;

	// reset vout to zero
	for(i=0;i<ord;i++){vout[i]=0.0;}

	// multiply
	for(i=0;i<len;i++)
	{
		rn =irn[i]-1;
		cn = jcn[i]-1;

		vout[rn] += A[i] * vin[cn];

		// symmetry
		if(rn != cn)
		{
			vout[cn] += A[i] * vin[rn];
		}
	}
}
