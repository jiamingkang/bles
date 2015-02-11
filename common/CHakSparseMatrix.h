/*
 * CHakSparseMatrix.h
 *
 *  Created on: Dec 19, 2014
 *      Author: Khalid Ismail, jeehang
 */

#ifndef CHAKSPARSEMATRIX_H_
#define CHAKSPARSEMATRIX_H_

#include <iostream>

class CHakSparseMatrix {

public:
	CHakSparseMatrix();
	virtual ~CHakSparseMatrix();

// Get/Set properties
public:
	void SetNumOfEntry(int numEntry)
	{
		m_numEntry = numEntry;
	}

	int GetNumOfEntry()
	{
		return m_numEntry;
	}

	void SetRowIndicator(int indRow)
	{
		m_indRow = indRow;
	}

	int GetRowIndicator()
	{
		return m_indRow;
	}

	void SetColIndicator(int indCol)
	{
		m_indCol = indCol;
	}

	int GetcolIndicator()
	{
		return m_indCol;
	}

	void SetMatrixEntryValue(double *pMatEntry)
	{
		m_pMatEntry = pMatEntry;
	}

	double* GetMatrixEntryValue()
	{
		if (m_pMatEntry != NULL)
			return m_pMatEntry;

		std::cout << "Matrix Entry Value is not ready";
		return NULL;

	}

// members
//private:
public:
	// number of entries (in arrays)
	int m_numEntry; //ne

	// indicators for row and column
	//int m_indRow; //*irn "To be used in encapsulation"
	//int m_indCol; //*jcn  "To be used in encapsulation"
	int *m_indRow; //*irn
	int *m_indCol; //*jcn

	// matrix entry value
	double *m_pMatEntry; //*A
};

#endif /* CHAKSPARSEMATRIX_H_ */
