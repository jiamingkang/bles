/*
 * CHakOptProblem.h
 *
 *  Created on: Dec 19, 2014
 *      Author: jeehang , Khalid Ismail
 */

#ifndef CHAKOPTPROBLEM_H_
#define CHAKOPTPROBLEM_H_

#include <iostream>
#include "CommonTypes.h"

class CHakOptProblem {
public:
	CHakOptProblem();
	virtual ~CHakOptProblem();

//
// Get/Set Properties
//
public:
	void SetObjectiveType(int idObjType) { m_idObjType = idObjType; }
	int GetObjectiveType() { return m_idObjType; }

	void SetNumOfConstraint(int numConstraint) { m_numConstraint = numConstraint; }
	int GetNumOfConstraint() { return m_numConstraint; }

	void SetConstraint(cnst *pConstraint) { m_pConstraint = pConstraint; }
	cnst *GetConstraint()
	{
		if (m_pConstraint != NULL)
			return m_pConstraint;

		std::cout << "Constraints are not ready";
		return NULL;
	}

private:
	// objective type identifier
	int m_idObjType; //obj

	// number of constraints
	int m_numConstraint; //num

	// pointer to array of constraint structs...
	cnst *m_pConstraint; //con
};

#endif /* CHAKOPTPROBLEM_H_ */
