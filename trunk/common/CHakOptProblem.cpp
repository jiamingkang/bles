/*
 * CHakOptProblem.cpp
 *
 *  Created on: Dec 19, 2014
 *      Author: jeehang
 */

#include "CHakOptProblem.h"

CHakOptProblem::CHakOptProblem() {
	// TODO Auto-generated constructor stub

	this->m_idObjType = 0;
	this->m_numConstraint = 0;
	this->m_pConstraint = NULL;
}

CHakOptProblem::~CHakOptProblem() {
	// TODO Auto-generated destructor stub
	if (this->m_pConstraint != NULL)
	{
		free(m_pConstraint); // should be replaced using delete (jeehanglee@gmail.com)
		this->m_pConstraint = NULL;
	}
}
