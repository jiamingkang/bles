/*
 * CExMinimiseCompliance.h
 *
 *  Created on: Dec 1, 2014
 *      Author: jeehang
 */

#ifndef CEXMINIMISECOMPLIANCE_H_
#define CEXMINIMISECOMPLIANCE_H_

#include <string>

using std;

class CExMinimiseCompliance {
public:
	CExMinimiseCompliance();
	virtual ~CExMinimiseCompliance();

public:
	void Solve(char* arg, char* filename);

private:
	std::string m_filename;	// input file name
};

#endif /* CEXMINIMISECOMPLIANCE_H_ */
