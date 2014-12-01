/*
 * CExMinimiseCompliance.h
 *
 *  Created on: Dec 1, 2014
 *      Author: jeehang
 */

#ifndef CEXMINIMISECOMPLIANCE_H_
#define CEXMINIMISECOMPLIANCE_H_

#include <string>

namespace topopt {

using std;

class CExMinimiseCompliance {
public:
	CExMinimiseCompliance();
	virtual ~CExMinimiseCompliance();

public:

private:
	std::string m_filename;	// input file name
};

} /* namespace topopt */

#endif /* CEXMINIMISECOMPLIANCE_H_ */
