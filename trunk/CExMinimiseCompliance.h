/*
 * CExMinimiseCompliance.h
 *
 *  Created on: Dec 1, 2014
 *      Author: jeehang
 */

#ifndef CEXMINIMISECOMPLIANCE_H_
#define CEXMINIMISECOMPLIANCE_H_

#include <string>

//using std;

class CExMinimiseCompliance {
public:
	CExMinimiseCompliance();
	virtual ~CExMinimiseCompliance();

//
// Interfaces
//
public:
	// Read the input file & assign initial values
	// returns -1 when starting process fails.
	int StartProcess(std::string inputfile);

	// Finite Element Anaysis
	// returns -1 when exceptions occurs thus process fails.
	int AnalyseFiniteElement();

	// Sensitivity 
	int AdjustSensitivity();

	// Optimisation
	int Optimisation();
	
	// main, dummy function for the compatibility
	void Solve(char* arg, char* filename);

//
// Utilities
//
protected:


private:
	std::string m_filename;	// input file name

	// Read input file and assign initial values specified in the file 
	CHakInput m_input;

	// handle output files
	CHakOutput m_output;
};

#endif /* CEXMINIMISECOMPLIANCE_H_ */
