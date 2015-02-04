/*
 * CHakOptControl.h
 *
 *  Created on: Dec 19, 2014
 *      Author: jeehang , Khalid Ismail
 */

#ifndef CHAKOPTCONTROL_H_
#define CHAKOPTCONTROL_H_

class CHakOptControl {
public:
	CHakOptControl();
	virtual ~CHakOptControl();

//
// Get/Set properties
//
public:
	void SetMaxIteration(int maxIter)
	{
		m_maxIter = maxIter;
	}

	int GetMaxIteration()
	{
		return m_maxIter;
	}

	void SetAmoutOfOutput(int outInfo)
	{
		m_outInfo = outInfo;
	}

	int GetAmountOfOutput()
	{
		return m_outInfo;
	}

	void SetNarrowBandWidth(double lband)
	{
		m_lband = lband;
	}

	double GetNarrowBandWidth()
	{
		return m_lband;
	}

	void SetConvergenceCriterion(double gmConv)
	{
		m_gmConv = gmConv;
	}

	double GetConvergenceCriterion()
	{
		return m_gmConv;
	}

	void SetMinAreaRatio(double minArea)
	{
		m_minArea = minArea;
	}

	double GetMinAreaRatio()
	{
		return m_minArea;
	}

	void SetMinMassRatio(double minMass)
	{
		m_minMass = minMass;
	}

	double GetMinMassRatio()
	{
		return m_minMass;
	}

//
// members
//
private:
	// maximum number of interations (should > 0)
	int m_maxIter;

	// amount of output (1->3, less->more)
	int m_outInfo;

	// convergence criterion (gamma)
	double m_gmConv;

	// narrow band width (2h -> large)
	double m_lband;

	// minimum area ratio (i.e. small or zero) - stiffness
	double m_minArea;

	// minimum area ratio for mass
	double m_minMass;
};

#endif /* CHAKOPTCONTROL_H_ */
