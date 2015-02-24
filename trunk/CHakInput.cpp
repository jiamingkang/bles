/*
 * CHakInput.cpp
 *
 *  Created on: Nov 24, 2014
 *      Author: Khalid Ismail, jeehang
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "CHakMesh.h"			// originally Numbering.h
#include "CHakLevelSet.h"		// originally Levels.h --> should be changed to COptimisation.h (jeehanglee@gmail.com)
#include "CHakMathUtility.h"	// originally ABFG.h
#include "CHakInput.h"
#include "CHakFiniteElement.h"

//
// Constructur / Destructor
//
CHakInput::CHakInput() {
	// TODO Auto-generated constructor stub
	m_pfilename = NULL;
	_SetDefault();
}

CHakInput::CHakInput(char *filename)
{ 
	m_pfilename = filename; 
	_SetDefault();
}

CHakInput::~CHakInput() {
	// TODO Auto-generated destructor stub
}

void CHakInput::_SetDefault() 
{
	// set default objective and number of constraint in problem definition
	m_problem.SetObjectiveType(0);
	m_problem.SetNumOfConstraint(0);

	// set default controls
	m_control.SetMaxIteration(200);
	m_control.SetAmoutOfOutput(1);
	m_control.SetConvergenceCriterion(1.0E-3);
	m_control.SetNarrowBandWidth(6.0);
	m_control.SetMinAreaRatio(1.0E-6);
	m_control.SetMinMassRatio(0.0);

	// initialize number load cases as 0
	m_numCase = 0;
}

//
// Implementation - Interfaces
//
int CHakInput::icmpfunc (const void * p1, const void * p2)
{
	int c = ( *(int *) p1 > *(int *) p2 ) ? 1 : -1;
	return c;
}

int CHakInput::ReadFile(char *datafile)
{
	_SetInputFile(datafile);

	CHakFiniteElement fem;
	int i, j, nd, nd2, ind, end;

	// temp data arrays
	// count number of materials
	int		matCnt = 0; 
	double	*load;
	bool	*fixdof_temp;
	double	*lm_temp=0;

	// temp for initial holes
	// number of circular holes
	int NumHole = 0;		
	// array to store circular hole data
	CirH *holes = (CirH *) malloc(MAX_HOLES * sizeof(CirH)); 
	// number of rectangular holes
	int NumRect=0;			
	// array to store rectangular hole data (max & min x&y coords)
	Coord *Rect = (Coord *) malloc(2 * MAX_HOLES * sizeof(Coord)); 

	//
	// Reading Start!
	//
	FILE *infile = fopen(datafile, "r"); // Try to open the file
	if (infile == NULL)
	{
		// If file does not exist then abort
		printf("\nCould not find the input file! - Abort\n");
		return -1;
	}

	// read an entire line into memory
	char buf[MAX_CHARS_PER_LINE];
	int reread = 0;
	int meshFound = 0;		// indicate when *mesh has been found
	m_material[0].m_e = -1.0;// used to indicate when *mat has been found

	// read each line of the file
	while (!feof(infile))
	{
		// array to store memory addresses of the tokens in buf
		const char* token[MAX_TOKENS_PER_LINE]; // = {}; // initialize to 0

		// get next line (unless we want to reread the line)
		if (reread == 0) 
			fgets(buf, MAX_CHARS_PER_LINE, infile);
		else 
			reread = 0;

		token[0] = strtok(buf, DELIMITER); // first token

		// read in all tokens
		if (token[0]) // zero if line is blank
		{
			for (i = 1; i < MAX_TOKENS_PER_LINE; i++)
			{
				token[i] = strtok(NULL, DELIMITER); // subsequent tokens
				if (!token[i]) 
					break; // no more tokens
			}
		}

		// search tokens for meaningful data
		if (buf[0] == '*' && buf[1] != '*') // keyword line
		{
			char *lead = strtok(buf, DELIMITER); // first token (possible keyword)

			// if keyword is mesh
			if(strncasecmp(lead, "*mesh", 5) == 0)
			{
				end = 0;
				if (meshFound == 1)
				{
					printf("\nWarning: Found second *mesh! - Ignoring");
					end = 1;
				}

				while (end == 0 && !feof(infile))
				{
					fgets(buf, MAX_CHARS_PER_LINE, infile); // get next line

					if(buf[0]!='*' && buf[1]!='*') // check for commented out line
					{
						// split line into tokens
						token[0] = strtok(buf, DELIMITER); // first token
						for (i = 1;i < 4; i++)
						{
							token[i] = strtok(NULL, DELIMITER); // subsequent tokens
							if (!token[i]) 
								break; // no more tokens
						}

						if(i < 2)
						{
							printf("\nError in input: Not enough data for *mesh! - Abort\n");
							return -1;
						}

						// read data into the inMesh structure
						sscanf(token[0], "%i", &m_mesh.m_elemX);
						sscanf(token[1], "%i", &m_mesh.m_elemY);
						if (i > 2)
						{
							sscanf(token[2], "%lf", &m_mesh.m_lenEdge);
						}
						else 
						{
							m_mesh.m_lenEdge = 1.0; // default to 1.0
						}
						
						if (i > 3)	{
							sscanf(token[3], "%lf", &m_mesh.m_thinkness);
						}
						else {
							m_mesh.m_thinkness = 1.0; // default to 1.0
						}

						// check data
						if(m_mesh.m_elemX < 1 || m_mesh.m_elemY < 1)
						{
							printf("\nError in input: elems %i x %i! - Abort\n", m_mesh.m_elemX, m_mesh.m_elemY);
							return -1;
						}
						if(m_mesh.m_lenEdge <= 0.0)
						{
							printf("\nError in input: h = %lf! - Abort\n", m_mesh.m_lenEdge);
							return -1;
						}
						if(m_mesh.m_thinkness <= 0.0)
						{
							printf("\nError in input: t = %lf! - Abort\n", m_mesh.m_thinkness);
							return -1;
						}

						// define all other data in CHakMesh
						m_mesh.m_numElem = m_mesh.m_elemX * m_mesh.m_elemY;
						m_mesh.m_nodeX = 3 + m_mesh.m_elemX;
						m_mesh.m_nodeY = 3 + m_mesh.m_elemY;	// one extra node is each direction
						m_mesh.m_numNodes = (1 + m_mesh.m_elemX) * (1 + m_mesh.m_elemY);
						m_mesh.m_maxX = m_mesh.m_elemX * m_mesh.m_lenEdge;
						m_mesh.m_maxY = m_mesh.m_elemY * m_mesh.m_lenEdge;
						m_mesh.m_pNodeCoord = (Coord *) calloc(m_mesh.m_numNodes, sizeof(Coord));
						m_mesh.m_tolerance = m_mesh.m_lenEdge * 1.0E-6; // rounding tollerance
						m_mesh.m_pMaterialType = (int *) calloc(m_mesh.m_numElem, sizeof(int)); // default to first material
						m_mesh.m_bars = false; // no bars (yet)
						m_mesh.m_bDesignableBc = false; // no designable bsc (yet)
						m_mesh.m_bDesignableMaterial = false; // no designable material (yet)

						// 2d array memory allocation (array of pointers to pointers)
						m_mesh.m_pNumber = (Elem **) malloc(m_mesh.m_elemX *sizeof(Elem *));
						for (i = 0; i< m_mesh.m_elemX; i++) 
						{
							m_mesh.m_pNumber[i] = (Elem *) malloc(m_mesh.m_elemY * sizeof(Elem));
						}

						m_mesh.m_pNodes2D = (int **) malloc(m_mesh.m_nodeX * sizeof(int *));
						for (i = 0;i < m_mesh.m_nodeX; i++)
						{
							m_mesh.m_pNodes2D[i] = (int *) malloc(m_mesh.m_nodeY * sizeof(int));
						}

						// call fucntions to compute data arrays in inMesh
						m_mesh.Numbering();
						m_mesh.Coordinates();
						m_mesh.NodeNums2D();

						// now mesh has been defined - set some memory
						fixdof_temp = (bool *) calloc(2 * m_mesh.m_numNodes, sizeof(bool));
						m_levelset.m_pFixedLsf = (bool *) calloc(m_mesh.m_numNodes, sizeof(bool));

						meshFound = 1;
						end = 1; // found data line
					}
					else if(buf[0]=='*'&& buf[1]!='*')
					{
						printf("\nError in input: No data line for *mesh! - Abort\n");
						return -1;
					}
				}
			}
			//
			// jeehanglee@gmail.com: applying OOD till here... (20141223)
			//

			// if keyword is bars
			// jeehanglee@gmail.com: may go to _ReadBarsInfo() ?
			else if (strncasecmp(lead, "*bars", 5)==0)
			{
				// check that mesh has defined
				if(meshFound == 0)
				{
					printf("\nError in input: *bars before *mesh! - Abort\n");
					return -1;
				}

				// set variables and memory for bar elements
				m_mesh.m_bars = true;	// problem contains bars
				nd = (m_mesh.m_elemX * (1 + m_mesh.m_elemY) + m_mesh.m_elemY * (1 + m_mesh.m_elemY));
				m_mesh.m_numBars = nd;	// total number of bars.
				m_mesh.m_pBarNums = (Bseg *) malloc(nd * sizeof(Bseg));
				m_mesh.m_pBarAreas = (double *) malloc(nd * sizeof(double));

				// read in remaining data from the dataline
				end = 0;
				while (end == 0 && !feof(infile))
				{
					fgets(buf, MAX_CHARS_PER_LINE, infile); // get next line

					if (buf[0] != '*' && buf[1] != '*') // check for commented out line
					{
						// split line into tokens
						token[0] = strtok(buf, DELIMITER); // first token
						for (i = 1;i < 3; i++)
						{
							token[i] = strtok(NULL, DELIMITER); // subsequent tokens
							if (!token[i])
								break; // no more tokens
						}

						if(i < 3)
						{
							printf("\nError in input: Not enough data for *bars! - Abort\n");
							return -1;
						}

						// read data into the inMesh structure
						sscanf(token[0], "%lf", &m_mesh.m_barMin);
						sscanf(token[1], "%lf", &m_mesh.m_barMax);
						sscanf(token[2], "%i", &j);
						m_mesh.m_pBarMaterialType = &m_material[j];

						// check data
						if (m_mesh.m_barMin < 1.0e-20)
						{
							printf("\nWarning: Minimum bar area too small or -ve = using default (1e-3)\n");
							m_mesh.m_barMin = 1.0e-3;
						}

						// call bar numbering function
						// cmesh.Bar_numbering(inMesh);
						m_mesh.NumberingBarElements();

						end = 1; // found data line
					}
					else if (buf[0] == '*' && buf[1] != '*')
					{
						printf("\nError in input: No data line for *mat! - Abort\n");
						return -1;
					}
				}
			}

			// if keyword is material
			// jeehanglee@gmail.com -> may go to _ReadMaterialInfo() ?
			else if (strncasecmp(lead, "*mat", 4)==0)
			{
				end = 0;

				if (matCnt == 5)
				{
					printf("\nError! Too many materials defined (max = 5) - Abort!");
					return -1;
				}

				// check that mesh has defined
				if (meshFound == 0)
				{
					printf("\nError in input: *mat before *mesh! - Abort\n");
					return -1;
				}

				while (end == 0 && !feof(infile))
				{
					fgets(buf, MAX_CHARS_PER_LINE, infile); // get next line

					if(buf[0]!='*' && buf[1]!='*') // check for commented out line
					{
						// split line into tokens
						token[0] = strtok(buf, DELIMITER); // first token
						for (i = 1; i < 3; i++)
						{
							token[i] = strtok(NULL, DELIMITER); // subsequent tokens
							if (!token[i])
								break; // no more tokens
						}

						if (i < 2)
						{
							printf("\nError in input: Not enough data for *mat! - Abort\n");
							return -1;
						}

						// read data into the m_material structure
						sscanf(token[0],"%lf",&m_material[matCnt].m_e);
						sscanf(token[1],"%lf",&m_material[matCnt].m_v);
						if (i > 2)
						{
							sscanf(token[2],"%lf",&m_material[matCnt].m_rho);
						}
						else
						{
							m_material[matCnt].m_rho = 1.0; // default to 1.0
						}

						// check data
						if (m_material[matCnt].m_e < 1.0e-20)
						{
							printf("\nError in input: Young's modulus = %12.4e! - Abort\n",m_material[matCnt].m_e);
							return -1;
						}
						else if (m_material[matCnt].m_v < 0.5001 && m_material[matCnt].m_v > 0.4999)
						{
							printf("\nError in input: Poisson's ratio = %12.4e! - Abort\n",m_material[matCnt].m_v);
							return -1;
						}
						else if(m_material[matCnt].m_rho < 1.0e-20)
						{
							printf("\nError in input: density = %12.4e! - Abort\n",m_material[matCnt].m_rho);
							return -1;
						}

						// compute plane stress constants
						double e11 = m_material[matCnt].m_e / ( 1.0-(m_material[matCnt].m_v * m_material[matCnt].m_v) );
						double v12 = e11 * m_material[matCnt].m_v;
						double g33 = 0.5 * e11 * (1.0 - m_material[matCnt].m_v);

						// compute the material matrix
						m_material[matCnt].m_mat[0] = e11;
						m_material[matCnt].m_mat[1] = v12;
						m_material[matCnt].m_mat[2] = 0.0;
						m_material[matCnt].m_mat[3] = v12;
						m_material[matCnt].m_mat[4] = e11;
						m_material[matCnt].m_mat[5] = 0.0;
						m_material[matCnt].m_mat[6] = 0.0;
						m_material[matCnt].m_mat[7] = 0.0;
						m_material[matCnt].m_mat[8] = g33;

                        // compute bulk and shear moduli
                        m_material[matCnt].m_k = m_material[matCnt].m_e / (3.0 * (1.0 - 2.0 * m_material[matCnt].m_v));
                        m_material[matCnt].m_g = m_material[matCnt].m_e / (2.0 * (1.0 + m_material[matCnt].m_v));

						matCnt++; // increase material count

						end = 1; // found data line
					}
					else if (buf[0] == '*' && buf[1] != '*')
					{
						printf("\nError in input: No data line for *mat! - Abort\n");
						return -1;
					}
				}
			}

			// if keyword is mat_def (material definition)
			// 	jeehanglee@gmail.com: _ReadMaterialDefinition()
			else if(strncasecmp(lead, "*def_mat", 8) == 0)
			{
				// check that mesh has defined
				if (meshFound == 0)
				{
					printf("\nError in input: *def_mat before *mesh! - Abort\n");
					return -1;
				}

				end = 0;
				while (end == 0 && !feof(infile))
				{
					fgets(buf, MAX_CHARS_PER_LINE, infile); // get next line

					double h = m_mesh.m_lenEdge;
					double tol = m_mesh.m_tolerance;
					double xmin,xmax,ymin,ymax;
					int mtype;

					if (buf[0] != '*' && buf[1] != '*') // check for commented out line
					{
						// split line into tokens
						token[0] = strtok(buf, DELIMITER); // first token
						for (i = 1 ;i < 5 ; i++)
						{
							token[i] = strtok(NULL, DELIMITER); // subsequent tokens
							if (!token[i])
								break; // no more tokens
						}

						if (i < 5)
						{
							printf("\nError in input: not enough data for *def_mat! - Abort\n");
							return -1;
						}

						sscanf(token[0],"%lf",&xmin); // x min coord
						sscanf(token[1],"%lf",&xmax); // x max coord
						sscanf(token[2],"%lf",&ymin); // y min coord
						sscanf(token[3],"%lf",&ymax); // y max coord
						sscanf(token[4],"%i",&mtype); // material type

						// find all elements completely enclosed in area

						// use coords to work out range of elements effected
						xmin -= tol;
						ymin -= tol;
						xmax += tol;
						ymax += tol;

						int ex_min = (int)ceil(xmin/h);
						int ex_max = (int)floor(xmax/h);
						int ey_min = (int)ceil(ymin/h);
						int ey_max = (int)floor(ymax/h);

						int n,m,num;
						for (m = ey_min; m < ey_max; m++)
						{
							for (n = ex_min; n < ex_max; n++)
							{
								num = m_mesh.m_pNumber[n][m].n; // Element number
								m_mesh.m_pMaterialType[num] = mtype; // material type
							}
						}
					}

					else if (buf[0] == '*' && buf[1] != '*')
					{
						end = 1; // stop when next keyword is found
						reread = 1; // nned to reread the keyword
					}
				}
			}

			// if keyword is designable material
			// jeehanglee@gmail.com: _ReadDesignableMaterial()
			else if (strncasecmp(lead, "*des_mat", 8) == 0)
			{
				// check that mesh has defined
				if (meshFound == 0)
				{
					printf("\nError in input: *des_mat before *mesh! - Abort\n");
					return -1;
				}

				// next two tokens are the material numbers
				int m1, m2;
				sscanf(token[1],"%i",&m1); // material 1
				sscanf(token[2],"%i",&m2); // material 2

                // check for simultaneous or sequential opt of materials
                if (i < 4)
                {
                	m_mesh.m_bDmSim = true;
                } // default is simultaneous
                else if (strncasecmp(token[3], "sim", 3) == 0)
                {
                	m_mesh.m_bDmSim = true;
                }
                else if (strncasecmp(token[3], "seq", 3) == 0)
                {
                	m_mesh.m_bDmSim = false;
                }
                else
                {
					printf("\nCannot determine sim or seq in *des_mat! - Abort\n");
					return -1;
				}

                // check for linear or H-S bound material model
                if (i < 5)
                {
                	m_mesh.m_bMaterialLin = true;
                } // default is linear model
                else if (strncasecmp(token[4], "lin", 3) == 0)
                {
                	m_mesh.m_bMaterialLin = true;
                }
                else if (strncasecmp(token[4], "h-s", 3) == 0)
                {
                	m_mesh.m_bMaterialLin = false;
                }
                else
                {
					printf("\nCannot determine linear of H-s bound material model in *des_mat! - Abort\n");
					return -1;
				}

				// check these materail exist
				if (m1 > matCnt || m2 > matCnt)
				{
					printf("\nMaterial not defined in *des_mat! - Abort\n");
					return -1;
				}

				// set materials
				m_mesh.m_materialOne = m1;
				m_mesh.m_materialTwo = m2;

				int numElem = m_mesh.m_numElem;
				double h = m_mesh.m_lenEdge;
				double tol = m_mesh.m_tolerance;
				double xmin,xmax,ymin,ymax;
				bool *exclude = (bool *) calloc(numElem, sizeof(bool));

				// now check for exclusion zones
				end = 0;
				while (end == 0 && !feof(infile))
				{
					fgets(buf, MAX_CHARS_PER_LINE, infile); // get next line

					if (buf[0] != '*' && buf[1] != '*') // check for commented out line
					{
						// split line into tokens
						token[0] = strtok(buf, DELIMITER); // first token
						for (i = 1; i < 4; i++)
						{
							token[i] = strtok(NULL, DELIMITER); // subsequent tokens
							if (!token[i])
								break; // no more tokens
						}

						if (i < 4)
						{
							printf("\nError in input: not enough data for *des_mat exclusion zone! - Abort\n");
							return -1;
						}

						sscanf(token[0], "%lf", &xmin); // x min coord
						sscanf(token[1], "%lf", &xmax); // x max coord
						sscanf(token[2], "%lf", &ymin); // y min coord
						sscanf(token[3], "%lf", &ymax); // y max coord

						// find all elements completely enclosed in area

						// use coords to work out range of elements effected
						xmin -= tol; ymin -= tol; xmax += tol; ymax += tol;
						int ex_min = (int) ceil(xmin/h);
						int ex_max = (int) floor(xmax/h);
						int ey_min = (int) ceil(ymin/h);
						int ey_max = (int) floor(ymax/h);

						int n,m,num;
						for (m = ey_min; m < ey_max; m++)
						{
							for (n = ex_min; n < ex_max; n++)
							{
								num = m_mesh.m_pNumber[n][m].n; // Element number
								exclude[num] = true;
							}
						}
					}

					else if (buf[0] == '*' && buf[1] != '*')
					{
						end=1; // stop when next keyword is found
						reread=1; // nned to reread the keyword
					}
				}

				// set array of elements with designable material
				int *enum_temp = (int *) malloc(numElem * sizeof(int));
				nd = 0; // count number of elements
				for (i = 0; i < numElem; i++)
				{
					if (!exclude[i])
					{
						enum_temp[nd++] = i;
					}
				}
				free(exclude);

				if (nd > 0)
				{
					m_mesh.m_numDesignableMaterialt = nd; // total number of variables
					if (nd < numElem)
					{
						enum_temp = (int *) realloc(enum_temp, nd * sizeof(int));
					}
					m_mesh.m_pMaterialElems = enum_temp;
					m_mesh.m_bDesignableMaterial = true;
				}
				else
				{
					m_mesh.m_numDesignableMaterialt = 0; // total number of variables
					free(enum_temp);
				}
			}

			// if keyword is mass - defines lumped mass at nodes
			else if (strncasecmp(lead, "*mass", 5) == 0)
			{
				// check that mesh has defined
				if (meshFound == 0)
				{
					printf("\nError in input: *mass before *mesh! - Abort\n");
					return -1;
				}

				end = 0;
				int nn = m_mesh.m_numNodes;
				double xtemp, ytemp, mass;
				if (!lm_temp)
				{
					// if not done already - set memory for lumped masses
					lm_temp = (double *) calloc(nn * NUM_DOF, sizeof(double));
				}

				while (end == 0 && !feof(infile))
				{
					fgets(buf,MAX_CHARS_PER_LINE,infile); // get next line

					if (buf[0] != '*' && buf[1] != '*') // check for commented out line
					{
						// split line into tokens
						token[0] = strtok(buf, DELIMITER); // first token
						for (i = 1;i < 3;i++)
						{
							token[i] = strtok(NULL, DELIMITER); // subsequent tokens
							if (!token[i])
								break; // no more tokens
						}

						if (i < 3)
						{
							printf("\nError in input: not enough data for *mass! - Abort\n");
							return -1;
						}

						sscanf(token[0], "%lf", &xtemp); // x coord
						sscanf(token[1], "%lf", &ytemp); // y coord
						sscanf(token[2], "%lf", &mass);  // mass

						// find closest node to the co-ordinates
						// nd = cmesh.closeNode(inMesh,xtemp,ytemp);
						nd = m_mesh.CloseNode(xtemp, ytemp);

						// Now add mass to specified dof
						nd *= NUM_DOF;
						for (i = 0; i < NUM_DOF; i++)
						{
							lm_temp[nd++] += mass;
						}
					}

					else if (buf[0]=='*' && buf[1]!='*')
					{
						end = 1; // stop when next keyword is found
						reread = 1; // nned to reread the keyword
					}
				}
			}

			// if keyword is boundary
			// jeehanglee@gmail.com: _ReadBoundaryInfo()
			else if(strncasecmp(lead, "*bound", 6)==0)
			{
				// check that mesh has defined
				if(meshFound==0)
				{
					printf("\nError in input: *bound before *mesh! - Abort\n");
					return -1;
				}

				int *flag = (int *) malloc(NUM_DOF * sizeof(int));
				int nn = m_mesh.m_numNodes;
				Coord *cp = m_mesh.m_pNodeCoord;

				// if type is point
				if(strncasecmp(token[1], "point", 5) == 0)
				{
					double xtemp,ytemp;
					end = 0;
					while (end == 0 && !feof(infile))
					{
						fgets(buf, MAX_CHARS_PER_LINE, infile); // get next line

						if (buf[0] != '*' && buf[1] != '*') // check for commented out line
						{
							// split line into tokens
							token[0] = strtok(buf, DELIMITER); // first token
							for (i = 1; i < 4; i++)
							{
								token[i] = strtok(NULL, DELIMITER); // subsequent tokens
								if (!token[i])
									break; // no more tokens
							}

							if (i < (2 + NUM_DOF))
							{
								printf("\nError in input: not enough data for *bound, point! - Abort\n");
								return -1;
							}

							sscanf(token[0],"%lf",&xtemp); // x coord
							sscanf(token[1],"%lf",&ytemp); // y coord
							for (i = 0; i < NUM_DOF; i++)
							{
								sscanf(token[2+i],"%i",&flag[i]); // is dof fixed
							}

							// find closest node to zero-displacement co-ordinates
							// nd = cmesh.closeNode(inMesh,xtemp,ytemp);
							nd = m_mesh.CloseNode(xtemp,ytemp);

							// first record node - for possible fixed lsf
							// bc_fix[nd] = true;

							// Now fix the specified dof
							nd *= NUM_DOF;
							for (i = 0; i < NUM_DOF; i++)
							{
								fixdof_temp[nd++] = (flag[i] == 0) ? false : true;
							}
						}

						else if(buf[0]=='*' && buf[1]!='*')
						{
							end=1; // stop when next keyword is found
							reread=1; // nned to reread the keyword
						}
					}
				}

				// if type is area
				else if(strncasecmp(token[1], "area", 4)==0)
				{
					end=0;
					while(end==0 && !feof(infile))
					{
						fgets(buf,MAX_CHARS_PER_LINE,infile); // get next line

						double tol = m_mesh.m_tolerance;
						double xmin,xmax,ymin,ymax;

						if(buf[0]!='*' && buf[1]!='*') // check for commented out line
						{
							// split line into tokens
							token[0] = strtok(buf, DELIMITER); // first token
							for(i=1;i<6;i++)
							{
								token[i] = strtok(NULL, DELIMITER); // subsequent tokens
								if (!token[i]) break; // no more tokens
							}

							if(i<(4+NUM_DOF))
							{
								printf("\nError in input: not enough data for *bound, area! - Abort\n");
								return -1;
							}

							sscanf(token[0],"%lf",&xmin); // x min coord
							sscanf(token[1],"%lf",&xmax); // x max coord
							sscanf(token[2],"%lf",&ymin); // y min coord
							sscanf(token[3],"%lf",&ymax); // y max coord
							for(i=0;i<NUM_DOF;i++)
							{
								sscanf(token[4+i],"%i",&flag[i]); // is dof fixed
							}

							for(j=0;j<nn;j++) // For all nodes - determine if the node lies within the fixed rectangular area
							{
							   if( ((cp[j].x-xmax) < tol) && ((xmin-cp[j].x) < tol) ) // If within x bounds
							   {
								   if( ((cp[j].y-ymax) < tol) && ((ymin-cp[j].y) < tol) ) // If also within y bounds
								   {
									   // first record node - for possible fixed lsf
									   // bc_fix[j] = true;

									   // Now fix the specified dof
									   nd = j*NUM_DOF;
									   for(i=0;i<NUM_DOF;i++)
									   {
										   if(flag[i] != 0)
										   {
											   fixdof_temp[nd] = true; // fix dof if flag != 0
										   }
										   nd++; // next dof in fixdof_temp
									   }
								   }
							   }
							}
						}

						else if(buf[0]=='*' && buf[1]!='*')
						{
							end=1; // stop when next keyword is found
							reread=1; // nned to reread the keyword
						}
					}
				}

				// if type is design (bc are design variables)
				else if(strncasecmp(token[1], "design", 6)==0)
				{
					end=0;
					while(end==0 && !feof(infile))
					{
						fgets(buf,MAX_CHARS_PER_LINE,infile); // get next line

						double tol = m_mesh.m_tolerance;
						double h = m_mesh.m_lenEdge;
						int m_elemX = m_mesh.m_elemX;
						int elemY = m_mesh.m_elemY;
						double xmin, xmax, ymin, ymax;

						if (buf[0] != '*' && buf[1] != '*') // check for commented out line
						{
							// split line into tokens
							token[0] = strtok(buf, DELIMITER); // first token
							for (i = 1; i < 4; i++)
							{
								token[i] = strtok(NULL, DELIMITER); // subsequent tokens
								if (!token[i])
									break; // no more tokens
							}

							if (i < 4)
							{
								printf("\nError in input: not enough data for *bound, design! - Abort\n");
								return -1;
							}

							sscanf(token[0],"%lf",&xmin); // x min coord
							sscanf(token[1],"%lf",&xmax); // x max coord
							sscanf(token[2],"%lf",&ymin); // y min coord
							sscanf(token[3],"%lf",&ymax); // y max coord

							// find all elements completely enclosed in area

							// use coords to work out range of elements effected
							xmin -= tol; ymin -= tol; xmax += tol; ymax += tol;
							int ex_min = (int)ceil(xmin/h);  if(ex_min < 0){ex_min=0;}
							int ex_max = (int)floor(xmax/h); if(ex_max > m_elemX){ex_max=m_elemX;}
							int ey_min = (int)ceil(ymin/h);  if(ey_min < 0){ey_min=0;}
							int ey_max = (int)floor(ymax/h); if(ey_max > elemY){ey_max=elemY;}

							// till here, jeehanglee@gmail.com, 150106
							int n,m;
							int num=(ex_max-ex_min)*(ey_max-ey_min);
							if (num > 0)
							{
								if(!m_mesh.m_bDesignableBc)
								{
									m_mesh.m_bDesignableBc=true;
									m_mesh.m_numBc=num;
									m_mesh.m_pBcNums = (int *) malloc(num*sizeof(int));
								}
								else
								{
									m_mesh.m_numBc += num;
									m_mesh.m_pBcNums = (int *) realloc(m_mesh.m_pBcNums,m_mesh.m_numBc*sizeof(int));
								}

								num = m_mesh.m_numBc - num; // start of new memory portion
								for(m=ey_min;m<ey_max;m++)
								{
									for(n=ex_min;n<ex_max;n++)
									{
										m_mesh.m_pBcNums[num++] = m_mesh.m_pNumber[n][m].n; // Element number
									}
								}
							}
						}

						else if(buf[0]=='*' && buf[1]!='*')
						{
							end=1; // stop when next keyword is found
							reread=1; // nned to reread the keyword
						}
					}
				}

				else
				{
					printf("\nError in input: no type for *bound! - Abort\n");
					return -1;
				}

				free(flag);
			}

			// if keyword is fix-lsf
			else if(strncasecmp(lead, "*fix-lsf", 8)==0)
			{
				// check that mesh has defined
				if(meshFound==0)
				{
					printf("\nError in input: *fix-lsf before *mesh! - Abort\n");
					return -1;
				}

				int nn = m_mesh.m_numNodes;
				Coord *cp = m_mesh.m_pNodeCoord;

				// if type is point
				if(strncasecmp(token[1], "point", 5)==0)
				{
					double xtemp,ytemp;
					end=0;
					while(end==0 && !feof(infile))
					{
						fgets(buf,MAX_CHARS_PER_LINE,infile); // get next line

						if(buf[0]!='*' && buf[1]!='*') // check for commented out line
						{
							// split line into tokens
							token[0] = strtok(buf, DELIMITER); // first token
							for (i=1;i<2;i++)
							{
								token[i] = strtok(NULL, DELIMITER); // subsequent tokens
								if (!token[i]) break; // no more tokens
							}

							if(i<2)
							{
								printf("\nError in input: not enough data for *unfix, point! - Abort\n");
								return -1;
							}

							sscanf(token[0],"%lf",&xtemp); // x coord
							sscanf(token[1],"%lf",&ytemp); // y coord

							// find closest node to zero-displacement co-ordinates
							nd = cmesh.closeNode(m_mesh,xtemp,ytemp);

							// first record node - for possible fixed lsf
							m_levelset.m_pFixedLsf[nd] = true;
						}

						else if(buf[0]=='*' && buf[1]!='*')
						{
							end=1; // stop when next keyword is found
							reread=1; // nned to reread the keyword
						}
					}
				}

				// if type is area
				else if(strncasecmp(token[1], "area", 4)==0)
				{
					end=0;
					while(end==0 && !feof(infile))
					{
						fgets(buf,MAX_CHARS_PER_LINE,infile); // get next line

						double tol = m_mesh.m_tolerance;
						double xmin,xmax,ymin,ymax;

						if(buf[0]!='*' && buf[1]!='*') // check for commented out line
						{
							// split line into tokens
							token[0] = strtok(buf, DELIMITER); // first token
							for(i=1;i<4;i++)
							{
								token[i] = strtok(NULL, DELIMITER); // subsequent tokens
								if (!token[i]) break; // no more tokens
							}

							if(i<4)
							{
								printf("\nError in input: not enough data for *fix-lsf, area! - Abort\n");
								return -1;
							}

							sscanf(token[0],"%lf",&xmin); // x min coord
							sscanf(token[1],"%lf",&xmax); // x max coord
							sscanf(token[2],"%lf",&ymin); // y min coord
							sscanf(token[3],"%lf",&ymax); // y max coord

							for(j=0;j<nn;j++) // For all nodes - determine if the node lies within the fixed rectangular area
							{
								if( ((cp[j].x-xmax) < tol) && ((xmin-cp[j].x) < tol) ) // If within x bounds
								{
									if( ((cp[j].y-ymax) < tol) && ((ymin-cp[j].y) < tol) ) // If also within y bounds
									{
										// first record node - for possible fixed lsf
										m_levelset.m_pFixedLsf[j] = true;
									}
								}
							}
						}

						else if(buf[0]=='*' && buf[1]!='*')
						{
							end=1; // stop when next keyword is found
							reread=1; // nned to reread the keyword
						}
					}
				}
				else
				{
					printf("\nError in input: no type for *fix-lsf! - Abort\n");
					return -1;
				}
			}

			// if keyword is load
			else if(strncasecmp(lead, "*load", 5)==0)
			{
				// check that mesh has defined
				if(meshFound==0)
				{
					printf("\nError in input: *load before *mesh! - Abort\n");
					return -1;
				}

				int nn = m_mesh.m_numNodes;
				Coord *cp = m_mesh.m_pNodeCoord;

				// defined number of cases (should only happen for first *load
				if(*numCase==0) // if numCase not been set
				{
					if(token[2]) // if number cases specified
					{
						sscanf(token[2],"%i",numCase);

						if(*numCase < 1 || *numCase > 100)
						{
							printf("\nWarning: num load cases = %i! - using 1 instead",*numCase);
							*numCase=1;
						}
					}
					else // otherwise set as 1
					{
						*numCase=1;
					}

					// set storage for load cases (dof x number cases)
					load = (double *) calloc(NUM_DOF*nn*(*numCase), sizeof(double));
				}

				// if type is point
				if(strncasecmp(token[1], "point", 5)==0)
				{
					int case_cnt = -1; // count number of cases found
					double *mag = (double *) malloc(NUM_DOF * sizeof(double));
					double xtemp,ytemp;
					end=0;
					while(end==0 && !feof(infile) && case_cnt<*numCase)
					{
						fgets(buf,MAX_CHARS_PER_LINE,infile); // get next line

						if(buf[0]!='*' && buf[1]!='*') // check for commented out line
						{
							// split line into tokens
							token[0] = strtok(buf, DELIMITER); // first token
							for(i=1;i<2;i++)
							{
								token[i] = strtok(NULL, DELIMITER); // subsequent tokens
								if (!token[i]) {token[i]="0.0";} // default to zero
							}

							if(case_cnt==-1) // if first line
							{
								sscanf(token[0],"%lf",&xtemp); // x coord
								sscanf(token[1],"%lf",&ytemp); // y coord

								// find closest node to zero-displacement co-ordinates
								nd = cmesh.closeNode(m_mesh,xtemp,ytemp);

								case_cnt++; // up case count (next lines should contain magnitudes)
							}
							else
							{
								for(i=0;i<NUM_DOF;i++)
								{
									sscanf(token[i],"%lf",&mag[i]); // magnitude
								}

								// add loads to the array
								j=(case_cnt*nn*NUM_DOF) + (nd*NUM_DOF);	// point to correct place in load
								for(i=0;i<NUM_DOF;i++)
								{
									load[j+i] += mag[i];	// next dof
								}
								case_cnt++;
							}
						}

						else if(buf[0]=='*' && buf[1]!='*')
						{
							end=1; // stop when next keyword is found
							reread=1; // nned to reread the keyword
						}
					}

					free(mag);
				}
				// if type is area
				else if(strncasecmp(token[1], "dist", 4)==0)
				{
					int case_cnt = -1; // count number of cases found
					double tol = m_mesh.m_tolerance;
					double h = m_mesh.m_lenEdge;
					double min,max,level,xmin,xmax,ymin,ymax;
					int flag;
					double *mag = (double *) calloc(*numCase*NUM_DOF,sizeof(double));
					end=0;
					while(end==0 && !feof(infile) && case_cnt<*numCase)
					{
						fgets(buf,MAX_CHARS_PER_LINE,infile); // get next line

						if(buf[0]!='*' && buf[1]!='*') // check for commented out line
						{
							// split line into tokens
							token[0] = strtok(buf, DELIMITER); // first token
							j = (case_cnt==-1) ? 4 : 2; // total number expected tokens
							for (i=1;i<j;i++)
							{
								token[i] = strtok(NULL, DELIMITER); // subsequent tokens
								if (!token[i]) {token[i]="0.0";} // default to zero
							}

							if(case_cnt==-1) // first line
							{
								sscanf(token[0],"%lf",&min); // min coord
								sscanf(token[1],"%lf",&max); // max coord
								sscanf(token[2],"%lf",&level); // level coord
								sscanf(token[3],"%i",&flag); // 1=x is level, 2=y is level
								case_cnt++; // up case count (next lines should contain magnitudes)

								// check data!!
								if(min > max)
								{
									printf("\nError in input: min > max for level in load.area! - Abort\n");
									return -1;
								}
							}
							// read in distributed load magnitudes
							else
							{
								ind = case_cnt * NUM_DOF; // point to start of next case
								for(i=0;i<NUM_DOF;i++)
								{
									sscanf(token[i],"%lf",&mag[ind]);
									mag[ind++] *= h; // multiply (force / length) by elem edge length
								}
								case_cnt++;
							}
						}

						else if(buf[0]=='*' && buf[1]!='*')
						{
							end=1; // stop when next keyword is found
							reread=1; // need to reread the keyword
						}
					}

					// now all data has been found (or assumed) - apply the distributed load

					// round values to be integer value of h
					level = round(level/h);
					level *= h;
					min = round(min/h);
					min *= h;
					max = round(max/h);
					max *= h;

					// check for bounds
					level = (level < 0.0) ? 0.0 : level;
					min = (min < 0.0) ? 0.0 : min;

					// level refers to an x coord, max & min are y coords
					if(flag==1)
					{
						level = (level > m_mesh.m_maxX) ? m_mesh.m_maxX : level;
						xmin = level - tol;
						xmax = level + tol;
						ymax = (max > m_mesh.m_maxY) ? m_mesh.m_maxY : max;
						ymin = min + tol;
						ymax = ymax - tol;
					}
					// level refers to an y coord, max & min are x coords
					else if(flag==2)
					{
						level = (level > m_mesh.m_maxY) ? m_mesh.m_maxY : level;
						ymin = level - tol;
						ymax = level + tol;
						xmax = (max > m_mesh.m_maxX) ? m_mesh.m_maxX : max;
						xmin = min + tol;
						xmax = xmax - tol;
					}
					else
					{
						printf("\nError in input: wrong flag (%i) for level in load.area! - Abort\n", flag);
						return -1;
					}

					// first find end nodes (required for consistent loading)
					nd = cmesh.closeNode(m_mesh, xmin, ymin); // end 1
					nd2 = cmesh.closeNode(m_mesh, xmax, ymax); // end 2

					for(i=0;i<*numCase;i++)
					{
						// add loads to the array
						ind=(i*nn + nd)*NUM_DOF;	// point to correct place in load (end 1)
						for(j=0;j<NUM_DOF;j++)
						{
							load[ind+j] += mag[(i*NUM_DOF)+j]*0.5;	 // next dof
						}

						ind=(i*nn + nd2)*NUM_DOF;	// point to correct place in load (end 2)
						for(j=0;j<NUM_DOF;j++)
						{
							load[ind+j] += mag[(i*NUM_DOF)+j]*0.5;	 // next dof
						}
					}

					// For all nodes - determine if the node lies on the distributed load line
					for(nd=0;nd<nn;nd++)
					{
						if( ((xmax-cp[nd].x) > 0.0) && ((cp[nd].x-xmin) > 0.0) ) // If within x bounds
						{
							if( ((ymax-cp[nd].y) > 0.0) && ((cp[nd].y-ymin) > 0.0) ) // If also within y bounds
							{
								for(i=0;i<*numCase;i++)
								{
									ind=(i*nn + nd)*NUM_DOF; // point to correct place in load
									for(j=0;j<NUM_DOF;j++)
									{
										load[ind+j] += mag[(i*NUM_DOF)+j]; // next dof
									}
								}
							}
						}
					}
					free(mag);
				}
                // if type is self-weight
                else if(strncasecmp(token[1], "self", 4)==0)
                {
                    if(*sw)
                    {
                        printf("\nError in input: found two *load, self-weight keywords! - Abort\n");
                        return -1;
                    }
                    else
                    {
                        *sw = true; // set self-weight flag to true
                        *acc = (Coord *) malloc(*numCase * sizeof(Coord));

                        int case_cnt = 0; // count number of cases found
                        double xtemp,ytemp;
                        end=0;
                        while(end==0 && !feof(infile) && case_cnt<*numCase)
                        {
                            fgets(buf,MAX_CHARS_PER_LINE,infile); // get next line

                            if(buf[0]!='*' && buf[1]!='*') // check for commented out line
                            {
                                // split line into tokens
                                token[0] = strtok(buf, DELIMITER); // first token
                                for(i=1;i<2;i++)
                                {
                                    token[i] = strtok(NULL, DELIMITER); // subsequent tokens
                                    if (!token[i]) {token[i]="0.0";} // default to zero
                                }

                                sscanf(token[0],"%lf",&xtemp); // x acc
                                sscanf(token[1],"%lf",&ytemp); // y acc

                                acc[0][case_cnt].x = xtemp;
                                acc[0][case_cnt].y = ytemp;

                                case_cnt++; // up case count
                            }

                            else if(buf[0]=='*' && buf[1]!='*')
                            {
                                end=1; // stop when next keyword is found
                                reread=1; // nned to reread the keyword
                            }
                        }
                    }
                }
				else
				{
					printf("\nError in input: no type for *bound! - Abort\n");
					return -1;
				}
			}

			// if keyword is objective
			else if(strncasecmp(lead, "*obj", 4)==0)
			{
				// second token should be a keyword relating to objective type
				// if type is compliance
				if(!token[1]){token[1] = "";}
				if(strncasecmp(token[1], "comp", 4)==0) { m_problem->obj=1; }
				else if(strncasecmp(token[1], "vol", 3)==0) { m_problem->obj=2; }
				else if(strncasecmp(token[1], "freq", 4)==0) { m_problem->obj=3; }
				//else if(strncasecmp(token[1], "disp", 4)==0) { m_problem->obj=4; }
				else if(strncasecmp(token[1], "mech", 4)==0) { m_problem->obj=5; }

				else
				{
					printf("\nError in input: could not determine objective! - Abort\n");
					return -1;
				}
			}

			// if keyword is constraint
			else if(strncasecmp(lead, "*const", 6)==0)
			{
				// second token should be the number of constraints
				if(token[1]) // if number cases specified
				{
					sscanf(token[1],"%i",&m_problem->num);
				}

				int numC = m_problem->num;

				if(numC > 0)
				{
					// create array to store constraint data
					cnst *cons = (cnst *) malloc(numC * sizeof(cnst));
					m_problem->con = cons;
					int count = 0;

					// read in keyword and associated dataline for each constraint
					end=0;
					while(end==0 && !feof(infile) && count<numC)
					{
						fgets(buf,MAX_CHARS_PER_LINE,infile); // get next line

						if(buf[0]!='*' && buf[1]!='*') // check for commented out line
						{
							// split line into tokens
							token[0] = strtok(buf, DELIMITER); // first token
							for(i=1;i<11;i++)
							{
								token[i] = strtok(NULL, DELIMITER); // subsequent tokens
								if (!token[i]) {token[i]="0.0";} // default to zero
							}

							// volume constraint
							if(strncasecmp(token[0], "vol", 3)==0)
							{
								cons[count].type = 1; // constraint type
								if(i>2)
								{
									//sscanf(token[1],"%i",&cons[count].sign);
									cons[count].sign = -1; // less than
									sscanf(token[2],"%lf",&cons[count++].data[0]);
								}
								else
								{
									printf("\nError in input: not enough data for volume constraint! - Abort\n");
									return -1;
								}
							}

							// mass constraint
							else if(strncasecmp(token[0], "mass", 4)==0)
							{
								cons[count].type = 2; // constraint type
								if(i>2)
								{
									//sscanf(token[1],"%i",&cons[count].sign);
									cons[count].sign = -1; // less than
									sscanf(token[2],"%lf",&cons[count++].data[0]);
								}
								else
								{
									printf("\nError in input: not enough data for mass constraint! - Abort");
									return -1;
								}
							}

							// compliance constraint
							else if(strncasecmp(token[0], "comp", 4)==0)
							{
								cons[count].type = 3; // constraint type
								if(i>2)
								{
									//sscanf(token[1],"%i",&cons[count].sign);
									cons[count].sign = -1; // less than
									sscanf(token[2],"%lf",&cons[count++].data[0]);
								}
								else
								{
									printf("\nError in input: not enough data for mass constraint! - Abort");
									return -1;
								}
							}

							/*
							// frequency constraint
							else if(strncasecmp(token[0], "freq", 4)==0)
							{
								cons[count].type = 4; // constraint type
								if(i>3)
								{
									sscanf(token[1],"%i",&cons[count].sign);
									sscanf(token[2],"%lf",&cons[count].data[0]);
									sscanf(token[3],"%lf",&cons[count++].data[1]);
								}
								else
								{
									printf("\nError in input: not enough data for frequency constraint! - Abort");
									return -1;
								}
							}
							*/

							// displacement constraint (on a node)
							else if(strncasecmp(token[0], "disp", 4)==0)
							{
								cons[count].type = 5; // constraint type
								if(i>5)
								{
									Coord p;
									int flag;
									//sscanf(token[1],"%i",&cons[count].sign);
									cons[count].sign = 0; // equality
									sscanf(token[2],"%lf",&p.x); // x coord
									sscanf(token[3],"%lf",&p.y); // y coord
									sscanf(token[4],"%i",&flag); // 1 = x, 2 = y direction
									sscanf(token[5],"%lf",&cons[count].data[0]); // magnitude
									sscanf(token[6],"%lf",&cons[count].data[2]); // load case

									if(flag < 1 || flag > 2)
									{
										printf("\nError in input: wrong direction flag for displacement constraint! - Abort");
										return -1;
									}

									// data[1] will be dof number
									nd = cmesh.closeNode(m_mesh, p.x, p.y); // closest node
									nd*=NUM_DOF; // x dof
									nd+=flag-1;  // (possible) y dof
									cons[count].data[1] = (double)nd;

									count++; // increase constraint count
								}
								else
								{
									printf("\nError in input: not enough data for displacement constraint! - Abort");
									return -1;
								}
							}

							// constraint on ratio of 1st and 2nd eigen-frequencies
							else if(strncasecmp(token[0], "eig_ratio", 9)==0)
							{
								cons[count].type = 6; // constraint type
								if(i>2)
								{
									//sscanf(token[1],"%i",&cons[count].sign);
									cons[count].sign = 1; // greater than
									sscanf(token[2],"%lf",&cons[count++].data[0]);
								}
								else
								{
									printf("\nError in input: not enough data for eig freq ratio constraint! - Abort");
									return -1;
								}
							}

							// input displacement constraint (for compliant mechanisms - obj = mesh)
							else if(strncasecmp(token[0], "in_disp", 7)==0)
							{
								// check that mesh has defined
								if(meshFound==0)
								{
									printf("\nError in input: displacement constrant defined before *mesh! - Abort\n");
									return -1;
								}

								cons[count].type = 7; // constraint type
								if(i>2)
								{
									//sscanf(token[1],"%i",&cons[count].sign);
									cons[count].sign = -1; // less than
									sscanf(token[2],"%lf",&cons[count].data[0]); // input disp value
									if(i>3){sscanf(token[3],"%lf",&cons[count].data[3]);} // spring stiffness
									else{cons[count].data[3] = 0.0;}
								}
								else
								{
									printf("\nError in input: not enough data for mechanism input disp constraint! - Abort");
									return -1;
								}
								count++; // increase constraint count
							}

							// optional output compliance constraint (for compliant mechanisms - obj = mesh)
							/*else if(strncasecmp(token[0], "out_comp", 8)==0)
							{
								cons[count].type = 8; // constraint type
								if(i>2)
								{
									sscanf(token[1],"%i",&cons[count].sign);
									sscanf(token[2],"%lf",&cons[count++].data[0]);
								}
								else
								{
									printf("\nError in input: not enough data for output compliance constraint! - Abort");
									return -1;
								}
							}*/

							// bc cost constraint
							else if(strncasecmp(token[0], "bc_cost", 7)==0)
							{
								cons[count].type = 10; // constraint type
								if(i>2)
								{
									sscanf(token[1],"%i",&cons[count].sign);
									sscanf(token[2],"%lf",&cons[count++].data[0]); // fraction of number of constraints
								}
								else
								{
									printf("\nError in input: not enough data for bc cost constraint! - Abort\n");
									return -1;
								}
							}

							else{ printf("\nWarning in input: unknown constraint type - %s",token[0]); }
						}

						else if(buf[0]=='*' && buf[1]!='*')
						{
							if(count < numC)
							{
								printf("\nError in input: not enough data for constraints! - Abort\n");
								return -1;
							}
							end=1; // stop when next keyword is found
							reread=1; // need to reread the keyword
						}
					}

				}
			}

			// if keyword is control
			else if(strncasecmp(lead, "*control", 8)==0)
			{
				end=0;
				while(end==0 && !feof(infile))
				{
					fgets(buf,MAX_CHARS_PER_LINE,infile); // get next line

					if(buf[0]!='*' && buf[1]!='*') // check for commented out line
					{
						// split line into tokens
						token[0] = strtok(buf, DELIMITER); // first token
						for(i=1;i<7;i++)
						{
							token[i] = strtok(NULL, DELIMITER); // subsequent tokens
							//if (!token[i]) {break;} // default to zero
						}
						for(j=i;j<6;j++) { token[j] = NULL; } // no more tokens

						if(token[0]) {
							sscanf(token[0],"%i",&m_control->maxItt);
						}
						if(token[1]) {
							sscanf(token[1],"%i",&m_control->pinfo);
						}
						if(token[3]) {
							sscanf(token[3],"%lf",&m_control->gm);
						}
						if(token[4]) {
							sscanf(token[4],"%lf",&m_control->lband);
						}
						if(token[5]) {
							sscanf(token[5],"%lf",&m_control->aMin);
						}
						if(token[6]) {
							sscanf(token[6],"%lf",&m_control->mMin);
						}
						end=1; // found data line
					}
					else if(buf[0]=='*' && buf[1]!='*')
					{
						end=1; // stop when next keyword is found
						reread=1; // nned to reread the keyword
					}
				}
			}

			// if keyword is hole
			else if(strncasecmp(lead, "*hole", 5)==0)
			{
				if(!token[1]){token[1] = "";}

				// if type is circle
				if(strncasecmp(token[1], "circ", 4)==0)
				{
					end=0;
					while(end==0 && !feof(infile))
					{
						fgets(buf,MAX_CHARS_PER_LINE,infile); // get next line

						if(buf[0]!='*' && buf[1]!='*') // check for commented out line
						{
							// split line into tokens
							token[0] = strtok(buf, DELIMITER); // first token
							for (i=1;i<3;i++)
							{
								token[i] = strtok(NULL, DELIMITER); // subsequent tokens
								if (!token[i]) break; // no more tokens
							}

							if(i<3)
							{
								printf("\nError in input: not enough data for *hole, circle! - Abort\n");
								return -1;
							}

							sscanf(token[0],"%lf",&holes[NumHole].x); // x coord
							sscanf(token[1],"%lf",&holes[NumHole].y); // y coord
							sscanf(token[2],"%lf",&holes[NumHole].r); // is x dof fixed
							NumHole++;

							if(NumHole == MAX_HOLES)
							{
								printf("\nWarning: Maximum number of circular initial holes reached - %i",MAX_HOLES);
								end=1;
							}
						}

						else if(buf[0]=='*' && buf[1]!='*')
						{
							end=1; // stop when next keyword is found
							reread=1; // need to reread the keyword
						}
					}
				}

				// if type is rectangle
				else if(strncasecmp(token[1], "rect", 4)==0)
				{
					end=0;
					while(end==0 && !feof(infile))
					{
						fgets(buf,MAX_CHARS_PER_LINE,infile); // get next line

						if(buf[0]!='*' && buf[1]!='*') // check for commented out line
						{
							// split line into tokens
							token[0] = strtok(buf, DELIMITER); // first token
							for (i=1;i<4;i++)
							{
								token[i] = strtok(NULL, DELIMITER); // subsequent tokens
								if (!token[i]) break; // no more tokens
							}

							if(i<4)
							{
								printf("\nError in input: not enough data for *hole, rectangle! - Abort\n");
								return -1;
							}

							j=NumRect*2; // indicator
							sscanf(token[0],"%lf",&Rect[j].x);	 // x min
							sscanf(token[1],"%lf",&Rect[j+1].x); // x max
							sscanf(token[2],"%lf",&Rect[j].y);   // y min
							sscanf(token[3],"%lf",&Rect[j+1].y); // y max
							NumRect++;

							if(NumRect == MAX_HOLES)
							{
								printf("\nWarning: Maximum number of rectangular initial holes reached - %i",MAX_HOLES);
								end=1;
							}
						}

						else if(buf[0]=='*' && buf[1]!='*')
						{
							end=1; // stop when next keyword is found
							reread=1; // need to reread the keyword
						}
					}
				}

				else
				{
					printf("\nError in input: no type for *hole! - Abort\n");
					return -1;
				}
			}

			else
			{
				printf("\nWarning unknown keyword in input file: %s",lead);
			}
		}
	}

	fclose(infile);

	// ------ check data is complete ------ //

	// check mesh is defined
	if(meshFound==0)
	{
		printf("\nError in input: mesh definition not found! - Abort\n");
		return -1;
	}

	// check material defined
	if(matCnt == 0)
	{
		printf("\nError in input: material definition not found! - Abort\n");
		return -1;
	}
	else
	{
		*numMat = matCnt;
		nd = m_mesh.NumElem;
		for(i=0;i<nd;i++)
		{
			if(m_mesh.mat_type[i] >= matCnt)
			{
				printf("\nError in input: element %i has undefined material %i! - Abort\n",i,m_mesh.mat_type[i]);
				return -1;
			}
		}
	}

	// use fixdof to create a dof map (saves time later)
	// also reduce load array
	int allDof = NUM_DOF * m_mesh.m_numNodes; // total dof
	*map = (int *) malloc(allDof * sizeof(int));
	nd = 0; // count free dof
    for(i=0;i<allDof;i++)
    {
        if(fixdof_temp[i])
		{
			(*map)[i] = -1; // fixed dof
		}
		else
		{
			(*map)[i] = nd;  // free dof
			nd++; // move to next dof
		}
    }
	free(fixdof_temp);

	*freeDof = nd; // total free dof

	// check for sufficient boundary conditions
	// at least 3 fixed dofs in 2D (2 translation & 1 rotation)
	if( (allDof-nd) < 3 )
	{
		printf("\nInsufficient fixed dof! - Abort\n");
		return -1;
	}

	// check for load definition (if obj or constraints require loads)
	end=0;
	if(m_problem->obj==1){end=1;} // if compliance objective
	else if(m_problem->obj==5){end=1;} // else if compliant mechanism
	else{
		for(i=0;i<m_problem->num;i++){
			if(m_problem->con[i].type==3 || m_problem->con[i].type==5) // if compliance or disp constraint
			{end=1; break;}
		}
	}
	if(end==1)
	{
		bool load_found = false; // check for at least 1 non-zero load
		double ltemp;
		if(*numCase>0)
		{
			// create reduced load array from map
			*load_in = (double *) calloc((*freeDof)*(*numCase), sizeof(double));

			for(i=0;i<allDof;i++) // for all dof
			{
				ind = (*map)[i];
				if(ind > -1) // if dof not fixed
				{
					for(j=0;j<*numCase;j++) // for all load cases
					{
						nd = allDof * j; // point to loaction in load
						nd2 = *freeDof * j; // point to loaction in load_in
						ltemp = load[nd+i];
						if(!load_found && fabs(ltemp) > 0.0){load_found = true;}
						(*load_in)[nd2+ind] = ltemp; // copy load accross
					}
				}
			}
			free(load);
		}
		else
		{
			printf("\nProblem requires loading, but number of load cases = 0! - Abort\n");
			return -1;
		}
		if(!load_found)
		{
			printf("\nProblem requires loading, but load vector is 0! - Abort\n");
			return -1;
		}
	}

	// consolodate the lumped mass vector
	if(lm_temp)
	{
		// set initial memory
		lump_mass->ne = allDof;
		lump_mass->irn=0; lump_mass->jcn=0; lump_mass->A=0;
		fem.set_sp_mat(lump_mass);

		j=0;
		for(i=0;i<allDof;i++) // for all dof
		{
			if(lm_temp[i] > 0.0) // if lumped mass at dof
			{
				lump_mass->irn[j] = i;
				lump_mass->jcn[j] = i;
				lump_mass->A[j++] = lm_temp[i];
			}
		}
		free(lm_temp);

		// reallocate memory
		lump_mass->ne = j;
		lump_mass->irn=(int *) realloc(lump_mass->irn,j*sizeof(int));
		lump_mass->jcn=(int *) realloc(lump_mass->jcn,j*sizeof(int));
		lump_mass->A=(double *) realloc(lump_mass->A,j*sizeof(double));
	}
	else {
		lump_mass->ne=0;
	}

	// check that obj & constraints are defined (and make sense)
	if(m_problem->obj == 1)
	{
		nd = m_problem->num; // number of constraints
		for(i=0;i<nd;i++)
		{
			// if constraint is volume fraction, change to real volume
			if(m_problem->con[i].type == 1)
			{
				m_problem->con[i].data[0] *= m_mesh.m_lenEdge * m_mesh.m_lenEdge * m_mesh.NumElem;
			}
			else if(m_problem->con[i].type == 5) { } // more data checks !!!
			else if(m_problem->con[i].type == 2) { } // mass constraint
			else if(m_problem->con[i].type == 10)
			{
				if(!m_mesh.m_bDesignableBc)
				{
					printf("\nBC cost constraint defoned without *bound, design ! - Abort\n");
					return -1;
				}
				m_problem->con[i].data[0] *= (double)m_mesh.m_numBc; // multiple by number of elem with designable bcs
			}
			else
			{
				printf("\nCompliance objective has unidentified or wrong constraint ! - Abort\n");
				return -1;
			}
		}
	}
	else if(m_problem->obj == 2)
	{
		nd = m_problem->num; // number of constraints
		for(i=0;i<nd;i++)
		{
			// if constraint is compliance
			if(m_problem->con[i].type == 3) { }
			else if(m_problem->con[i].type == 5) { } // more data checks !!!
			else
			{
				printf("\nVolume objective has unidentified or wrong constraint ! - Abort\n");
				return -1;
			}
		}
	}
	else if(m_problem->obj == 3)
	{
		//if(*numCase==0){*numCase=1;} // ensure this is at least 1
		nd = m_problem->num; // number of constraints
		for(i=0;i<nd;i++)
		{
			// if constraint is volume fraction, change to real volume
			if(m_problem->con[i].type == 1)
			{
				m_problem->con[i].data[0] *= m_mesh.m_lenEdge * m_mesh.m_lenEdge * m_mesh.NumElem;
			}
			else if(m_problem->con[i].type == 2) { } // mass constraint
			else if(m_problem->con[i].type == 6) { } // eig ratio constraint
			else
			{
				printf("\n1st Eigenvalue objective has unidentified or wrong constraint ! - Abort\n");
				return -1;
			}
		}
        if(*sw)
        {
            printf("\nSelf-weight specified for eigenvalue problem! - Ignoring\n");
            *sw = false;
        }
	}
	else if(m_problem->obj == 5)
	{
		nd = m_problem->num; // number of constraints
		for(i=0;i<nd;i++)
		{
			// if constraint is volume fraction, change to real volume
			if(m_problem->con[i].type == 1)
			{
				m_problem->con[i].data[0] *= m_mesh.m_lenEdge * m_mesh.m_lenEdge * m_mesh.NumElem;
			}
			// if constraint is input displacement or compliance
			else if(m_problem->con[i].type == 7){ }
			else if(m_problem->con[i].type == 8){ }

			else
			{
				printf("\nCompliant mechnaism has unidentified or wrong constraint ! - Abort\n");
				return -1;
			}
		}
        if(*sw)
        {
            printf("\nSelf-weight specified for compliant mechanism problem! - Ignoring\n");
            *sw = false;
        }
	}

	else
	{
		printf("\nObjective not specified! - Abort\n");
		return -1;
	}

	// check bars (and that problem is min compliance s.t. mass constraint)
	if(m_mesh.m_bars)
	{
		if(m_problem->obj != 1)
		{
			printf("\nWarning! Bars specified, but objective is not compliance - ignoring bars\n");
			m_mesh.m_bars = false;
		}
		else if(m_problem->num > 1 || m_problem->con[0].type != 2)
		{
			printf("\nWarning! Bars specified, but constraint is not mass - ignoring bars\n");
			m_mesh.m_bars = false;
		}
	}

	if(m_mesh.m_bars)
	{
		if(m_mesh.m_barMin < 0.0 || m_mesh.m_barMax < m_mesh.m_barMin)
		{
			printf("\nError in bar area side constraints ! - Abort\n");
			return -1;
		}

		// initialize bar areas as average of max & min
		double ftemp = 0.5*(m_mesh.m_barMax + m_mesh.m_barMin);

		for(j=0;j<m_mesh.m_numBars;j++)
		{
			m_mesh.m_pBarAreas[j] = ftemp;
		}
	}

	if(m_mesh.m_bDesignableBc)
	{
		if(m_mesh.m_bars)
		{
			printf("\nError *bc, design and *bars are mutually exclusive! - Abort\n");
			return -1;
		}

		nd = m_problem->num; // number of constraints
		for(i=0;i<nd;i++)
		{
			if(m_problem->con[i].type == 10){break;}
		}
		if(i==nd)
		{
			printf("\nError no cost constraint for designable BC ! - Abort\n");
			return -1;
		}

		// sort bc numbers in asscening order
		qsort(m_mesh.m_pBcNums, m_mesh.m_numBc, sizeof(int), icmpfunc);

		// initialize designable bc variables
		m_mesh.K_bc = (double *) malloc(m_mesh.m_numBc * sizeof(double));
		m_mesh.K0_bc = m_material[0].e * 1.0e6; // max stiffness
		double ftemp = m_problem->con[i].data[0] / (double)m_mesh.m_numBc; // initial design is feasible
		for(j=0;j<m_mesh.m_numBc;j++)
		{
			//printf("\n%i",m_mesh.m_pBcNums[j]);
			m_mesh.K_bc[j] = ftemp;
		}
	}

	if(m_mesh.des_mat)
	{
		if(m_problem->obj != 3 && m_problem->obj != 1)
		{
			printf("\nWarning! Designable material specified, but objective is not compliance or frequency - Aborting\n");
			return -1;
		}

		m_mesh.m_pMaterialVars = (double *) malloc(m_mesh.NumDesMat * sizeof(double));
		for(j=0;j<m_mesh.NumDesMat;j++)
		{
			m_mesh.m_pMaterialVars[j] = 0.5; // start with 50/50 mix
		}

        // FD check
        //m_mesh.m_pMaterialVars[17700] = 0.49;

		// redefine material for all elems with des mat to the first material
		for(j=0;j<m_mesh.NumDesMat;j++)
		{
			nd = m_mesh.mat_elems[j];
			m_mesh.mat_type[nd] = m_mesh.mat1;
		}
	}

	// check controls
	if(m_control->maxItt < 1)
	{
		printf("\nWarning! Specified max number of iterations < 1 - Using default 200");
		m_control->maxItt = 200;
	}
	if(m_control->pinfo > 3 || m_control->pinfo < 1)
	{
		printf("\nWarning! info printout flag out of range! - Using default 2");
		m_control->pinfo = 2;
	}
	if(m_control->gm < 1.0e-5 || m_control->gm > 1.0e-2)
	{
		printf("\nWarning! covergence criterion out of range! - Using default 1e-3");
		m_control->gm = 1.0e-3;
	}
	if(m_control->lband < 2.0)
	{
		printf("\nWarning! narrow band too small! - Using default 6h");
		m_control->lband = 6.0;
	}
	// multiply lband by h
	m_control->lband = m_control->lband * m_mesh.m_lenEdge;

	if(m_control->aMin < 0.0 || m_control->aMin > 1.0e-3)
	{
		printf("\nWarning! min area ratio out of range! - Using default 1.0e-6");
		m_control->aMin = 1.0e-6;
	}

	if(m_control->mMin < 0.0 || m_control->mMin > 1.0e-3)
	{
		printf("\nWarning! min area ratio (mass) out of range! - Using default 0.0");
		m_control->mMin = 0.0;
	}

	// check for initial holes (warn if none found)
	if(NumHole < 1 && NumRect < 1)
	{
		printf("\nWarning! Structure has no initial holes!");
	}

	// populate the level set struct
	m_levelset.m_numNode = m_mesh.m_numNodes;
	m_levelset.m_pNodalLsf =  (double *) malloc(m_mesh.m_numNodes * sizeof(double));
	m_levelset.active = (bool *) malloc(m_mesh.m_numNodes * sizeof(bool));
	m_levelset.numMine = 0;
	m_levelset.mine = (int *) malloc(m_mesh.m_numNodes * sizeof(int));

	// initalize signed distance function
	m_levelset.initialLsf(m_mesh, m_levelset, NumHole, holes, NumRect, Rect, m_control->lband);

	// clean up memory
	free(holes);
	free(Rect);

	return 0;
}

// function to read new-style input file
