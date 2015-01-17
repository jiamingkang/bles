/*
 * CExMinimiseCompliance.h
 *
 *  Created on: Dec 1, 2014
 *      Author: jeehang, Khalid Ismail
 */

#ifndef CEXMINIMISECOMPLIANCE_H_
#define CEXMINIMISECOMPLIANCE_H_

#include <string>
#include "CHakInput.h"
#include "CHakOutput.h"

#include "CommonTypes.h"
#include "CHakFiniteElement.h"
#include "CHakMathUtility.h"
#include "CHakSolver.h"
#include "CHakSensitivity.h"
#include "CHakLevelSet.h"
#include "CHakMesh.h"
#include "CHakBoundary.h"

//using std;




class CExMinimiseCompliance {
public:
	CExMinimiseCompliance();
	virtual ~CExMinimiseCompliance();

 
//
// Interfaces
//
//knai20@bath.ac.uk: OOD implementations
public:
    int InitialiseOOP(char *path);
	// Finite Element Anaysis
	// returns -1 when exceptions occurs thus process fails.
	int AnalyseOOP(int itt);

	// Sensitivity
	int SensitivityOOP();

	// Optimisation
	int OptimiseOOP();

    void OutputOOP();

//knai20@bath.ac.uk: Non OOD implementations
public:
	// Read the input file & assign initial values
	// returns -1 when starting process fails.
	int StartProcess(std::string inputfile);
    
    int Initialise(char *path);
    
	// Finite Element Anaysis
	// returns -1 when exceptions occurs thus process fails.
	int Analyse(int itt);

	// Sensitivity 
	int Sensitivity();

	// Optimisation
	int Optimise();
    
    void Output();
	
	// main, dummy function for the compatibility
	void Solve(char* arg, char* filename);

//
// Utilities
//
protected:
	void _SetFilePath(char *path);

//
// Attributes
//
public:
    int i,j,i2,j2,k,temp,temp2;	// incrementors
    double ftemp;
    char filename[100];
    
    // inital data arrays & structs

    // Structures
    //mesh inMesh;  // struct to hold mesh data //Khalid: remove when OOD
    //levSet levelset; // struct to hold level set info
    //prob lsprob;  // struct to hold problem defintion
    //ctrl control; // struct to hold control data
    //isoMat inMat[5]; // isotropic material - maximum of 5 different materials
    //sp_mat lump_mass; // lumped mass matrix
    //boundary Bndry;	// boundry discretization
    //sp_mat Kg, Mg; // global stiffness and mass matrices
    
    Coord *gCoord;
    Coord *acc;   // acceleration vector for self-weight loading

    int numMat; // numberof materials
    int numCase;  // number load cases
    double *load; // load vector (rhs)
    bool sw = false;      // self-weight loading flag
    int *fixDof;  // fixed dof (turn into map)
    int freeDof;  // number of free dof

    double AreaElem;
    
    double **KE;
    double **ME;
    
    int itt = 0; // Initalise number of itterations to 0
    int itt0 = 0;
    double oMax, oMin, obj_val;	 // varibale to compute convergence criterion
    int ReCount = 0;  // initialize reinitialization count
    double cfl = 0.5; // set time step modifier (for CFL condition)
    double u11, u12, u21, u22, p1, p2; // displacement & load values - used for compliant mechanism design
    double fact[6];
    
    // Variables and arrays to store additional mesh data
    int Ntot;		 // total nodes, including auxillary ones
    double *alpha;  // Array to store element areas
    
    // Arrays to store node data related to the optimisation
    double *Nsens;  // pointer for node (+ aux node) sensitivity array
    double *vol_sens=0; // volume sensitivity array (will be all 1's)
    double *mass_sens=0; // mass sensitvitiy array
    double *zero_sens=0; // when lsf does not influence a constraint
    double **sens_ptr; // pointer to senstivity arrays for each fucntion (obj + constraints)
    double *Vnorm;	// pointer for node (+ aux node) normal velocity array
    double *Grad;   // pointer for lsf gradient info (using upwind scheme)
    int *active;  //array for active constraints
    
    // boundary integral variables
    double *Lbound;   // array to store boundary intergral coefficents
    int *Lbound_nums; // global node nums corresponding to Lbound entries
    int num_Lbound;   // length of Lbound & Lbound_nums
    double *lambda;
    
    // stuff for Ku = f solve using MA57 solver (& eigenvalue solver)
    int num_eig; // number of eigenvalues to compute
    int comp_eig; // flag to decide if eigenvalues are to be computed
    double *eig_vals, *eig_vecs; // arrays for eigenvalues and vectors
    double *save_freq; // used to store freq values for output
    int order; // full matrix order
    int numEnt; // max number of index entries (symmetry!)
    
    // bar additonal design variables
    double *bar_sens=0; // bar senstivities
    double *bar_step=0; // bar update step
    double *bar_max=0, *bar_min=0; // move limits on bar areas
    
    // designable bc addtional variables
    double *bc_sens=0; // bc senstivities
    double *bc_step=0; // bc update step
    double *bc_max=0, *bc_min=0; // move limits on bc variables
    
    // designable material varibles
    double *mat_sens=0; // material sensitivites
    double *mat_step=0; // bc update step
    double *mat_max=0, *mat_min=0; // move limits on bc variables
    bool mat_opt;

    double *disp; // displacement array
    double *load_sw; // load vector, including self_weight load vector
    double *adjont; // adjoint load (then disp) vector - for non-self adjoint problems
    
    int dispLen; // length of disp array (inc fixed dof)
    int loadLen; // length of load array (excld fixed dof)
    
    int num_adj = 0; // number of adjoint cases to solve (objective, then in order of constraint list)

    int num_sens = 0; // number of senstivity calculations

    double **prim_ptr; // vector, length = numCase (+ num_eig)
    double **dual_ptr; // matrix row num = load case num, col num = dual state num
    // or row = eigenvalue num, col = dual state num
    double *wgt; // weights for load cases
    
    // Arrays to store objective & constraint data for each iteration
    double *obj;
    double *cnstr;
    double *delcon; // array to store distance from constraint (constraint violation)
    double *relax; // relaxation factor
    double *pred; // predicted obj & constraint change
    double *pred_temp; // temp to feed SLPsubSol

    //OOD implementation knai20@bath.ac.uk
private:
	// Input file path
	char *m_pFileIn;

	// Output file path
	char m_fileOut[128];

	// Read input file and assign initial values specified in the file 
	CHakInput m_input;

	// handle output files
	CHakOutput m_output;

    CHakMesh m_mesh;

    CHakMathUtility m_mathUtility;

    CHakFiniteElement m_fem;

    CHakLevelSet m_levelset;

    CHakMaterial m_material;

    CHakSolver m_solver;

    CHakSensitivity m_sensitivity;

    CHakBoundary m_boundary;
};

#endif /* CEXMINIMISECOMPLIANCE_H_ */
