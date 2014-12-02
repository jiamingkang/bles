/*
 * CExMinimiseCompliance.cpp
 *
 *  Created on: Dec 1, 2014
 *      Author: jeehang
 */

#include "CommonTypes.h"
#include "CFiniteElement.h"
#include "CMathUtility.h"
#include "CSolver.h"
#include "CSensitivity.h"
#include "CLevelSet.h"
#include "CInput.h"
#include "COutput.h"
#include "CMesh.h"
#include "CBoundary.h"

#include "CExMinimiseCompliance.h"

CExMinimiseCompliance::CExMinimiseCompliance() {
	// TODO Auto-generated constructor stub

}

CExMinimiseCompliance::~CExMinimiseCompliance() {
	// TODO Auto-generated destructor stub
}


//
// implementation
//

void CExMinimiseCompliance::Solve(char* arg, char* filename)
{
	int i,j,i2,j2,k,temp,temp2;	// incrementors
	double ftemp;

	// inital data arrays & structs
	mesh inMesh;  // struct to hold mesh data
	int numMat; // numberof materials
	isoMat inMat[5]; // isotropic material - maximum of 5 different materials
	levSet levelset; // struct to hold level set info
	prob lsprob;  // struct to hold problem defintion
	ctrl control; // struct to hold control data
	int numCase;  // number load cases
	double *load; // load vector (rhs)
	bool sw = false;      // self-weight loading flag
	Coord *acc;   // acceleration vector for self-weight loading
	int *fixDof;  // fixed dof (turn into map)
	int freeDof;  // number of free dof
	sp_mat lump_mass; // lumped mass matrix

	// read the input file
	CInput cinput;
	temp = cinput.read_input(arg, &inMesh, &numMat, inMat, &levelset, &lsprob, &control, &fixDof,
			&numCase, &load, &freeDof, &lump_mass, &sw, &acc);
	if(temp==-1)//{return -1;} // exit on error
	{
		printf("Failure in reading input file....\n");
		return;
	}

	COutput coutput;
	if(control.pinfo==3)
	{
		coutput.OutNumber(&inMesh, filename);
		coutput.OutNodeCoord(&inMesh, filename);
	}

	// compute gauss point coords (used in sensitivity smoothing)
	CMathUtility cmu;
	Coord *gCoord = malloc(4 * inMesh.NumElem * sizeof(Coord));
	cmu.Gauss_Coord(&inMesh, gCoord);

	// calculate IN element stiffness matrix (& mass matrix)
	double AreaElem = inMesh.h * inMesh.h; // Area of an element
	double **KE = malloc(numMat*sizeof(double*));
	double **ME = malloc(numMat*sizeof(double*));
	for(i=0;i<numMat;i++)
	{
		KE[i] = malloc(KE_SIZE*sizeof(double));
		ME[i] = malloc(KE_SIZE*sizeof(double));
		// use 1.0 for thickness, as E & rho already multipled by thickness
		KEMatrix(KE[i], &inMat[i], inMesh.t);	// Element Stiffness Matrix for an IN element
		MEMatrix(ME[i], inMat[i].rho, AreaElem, inMesh.t); // Element Mass Matrix for an IN element
	}

	/*------------------------------------------------------------------/
	/																	/
	/		Optimise structure boundary by the level set method			/
	/																	/
	/------------------------------------------------------------------*/

	int itt = 0; // Initalise number of itterations to 0
	int itt0 = 0;
	double oMax, oMin, obj_val;	 // varibale to compute convergence criterion
	int ReCount = 0;  // initialize reinitialization count
	double cfl = 0.5; // set time step modifier (for CFL condition)
	double u11, u12, u21, u22, p1, p2; // displacement & load values - used for compliant mechanism design
	double fact[6];

	// Variables and arrays to store additional mesh data
	int Ntot;		 // total nodes, including auxillary ones
	double *alpha = malloc(inMesh.NumElem * sizeof(double)); // Array to store element areas

	// Arrays to store node data related to the optimisation
	double *Nsens;  // pointer for node (+ aux node) sensitivity array
	double *vol_sens=0; // volume sensitivity array (will be all 1's)
	double *mass_sens=0; // mass sensitvitiy array
	double *zero_sens=0; // when lsf does not influence a constraint
	double **sens_ptr = malloc( (1+lsprob.num) * sizeof(double *) ); // pointer to senstivity arrays for each fucntion (obj + constraints)
	double *Vnorm;	// pointer for node (+ aux node) normal velocity array
	double *Grad;   // pointer for lsf gradient info (using upwind scheme)
	int *active = malloc(lsprob.num * sizeof(int)); // array for active constraints

	// boundary integral variables
	boundary Bndry;	// boundry discretization
	// Initialise length of AuxNodes & Boundary segment arrays - generous initial size
	Bndry.AuxNodes = malloc(inMesh.NumElem * sizeof(Coord));
	Bndry.Bound = malloc(inMesh.NumElem * sizeof(Bseg));
	Bndry.BsegLen = malloc(sizeof(double)); // initial memory allocation
	Bndry.Bwgt = malloc(sizeof(double)); // initial memory allocation
	Bndry.na_conn = malloc(sizeof(int)); // initial memory allocation
	Bndry.na_conn_ind = malloc( (inMesh.NumNodes+1) * sizeof(int)); // memory allocation
	double *Lbound;   // array to store boundary intergral coefficents
	int *Lbound_nums; // global node nums corresponding to Lbound entries
	int num_Lbound;   // length of Lbound & Lbound_nums
	double *lambda = calloc(1+lsprob.num, sizeof(double)); // lambda values

	// stuff for Ku = f solve using MA57 solver (& eigenvalue solver)
	sp_mat Kg, Mg; // global stiffness and mass matrices
	Kg.irn=0; Kg.jcn=0; Kg.A=0; Mg.irn=0; Mg.jcn=0; Mg.A=0;
	int num_eig = (lsprob.obj == 3) ? 1 : 0; // number of eigenvalues to compute
	int comp_eig = (lsprob.obj == 3) ? 1 : 0; // flag to decide if eigenvalues are to be computed
	double *eig_vals, *eig_vecs; // arrays for eigenvalues and vectors
	double *save_freq; // used to store freq values for output
	int order = (inMesh.NumNodes * NUM_DOF); // full matrix order
	int numEnt = inMesh.NumElem * ((NUM_DOF*4)+1)*(NUM_DOF*2); // max number of index entries (symmetry!)

	// bar additonal design variables
	if(inMesh.bars){numEnt += inMesh.NumBars * 4;} // additional entries for bar elements
	double *bar_sens=0; // bar senstivities
	double *bar_step=0; // bar update step
	double *bar_max=0, *bar_min=0; // move limits on bar areas
	if(inMesh.bars){j=inMesh.NumBars; bar_max=malloc(j*sizeof(double));
		bar_min=malloc(j*sizeof(double)); bar_step=malloc(j*sizeof(double));
		bar_sens = malloc(j * (1+lsprob.num) * sizeof(double));}

	// designable bc addtional variables
	if(inMesh.des_bc){numEnt += inMesh.NumBC * 4 * NUM_DOF;} // additional entries for designable bcs
	double *bc_sens=0; // bc senstivities
	double *bc_step=0; // bc update step
	double *bc_max=0, *bc_min=0; // move limits on bc variables
	if(inMesh.des_bc){j=inMesh.NumBC; bc_max=malloc(j*sizeof(double));
		bc_min=malloc(j*sizeof(double)); bc_step=malloc(j*sizeof(double));
		bc_sens = malloc(j * (1+lsprob.num) * sizeof(double));}
	int *bc_con=0, *bc_con_ind=0;
	double *bc_dist=0;

	// designable material varibles
	double *mat_sens=0; // material sensitivites
	double *mat_step=0; // bc update step
	double *mat_max=0, *mat_min=0; // move limits on bc variables
	bool mat_opt = inMesh.dm_sim; // flag to switch between seq & sim opt
		if(inMesh.des_mat){j=inMesh.NumDesMat; mat_max=malloc(j*sizeof(double));
			mat_min=malloc(j*sizeof(double)); mat_step=malloc(j*sizeof(double));
			mat_sens = malloc(j * (1+lsprob.num) * sizeof(double));}

	double *disp; // displacement array
	double *load_sw; // load vector, including self_weight load vector
	double *adjont; // adjoint load (then disp) vector - for non-self adjoint problems
	int dispLen = order * numCase; // length of disp array (inc fixed dof)
	int loadLen = freeDof * numCase; // length of load array (excld fixed dof)

	int num_adj = 0; // number of adjoint cases to solve (objective, then in order of constraint list)
	for(i=0;i<lsprob.num;i++)
	{
		if(lsprob.con[i].type==5){ num_adj++; } // displacement constraint
	}

	int num_sens = 0; // number of senstivity calculations
	// compliance, 1st eigenvalue or mechanism objective
	if(lsprob.obj==1 || lsprob.obj==3 || lsprob.obj==5){ num_sens++; }
	for(i=0;i<lsprob.num;i++)
	{
		j = lsprob.con[i].type;
		if(j==3 || j==5 || j==7 || j==8){ num_sens++; } // compliance or displacement constraint
		else if(j==6){num_sens++; num_eig++;} // eigenvalue ratio constraint
	}
	if(lsprob.obj==5){num_sens++;}
	if(lsprob.obj==3){numCase += num_eig;} // 1 primary state for each eigenvector

	double **prim_ptr; // vector, length = numCase (+ num_eig)
	double **dual_ptr; // matrix row num = load case num, col num = dual state num
						// or row = eigenvalue num, col = dual state num
	double *wgt; // weights for load cases

	if(lsprob.obj==5){ prim_ptr = malloc(3 * sizeof(double*));
						dual_ptr = malloc(3 * sizeof(double*));
						wgt = malloc(sizeof(double));
						wgt[0] = 1.0; } // special cosideration for compliant mech
	else{ prim_ptr = malloc(numCase * sizeof(double*));
			dual_ptr = malloc(numCase * num_sens * sizeof(double*));
			wgt = malloc(numCase * sizeof(double));
			for(i=0;i<numCase;i++){wgt[i]=1.0;} } // set to 1.0 (for now)

	// compute load magnitudes for compliant mech design
	if(lsprob.obj==5)
	{
		for(i=0;i<lsprob.num;i++)
		{
			if(lsprob.con[i].type==7)
			{
				// sum total applied loads
				p1=0.0; p2=0.0; k=freeDof;
				for(j=0;j<freeDof;j++)
				{
					p1 += fabs(load[j]); // 1st load case
					p2 += fabs(load[k++]); // 2nd load case
				}
			}
		}
	}

	// Arrays to store objective & constraint data for each iteration
	double *obj = malloc(control.maxItt*sizeof(double));
	double *cnstr = malloc(lsprob.num*control.maxItt*sizeof(double));
	double *pred = malloc((1+lsprob.num)*control.maxItt*sizeof(double)); // predicted obj & constraint change
	double *pred_temp = malloc((1+lsprob.num)*sizeof(double)); // temp to feed SLPsubSol
	double *delcon = malloc(lsprob.num*sizeof(double)); // array to store distance from constraint (constraint violation)
	double *relax = malloc(lsprob.num*sizeof(double)); // relaxation factor
	if(comp_eig==1){save_freq = malloc(num_eig*2*control.maxItt*sizeof(double)); } // frequencies
	for(i=0;i<lsprob.num;i++){relax[i] = 1.0;} // initialize to 1

	// jeehanglee@gmail.com -- temp code
	CMesh cmesh;
	CMaterial cmat;
	CSolver csolver;
	CSensitivity csensitivity;
	CBoundary cboundary;
	CLevelSet clevelset;

	do {
		printf("\n\n------------\nStart iteration %i ...",itt);
		/*--------------------------------------------------------------/
		/																/
		/	   Analyse Structure, Solve Primary & Adjoint Equations		/
		/																/
		/--------------------------------------------------------------*/

		// Discretize the boundary & compute area ratios for AFG method
		cmesh.Find_struct(&inMesh, &levelset, &Bndry, alpha, control.aMin);
		Ntot = inMesh.NumNodes + Bndry.NumAux; // total number of nodes

		if(control.pinfo > 1 || itt==0)	{
			// output lsf (& alpha) values
			coutput.OutPLotShapeVTK2(&inMesh, levelset.lsf, alpha, control.pinfo, itt, filename);
		}

		if(control.pinfo > 1) {
			printf("\nTotal nodes (inc. auxillary) = %i",Ntot);
		}

		// Assemble global matrices (in triplet form)

		// Initialize Index Arrays using dynamic memory allocation
		Kg.ne = numEnt; set_sp_mat(&Kg);

		if(comp_eig==1)
		{
			Mg.ne = numEnt+lump_mass.ne; set_sp_mat(&Mg); // mass matrix
			eig_vals = malloc(num_eig * 2 * sizeof(double)); // eigenvalues
			eig_vecs = malloc(num_eig * 2 * order * sizeof(double)); // eigenvectors
		}

		// calculate total structure volume (area) by summing element area array
		ftemp = 0.0;
		for(i=0;i<inMesh.NumElem;i++) {
			ftemp += alpha[i];
		}
		ftemp *= AreaElem;
		printf("\nTotal volume = %12.4e",ftemp);
		if(lsprob.obj==2) // if volume is objective or constraint
		{
			obj[itt] = ftemp; // store value
		}
		else
		{
			for(i=0;i<lsprob.num;i++){
				if(lsprob.con[i].type==1){
					j = lsprob.num * itt + i; // point to place in cnstr
					cnstr[j] = ftemp; // store volume as constraint value
				}
			}
		}

		// compute total mass (if required)
		for(i=0;i<lsprob.num;i++)
		{
			if(lsprob.con[i].type==2)
			{
				ftemp = 0.0;
				if(inMesh.des_mat){temp=0; oMin = inMat[inMesh.mat2].rho / inMat[inMesh.mat1].rho;}
				for(j=0;j<inMesh.NumElem;j++) {
					oMax = alpha[j] * inMat[inMesh.mat_type[j]].rho * inMesh.t; // area ratio x density
					if(inMesh.des_mat && inMesh.mat_elems[temp]==j) {
						//oMax *= inMesh.mat_vars[temp] + (1.0-inMesh.mat_vars[temp])*oMin; // modify for current variable
	                    oMax *= inMesh.mat_vars[temp]*oMin + (1.0-inMesh.mat_vars[temp]); // modify for current variable
						temp++; // next variable
					}
					ftemp += oMax;
				}
				ftemp *= inMesh.h * inMesh.h; // multiply by side length sqrd to get total mass
				printf("\nTotal mass (plane) = %12.4e",ftemp);

				// add bar mass (if there are some)
				if(inMesh.bars)
				{
					oMax = 0.0;
					oMin = inMesh.bar_mat->rho * inMesh.h; // density x length
					for(j=0;j<inMesh.NumBars;j++) {
						oMax += inMesh.bar_areas[j] * oMin; // bar mass
					}
					printf("\nTotal mass (bar) = %12.4e",oMax);
					ftemp += oMax;
				}

				j = lsprob.num * itt + i; // point to place in cnstr
				cnstr[j] = ftemp; // store mass as constraint value
				printf("\nTotal mass = %12.12e",ftemp);
			}
		}

	    // compute % of each material if using designable materials
	    if(inMesh.des_mat)
	    {
	        temp = 0;
	        oMax = 0.0; // % material 2
	        oMin = 0.0; // sum alpha
	        for(j=0;j<inMesh.NumElem;j++)
	        {
	            if(inMesh.des_mat && inMesh.mat_elems[temp]==j)
	            {
	                oMin += alpha[j];
	                oMax += inMesh.mat_vars[temp]*alpha[j]; // multiply by area ratio
	                temp++; // next variable
	            }
	        }
	        oMax /= oMin; // % material 2
	        printf("\nPcnt material 1 = %f",1.0-oMax);
	        printf("\nPcnt material 2 = %f",oMax);
	    }

		// asemble gloabl stiffness (& maybe mass) matrix
		AFG_Matrix(comp_eig, KE, ME, &Kg, &Mg, &lump_mass, alpha, &inMesh, inMat, control.aMin, control.mMin);

		// Solve the FE equation using the MA57 solver

		// solve Ku = f (for all cases)
		// remove fixed dof from global matrix
	    rem_sp_mat(&Kg, fixDof, 1);

		if(lsprob.obj != 3)
		{
	        disp = malloc(dispLen * sizeof(double));
	        // compute total load vector (in case of self-weight)
	        if(sw)
	        {
	            load_sw = malloc(dispLen * sizeof(double));
	            cmat.self_weight(&inMesh, inMat, control.aMin, control.mMin, alpha, freeDof, fixDof, numCase, load, load_sw, acc);
	            cblas_dcopy(loadLen, load_sw, 1, disp, 1); // copy in load to send to FE_solve
	        }
	        else
	        {
	            cblas_dcopy(loadLen, load, 1, disp, 1); // copy in load to send to FE_solve
	        }
			csolver.FE_Solve(&Kg, fixDof, disp, freeDof, order, numCase);
		}

		// compute eigenvalues (if required)
		if(comp_eig==1)
		{
			rem_sp_mat(&Mg, fixDof, 1); // remove fixed dof from global matrix
			csolver.eig_solve(num_eig, &Kg, &Mg, freeDof, order, fixDof, eig_vals, eig_vecs, control.pinfo);
			if(lsprob.obj==3){ obj[itt] = sqrt(eig_vals[0]); } // 1st freq
			for(i=0;i<lsprob.num;i++){
				if(lsprob.con[i].type==6){
					j = lsprob.num * itt + i; // point to place in cnstr
					cnstr[j] = sqrt(eig_vals[1] / eig_vals[0]); // store eigenvalue ratio
					printf("\nEigen freq ratio = %f",cnstr[j]);
				}
			}

			// save frequencies
			j = num_eig*2; k = j*itt;
			for(i=0;i<j;i++){ save_freq[i+k] = sqrt(eig_vals[i]); }
		}

		// calculate total compliance (if required)
		temp=0;
		if(lsprob.obj==1){temp=1;}
		else {
			for(i=0;i<lsprob.num;i++){ if(lsprob.con[i].type==3){temp=1; break;} }
		}
		if(temp==1) // if compliance is the objective or a constraint
		{
			obj_val = 0.0;
			for(j=0;j<numCase;j++)
			{
				ftemp = 0.0;
				k = freeDof*j;  // point to start of load array for case j
				j2 = order*j;   // point to start of disp array for case j
				for(i=0;i<order;i++) {
					if(fixDof[i] > -1) {
	                    if(sw){ ftemp += disp[j2+i]*load_sw[k+fixDof[i]]; }
						else{ ftemp += disp[j2+i]*load[k+fixDof[i]]; }
					}
				}
				printf("\ncompliance for case %i = %12.4e",j+1,ftemp);
				obj_val += ftemp * wgt[j]; // multiply by load case weight
			}

			printf("\nTotal compliance = %12.12e",obj_val);

			if(lsprob.obj==1){obj[itt]=obj_val;}
			for(i=0;i<lsprob.num;i++){
				if(lsprob.con[i].type==3){
					j = lsprob.num * itt + i; // point to place in cnstr
					cnstr[j] = obj_val; // store compliance as constraint value
				}
			}
		}

	    if(sw){free(load_sw);} // no longer needed

		// calculate constraint values and adjoints for displacement constraints
		if(num_adj>0){adjont = calloc(order*num_adj, sizeof(double));}
		temp=0; // count adjoint load cases
		for(i=0;i<lsprob.num;i++)
		{
			if(lsprob.con[i].type==5)
			{
				j = (lsprob.num * itt) + i; // point to place in cnstr
				k = (order * (int)lsprob.con[i].data[2]) + (int)lsprob.con[i].data[1]; // place in disp vector
				// compute constraint value
				ftemp = disp[k] - lsprob.con[i].data[0];
				//printf("\nDisplacement constraint %i violated by %12.4e\n",i+1,ftemp);

				i2 = (int)(lsprob.con[i].data[1]);
				k = (temp * freeDof) + fixDof[i2]; // place in adjoint load vector

				cnstr[j] = fabs(ftemp); // quicker than ^2 then sqrt!
				adjont[k] = ftemp/cnstr[j]; // adjoint load
				printf("\nAdjoint load %i = %12.4e\n",i+1,adjont[k]);
				temp++; // next adjoint
			}

			// compliant mechanism disp constraint
			else if(lsprob.con[i].type==7)
			{
				// obtain the displacement values
				u11=0.0; u12=0.0; u21=0.0; u22=0.0;

				// compute u11, etc ...
				for(j=0;j<order;j++)
				{
					if(fixDof[j] != -1)	// if dof not fixed copy accross the value
					{
						u11 += load[fixDof[j]] * disp[j];
						u12 += load[fixDof[j]] * disp[order+j];
						u21 += load[fixDof[j]+freeDof] * disp[j];
						u22 += load[fixDof[j]+freeDof] * disp[order+j];
					}
				}

				u11 /= p1;
				u22 /= p2;
				u21 /= p2;
				u12 /= p1;

				printf("\nu11=%12.4e, u21=%12.4e, u12=%12.4e, u22=%12.4e",u11,u21,u12,u22);

				// compute reaction force
				ftemp = u22;
				if(lsprob.con[i].data[3] > 0.0){ ftemp += p2/lsprob.con[i].data[3]; } // add spring stiffness
				ftemp = (p2*u21)/ftemp;

				// objective
				obj[itt] = ftemp / p1; // output / input force
				printf("\nMechanical advantage = %12.4e",obj[itt]);

				// input displacement
				j = (lsprob.num * itt) + i; // point to place in cnstr
				if(lsprob.con[i].data[3] > 0.0)
				{
					ftemp = u21 - (ftemp/lsprob.con[i].data[3]);
					ftemp /= u22;
					cnstr[j] = u11 - ftemp*u12;
				}
				else { cnstr[j] = u11 - (u21*u12)/u22; }
				printf("\nInput disp = %12.4e",cnstr[j]); // current input displacement
			}
		}

		if(num_adj>0)
		{
			// solve Kp = g (for all adjoint cases)
			csolver.FE_Solve(&Kg, fixDof, adjont, freeDof, order, num_adj);
			// this could be more efficient (by saving the factorization)
		}

		// clear memory for K (&M) matrices
		free_sp_mat(&Kg); if(comp_eig==1){free_sp_mat(&Mg);}

		if(control.pinfo == 3) {
			// write displacement & mode shape info file
			temp = (lsprob.obj==3) ? 0 : numCase;
			coutput.OutDispVTK(&inMesh, temp, disp, 2*num_eig, eig_vecs, itt, filename);
		}

		/*--------------------------------------------------------------/
		/																/
		/		Calculate shape sensitivities and velocity function		/
		/																/
		/--------------------------------------------------------------*/

		// Calculate boundary node sensitivities
		Nsens = calloc(Ntot*num_sens,sizeof(double));

		// point to correct senstivities etc ...

		// setup array of pointers to correct places in disp (for each case)
		if(lsprob.obj != 3 && lsprob.obj != 5)
		{
			for(i=0;i<numCase;i++)
			{
				j = order * i;
				prim_ptr[i] = &disp[j]; // primary variables always displacement
			}
		}

		// objective
		temp = 0; // count senstivities (should = num_sens by end)
		temp2 = 0; // count adjoint states (should = num_adj by end)
		if(lsprob.obj==1){
			sens_ptr[0] = &Nsens[temp++]; // compliance sens
			for(i=0;i<numCase;i++){dual_ptr[i*num_sens] = prim_ptr[i];} // compliance is self-adjoint
		}
		else if(lsprob.obj==2){vol_sens = malloc(Ntot*sizeof(double));
								for(k=0;k<Ntot;k++){ vol_sens[k] = -1.0; }
								sens_ptr[0] = vol_sens;}   // volume sens
		else if(lsprob.obj==3){
			sens_ptr[0] = &Nsens[temp++]; // 1st eigenvalue sensitivity
			prim_ptr[0] = eig_vecs; // eigenvector
			dual_ptr[0] = prim_ptr[0];
		}
		else if(lsprob.obj==5) {
			sens_ptr[0] = Nsens; // first array of Nsens for objective
		}

		// constraints
		for(i=0;i<lsprob.num;i++)
		{
			j = lsprob.con[i].type;
			if(j==1)
			{
				vol_sens = malloc(Ntot*sizeof(double));
				for(k=0;k<Ntot;k++){ vol_sens[k] = -1.0; } // sensitvity for volume is -1 (geometric constraint)
				sens_ptr[i+1] = vol_sens;
			}
			else if(j==2)
			{
				mass_sens = calloc(Ntot,sizeof(double));

				if(inMesh.des_mat)
				{
					i2=0; j2=0; // i2 = material element, j2 = boundary segment
					for(k=0;k<inMesh.NumElem;k++)
					{
						if(k==Bndry.Bound[j2].e)
						{
							ftemp = inMesh.mat_vars[i2];
							//ftemp = ftemp*inMat[inMesh.mat1].rho + (1.0-ftemp)*inMat[inMesh.mat2].rho;
	                        ftemp = (1.0-ftemp)*inMat[inMesh.mat1].rho + ftemp*inMat[inMesh.mat2].rho;
							ftemp *= -0.5 * inMesh.t;
							mass_sens[Bndry.Bound[j2].n1] += ftemp;
							mass_sens[Bndry.Bound[j2].n2] += ftemp;

							j2++; // update count

							if(j2 < Bndry.NumBound && k==Bndry.Bound[j2].e) // 2 segments / elem
							{
								mass_sens[Bndry.Bound[j2].n1] += ftemp;
								mass_sens[Bndry.Bound[j2].n2] += ftemp;

								j2++; // update count
							}
						}

						if(k==inMesh.mat_elems[temp]){i2++;} // update count of material element
					}
				}

				else
				{
					ftemp = -inMat[0].rho * inMesh.t; // mass senstivity (assumed) density x thickness
					for(k=0;k<Ntot;k++){ mass_sens[k] = ftemp; }
				}

				sens_ptr[i+1] = mass_sens;
			}
			else if(j==3)
			{
				sens_ptr[i+1] = &Nsens[temp * Ntot]; // compliance sens
				for(k=0;k<numCase;k++){dual_ptr[(k*num_sens)+temp] = prim_ptr[k];} // compliance is self-adjoint
				temp++;
			}
			else if(j==5)
			{
				sens_ptr[i+1] = &Nsens[temp * Ntot]; // displacement sens
				i2 = (int)lsprob.con[i].data[2]; // load case number
				for(k=0;k<numCase;k++)
				{
					if(k==i2){dual_ptr[(k*num_sens)+temp] = &adjont[temp2 * order]; temp2++;}
					else{dual_ptr[(k*num_sens)+temp] = NULL;} // no adjoint state (not realted to load case)
				}
				temp++;
			}
			else if(j==6)
			{
				sens_ptr[i+1] = &Nsens[temp * Ntot]; // 2nd eigenvalue sensitivity
				prim_ptr[1] = &eig_vecs[order]; // 2nd eigenvector
				dual_ptr[1] = prim_ptr[1];
				temp++;
			}
			else if(j==7)
			{
				sens_ptr[i+1] = &Nsens[Ntot]; // 2nd array of Nsens
			}
			else if(j==8)
			{
				sens_ptr[i+1] = &Nsens[2*Ntot]; // 3rd array of Nsens
			}
			else if(j==10)
			{
				zero_sens = calloc(Ntot,sizeof(double));
				sens_ptr[i+1] = zero_sens;
			}
		}

		if(lsprob.obj == 3)
		{
			// compute boundary shape sensitivies (for eigenvalues)
			// MAT VAR EDIT
			csensitivity.AFG_Sens(&inMesh, &Bndry, alpha, inMat, Nsens, prim_ptr, dual_ptr, num_sens, num_eig, eig_vals, gCoord, control.aMin, 1,fact,sw,acc);
			for(i=0;i<lsprob.num;i++)
			{
				if(lsprob.con[i].type==6)
				{
					// process individual senstivities for constraint
					ftemp = cnstr[lsprob.num * itt + i] / (2.0*eig_vals[0]);
					oMax = 1.0 / (2.0*sqrt(eig_vals[0]*eig_vals[1]));
					//oMin = 1.0 / (2.0*sqrt(eig_vals[0]));  // factors
					k = Ntot;
					for(j=0;j<Ntot;j++)
					{
						Nsens[k] = oMax*Nsens[k] - ftemp*Nsens[j]; // swap sign for reduction
						//Nsens[j] *= oMin; // do not scale 1st eigenvalue senstivity
						k++;
					}
				}
			}
		}
		else if(lsprob.obj == 5)
		{
			// combine senstivities together
			for(i=0;i<lsprob.num;i++)
			{
				if(lsprob.con[i].type==7){temp = i; break;} // displacement constraint
			}

			// compute reaction force
			ftemp = u22;
			if(lsprob.con[temp].data[3] > 0.0){ ftemp += p2/lsprob.con[temp].data[3];} // add spring stiffness

			fact[0] = (1.0/p1)*ftemp;
			fact[1] = -fact[0] * (u21 / ftemp);

			ftemp = (lsprob.con[temp].data[3] > 0.0) ? p2/(ftemp*lsprob.con[temp].data[3]) : 0.0;
			fact[2] = (u12/u22)*(ftemp-1.0)/p2;
			fact[2] += (u21/u22)*(ftemp-1.0)/p1;

			if(lsprob.con[temp].data[3] > 0.0)
			{
				oMax = 1.0 - ftemp - (u22/( u22 + (p2/lsprob.con[temp].data[3]) ))*ftemp;
				ftemp = oMax/p2;
			}
			else { ftemp = 1.0/p2; }
			fact[3] = (u21*u12)/(u22*u22);
			fact[4] = 1.0/p1;
			fact[5] = 1.0/p2;

			prim_ptr[0] = disp;

			dual_ptr[0] = prim_ptr[0];
			dual_ptr[1] = &disp[order];

			// compute sensitivities
			csensitivity.AFG_Sens(&inMesh, &Bndry, alpha, inMat, Nsens, prim_ptr, dual_ptr, 3, 1, wgt, gCoord, control.aMin, 2,fact,sw,acc);
		}
		else
		{
			// compute boundary shape sensitivies (for compliance or displacement functions)
			csensitivity.AFG_Sens(&inMesh, &Bndry, alpha, inMat, Nsens, prim_ptr, dual_ptr, num_sens, numCase, wgt, gCoord, control.aMin, 0,fact,sw,acc);
			free(disp);
		}

		// bar senstivities
		if(inMesh.bars && lsprob.obj == 1)
		{
			// call routine for computing bar compliance senstivites
			// inlcuding multi-load cases
			for(j=0;j<inMesh.NumBars;j++) {
				bar_sens[j] = 0.0; // zero initial
			}

			csensitivity.barSens(&inMesh, bar_sens, prim_ptr, numCase, wgt);

			// mass senstivites
			for(i=0;i<lsprob.num;i++)
			{
				if(lsprob.con[i].type==2)
				{
					ftemp = inMesh.bar_mat->rho * inMesh.h; // bar mass senstivity (wrt area)
					k = inMesh.NumBars * (i+1);
					for(j=0;j<inMesh.NumBars;j++) {
						bar_sens[k++] = ftemp; // mass constraint
					}
				}
			}

			// also set move limits on bar areas here
			ftemp = 0.05 * (inMesh.bar_max - inMesh.bar_min);
			for(j=0;j<inMesh.NumBars;j++)
			{
				oMax = inMesh.bar_areas[j] + ftemp;
				oMin = inMesh.bar_areas[j] - ftemp;
				bar_max[j] = (oMax < inMesh.bar_max) ? ftemp : inMesh.bar_max - inMesh.bar_areas[j];
				bar_min[j] = (oMin > inMesh.bar_min) ? -ftemp : inMesh.bar_min - inMesh.bar_areas[j];
			}

			// output bar info
			coutput.OutBars(&inMesh, (lsprob.num+1), bar_sens, control.pinfo, itt, filename);
		}

		// additional bc sensitivities
		if(inMesh.des_bc && lsprob.obj == 1)
		{
			for(j=0;j<inMesh.NumBC;j++) {
				bc_sens[j] = 0.0; // zero initial
			}

			// CALL SENSITIVITY FUNCTION
			csensitivity.bcSens(&inMesh, bc_sens, prim_ptr, numCase, wgt);

			// Sensitivity filtering
			if(itt==0)
			{
				// analyse connectivity of des bc elements
				int xMin, xMax, yMin, yMax, n, m;
				temp = inMesh.NumBC;
				bc_con = malloc(8*temp*sizeof(int));
				bc_dist = malloc(8*temp*sizeof(double));
				bc_con_ind = malloc((temp+1)*sizeof(int));
				int count = 0;
				for(i=0;i<temp;i++)
				{
					// compute element indices
					ftemp = (double)inMesh.BC_nums[i] / (double)inMesh.elemX;
					m = floor(ftemp);
					n = inMesh.BC_nums[i] - m*inMesh.elemX;
					bc_con_ind[i] = count;

					// search for neighbouring elements with designable bc
					xMin = n-1; xMax=n+2;
					yMin = m-1; yMax=m+2; // search limits
					xMin = (xMin < 0) ? 0 : xMin;
					yMin = (yMin < 0) ? 0 : yMin;
					xMax = (xMax > inMesh.elemX) ? inMesh.elemX : xMax;
					yMax = (yMax > inMesh.elemY) ? inMesh.elemY : yMax; // adjust search limits

					for(j2=yMin;j2<yMax;j2++)
					{
						for(i2=xMin;i2<xMax;i2++)
						{
							temp2 = j2*inMesh.elemX + i2; // elem number
							for(j=0;j<temp;j++)
							{
								if(temp2 == inMesh.BC_nums[j] && j != i)
								{
									// add to conn array
									if(m==j2 || n==i2) { ftemp = 0.5; }
									else { ftemp = 0.085786437626905; } // 1.5 - distance
									bc_con[count] = j;
									bc_dist[count++] = ftemp;
								}
							}
						}
					}
				}
				bc_con_ind[temp] = count;
				bc_con = realloc(bc_con, count * sizeof(int));
				bc_dist = realloc(bc_dist, count * sizeof(double));
			}

			// Filter!
			{
				temp = inMesh.NumBC;
				double *bc_temp = malloc(temp*sizeof(double));
				cblas_dcopy(temp, bc_sens, 1, bc_temp, 1); // copy original

				for(i=0;i<temp;i++) // for all des bcs
				{
					temp2 = bc_con_ind[i+1];
					oMax = 1.5*inMesh.K_bc[i]*bc_temp[i];
					oMin = 1.5;
					for(j=bc_con_ind[i];j<temp2;j++) // for all connected elems
					{
						k = bc_con[j]; // connectd bc number
						oMax += bc_dist[j] * inMesh.K_bc[k] * bc_temp[k];
						oMin += ftemp;
					}

					bc_sens[i] = oMax / (oMin * inMesh.K_bc[i]); // filtered sensitivity
				}

				free(bc_temp);
			}

			// scale maximum sensitiviy to 1
			oMax = 0.0;
			for(j=0;j<inMesh.NumBC;j++) {
				ftemp = fabs(bc_sens[j]);
				oMax = (ftemp > oMax) ? ftemp : oMax;
			}
			for(j=0;j<inMesh.NumBC;j++) {
				bc_sens[j] /= oMax;
			}

			// other senstivities
			for(i=0;i<lsprob.num;i++)
			{
				// constraint cost sensitivities are 1.0
				// also compute the constraint value
				if(lsprob.con[i].type==10)
				{
					ftemp = 0.0;
					k = inMesh.NumBC;
					for(j=0;j<inMesh.NumBC;j++) {
						bc_sens[k++] = 1.0;
						ftemp += inMesh.K_bc[j];
					}

					j = (lsprob.num * itt) + i; // point to place in cnstr
					cnstr[j] = ftemp;
				}
			}

			// also set move limits on bc variables here
			for(j=0;j<inMesh.NumBC;j++)
			{
				oMax = inMesh.K_bc[j] + 0.05;
				oMin = inMesh.K_bc[j] - 0.05;
				bc_max[j] = (oMax < 1.0) ? 0.05 : 1.0 - inMesh.K_bc[j];
				bc_min[j] = (oMin > 1.0e-4) ? -0.05 : 1.0e-4 - inMesh.K_bc[j];
			}

			// output designable BC info
			coutput.OutDesBC(&inMesh, bc_sens, control.pinfo, itt, filename);
		}

		// material sensitvities
		if(inMesh.des_mat)
		{
	        double *mat_sens_temp = calloc(numCase * inMesh.NumDesMat, sizeof(double));

	        // call routine for computing material complinace senstivites
	        if(lsprob.obj == 1)
	        {
	            // compliance objective sensitivity
	            if(inMesh.mat_lin){ csensitivity.matSens_comp(&inMesh, inMat, KE[0], mat_sens_temp, numCase, wgt, prim_ptr, dual_ptr, alpha, control.aMin, sw, acc); }
	            //else{ HS_Sens_eig(&inMesh, inMat, KE[0], ME[0], mat_sens_temp, num_eig, eig_vals, prim_ptr, alpha, control.aMin); }
	            cblas_dcopy(inMesh.NumDesMat, mat_sens_temp, 1, mat_sens, 1); // copy objective sensitivities
	        }

			// call routine for computing material eig senstivites
	        if(lsprob.obj == 3)
	        {
	            // eigenvalue objective sensitivity
	            if(inMesh.mat_lin){ csensitivity.matSens_eig(&inMesh, inMat, KE[0], ME[0], mat_sens_temp, num_eig, eig_vals, prim_ptr, alpha, control.aMin); }
	            else{ csensitivity.HS_Sens_eig(&inMesh, inMat, KE[0], ME[0], mat_sens_temp, num_eig, eig_vals, prim_ptr, alpha, control.aMin); }
	            cblas_dcopy(inMesh.NumDesMat, mat_sens_temp, 1, mat_sens, 1); // copy objective sensitivities
	        }

	        // constraint senstivites
	        for(i=0;i<lsprob.num;i++)
	        {
	            // volume constraint sensitivities
	            if(lsprob.con[i].type==1)
	            {
	                k = inMesh.NumDesMat * (i+1);
	                for(j=0;j<inMesh.NumDesMat;j++) {
	                    mat_sens[k++] = 0.0; // material choice does not effect volume
	                }
	            }

	            else if(lsprob.con[i].type==2)
	            {
	                // material mass senstivity
	                ftemp = (inMat[inMesh.mat2].rho - inMat[inMesh.mat1].rho) * inMesh.h * inMesh.h * inMesh.t;
	                k = inMesh.NumDesMat * (i+1);
	                for(j=0;j<inMesh.NumDesMat;j++) {
	                    mat_sens[k++] = alpha[inMesh.mat_elems[j]] * ftemp; // mass sensitivity
	                }
	            }

	            // freq ratio sensitivites
	            else if(lsprob.con[i].type==6)
	            {
	                // process individual senstivities for constraint
	                ftemp = cnstr[lsprob.num * itt + i] / (2.0*eig_vals[0]);
	                oMax = 1.0 / (2.0*sqrt(eig_vals[0]*eig_vals[1]));
	                //oMin = 1.0 / (2.0*sqrt(eig_vals[0]));  // factors
	                k = inMesh.NumDesMat * (i+1); // place in mat_sens for freq ratio sens
	                temp = inMesh.NumDesMat; // start of 2nd eigenvalue sens
	                for(j=0;j<inMesh.NumDesMat;j++)
	                {
	                    mat_sens[k++] = oMax*mat_sens_temp[temp++] - ftemp*mat_sens_temp[j]; // swap sign for reduction
	                    //Nsens[j] *= oMin; // do not scale 1st eigenvalue senstivity
	                }
	            }
	        }

	        // find max for FD check
	        /*oMax=0.0;
	        j=-1;
	        for(i=0;i<inMesh.NumDesMat;i++)
	        {
	            ftemp = fabs(mat_sens[i]);
	            if(ftemp > oMax){ oMax=ftemp; j=i;}
	        }

	        printf("\nMax obj sens %i = %12.12e",j,mat_sens[j]);
	        for(i=0;i<lsprob.num;i++)
	        {
	            k = inMesh.NumDesMat * (i+1);
	            printf("\nConst sens %i = %12.12e",i,mat_sens[k+j]);
	        }*/

	        free(mat_sens_temp);

			// also set move limits on material variables
			ftemp = 0.02;
			for(j=0;j<inMesh.NumDesMat;j++)
			{
				oMax = inMesh.mat_vars[j] + ftemp;
				oMin = inMesh.mat_vars[j] - ftemp;
				mat_max[j] = (oMax < 1.0) ? ftemp : 1.0 - inMesh.mat_vars[j];
				mat_min[j] = (oMin > 0.0) ? -ftemp : 0.0 - inMesh.mat_vars[j];
			}

			// output material info
			coutput.OutDesMat(&inMesh, alpha, control.aMin, (lsprob.num+1), mat_sens, control.pinfo, itt, filename);
		}

		// free required memory
		if(comp_eig==1) { free(eig_vals); free(eig_vecs); }
		if(num_adj>0){ free(adjont); }

		if(control.pinfo == 3)	{
			// Write node sensitvity information file (as a boundary mesh)
			coutput.OutBoundVTK(&inMesh, &Bndry, (lsprob.num+1), sens_ptr, itt, filename);
		}

		// Check convergence criteria
		// calculate difference between current and target constraint
		temp = 1; // check for violated constraints
		j = lsprob.num;
		for(i=0;i<j;i++)
		{
			k = j*itt + i; // pointer for cnstr
			ftemp = cnstr[k]; // current constraint value
			if(lsprob.con[i].type==1)
			{
				delcon[i] = lsprob.con[i].data[0] - ftemp;
				printf("\nVolume constraint violated by %12.4e\n",delcon[i]);
				if(delcon[i]/lsprob.con[i].data[0] < -control.gm){temp=0;}
			}
			else if(lsprob.con[i].type==2)
			{
				delcon[i] = lsprob.con[i].data[0] - ftemp;
				printf("\nMass constraint violated by %12.4e\n",delcon[i]);
				if(delcon[i]/lsprob.con[i].data[0] < -control.gm){temp=0;}
			}
			else if(lsprob.con[i].type==3)
			{
				delcon[i] = lsprob.con[i].data[0] - ftemp;
				printf("\nCompliance constraint violated by %12.4e\n",delcon[i]);
				if(delcon[i] < 0.0){temp=0;}
			}
			else if(lsprob.con[i].type==5)
			{
				delcon[i] = -ftemp;
				printf("\nDisplacement constraint %i violated by %12.4e\n",i+1,delcon[i]);
				if(fabs(delcon[i]/lsprob.con[i].data[0]) > 0.01){temp=0;} // within 1%
			}
			else if(lsprob.con[i].type==6)
			{
				delcon[i] = ftemp - lsprob.con[i].data[0]; // reverse sign (as greater than constraint)
				printf("\nEig freq ratio constraint %i violated by %12.4e\n",i+1,-delcon[i]);
				if(delcon[i] < 0.0){temp=0;}
			}
			else if(lsprob.con[i].type==7)
			{
				delcon[i] = lsprob.con[i].data[0] - ftemp;
				printf("\nInput disp constraint %i violated by %12.4e\n",i+1,delcon[i]);
				if(delcon[i] < 0.0){temp=0;}
			}
			else if(lsprob.con[i].type==8)
			{
				delcon[i] = lsprob.con[i].data[0] - ftemp;
				printf("\nOutput comp constraint %i violated by %12.4e\n",i+1,delcon[i]);
				if(delcon[i] < 0.0){temp=0;}
			}
			else if(lsprob.con[i].type==10)
			{
				delcon[i] = lsprob.con[i].data[0] - ftemp;
				printf("\nConstraint cost constraint %i violated by %12.4e\n",i+1,delcon[i]);
				delcon[i] = 1.0e+20; // fake large number
			}
		}

		if(itt > 10+itt0 && temp==1)
		{
			oMax = obj[itt]; // initialize variable for max objective value
			oMin = obj[itt];	// initialize variable for min objective value
			for(i=itt-1;i>itt-10;i--) // work out maximum objective change over last 10 itteration
			{
				ftemp = obj[i];
				oMax = (ftemp > oMax) ? ftemp : oMax; // update maximum objective value
				oMin = (ftemp < oMin) ? ftemp : oMin; // update minimum objective value
			}
			printf("Max obj = %12.4e, Min obj = %12.4e\n",oMax,oMin);
			ftemp = (oMax - oMin) / (oMax + oMin); // normalized max change
			// if max average over last 10 itteration is < gamma then stop!!
			if(ftemp < control.gm) {
	            if(inMesh.des_mat && !inMesh.dm_sim && !mat_opt)
	            {
	                printf("Level-set Converged! Now optimizing material ...");
	                mat_opt = true; // now start optimizing material
	                itt0 = itt;
	            }
	            else
	            {
	                printf("Converged!!! Criteria = %12.4e",ftemp);
	                break;
	            }
			}
			else {
				printf("Convergence criteria = %12.4e",ftemp);
			}
		}

		// relaxation for constraint change estimates & cfl (on objective)
		if(itt>0)
		{
			// for each constraint
			for(i=0;i<lsprob.num;i++)
			{
				if(active[i] != 0)
				{
					// multiply target constraint change by: actual / predicted from last itt
					j=(1+lsprob.num)*(itt-1) + i+1; // place in predicted (last itt)

					ftemp = cnstr[itt*lsprob.num + i];
					ftemp -= cnstr[(itt-1)*lsprob.num + i];  // actual constraint change
					ftemp = pred[j] / ftemp; // relaxation factor

					// compute constraint relaxation factors
					if(ftemp > 1.0){relax[i] *= 2.0;} // increase
					if(ftemp < 0.5){relax[i] *= 0.5;} // decrease
					if(relax[i] > 1.0){relax[i] = 1.0;} // upper limit
					else if(relax[i] < 0.125){relax[i] = 0.125;} // lower limit

					printf("\nRelaxation factor for %i = %f",i+1,relax[i]);
					delcon[i] *= relax[i];
				}
			}
		}

		// Compute boundary intergral of shape senstivities
		// then use LP sub-problem to obtain the velocity function
		temp = lsprob.num + 1; // number of "functions"
		Vnorm = calloc(Ntot,sizeof(double));		// normal velocity array
		Lbound = malloc(Ntot*temp*sizeof(double));	// boundary intergral data array
		Lbound_nums = malloc(Ntot*sizeof(int)); // global boundary point numbers

		// compute variabe weights for boundary intergral
		cboundary.BsegWgt(&Bndry, &inMesh);

		// compute coefficients for discretized boundary integral
		cboundary.BoundInt(&inMesh, &levelset, &Bndry, temp, sens_ptr, Lbound_nums, &num_Lbound, Lbound);

		if(control.pinfo == 3)	{
			// Write boundary integration information file
			coutput.OutBoundInt((lsprob.num+1), num_Lbound, Lbound_nums, Lbound, itt, filename);
		}

		// realloc & free memory
		Lbound = realloc(Lbound, num_Lbound*temp*sizeof(double));
		Lbound_nums = realloc(Lbound_nums, num_Lbound*sizeof(int));

		// obtian Vnorm using SLP subproblem
		if(inMesh.bars)
		{
			temp = clevelset.SLPsubSol4(&inMesh, &levelset, &Bndry, cfl, delcon, lsprob.num, sens_ptr, Lbound, num_Lbound, Lbound_nums, lambda, active, Vnorm, pred_temp, inMesh.NumBars, bar_sens, bar_min, bar_max, bar_step, control.pinfo);
		}
		else if(inMesh.des_bc)
		{
			// solve for level set velocity
			temp = clevelset.SLPsubSol4(&inMesh, &levelset, &Bndry, cfl, delcon, lsprob.num, sens_ptr, Lbound, num_Lbound, Lbound_nums, lambda, active, Vnorm, pred_temp, 0, 0, 0, 0, 0, control.pinfo);

			// designable boundary condition cost constraint
			// shift for LP sub-solve
			for(i=0;i<lsprob.num;i++)
			{
				if(lsprob.con[i].type==10)
				{
					j = (lsprob.num * itt) + i; // point to place in cnstr
					ftemp = lsprob.con[i].data[0] - cnstr[j];
					ftemp -= cblas_ddot(inMesh.NumBC, &bc_sens[inMesh.NumBC], 1, bc_min, 1);
					break;
				}
			}

			for(i=0;i<inMesh.NumBC;i++)
			{
				bc_max[i] -= bc_min[i]; // adjust so min is zero
			}

			// Can directly solve using LP sub-problem
			csolver.LPsolve(inMesh.NumBC, 1, inMesh.NumBC, bc_step, bc_sens, &bc_sens[inMesh.NumBC], &ftemp, bc_max, control.pinfo);

			for(i=0;i<inMesh.NumBC;i++)
			{
				bc_step[i] += bc_min[i]; // re-adjust
			}
		}
		else if(inMesh.des_mat)
		{
	        if(inMesh.dm_sim || mat_opt)
	        {
	            // reduce number of mat vars using a map
	            bool *inc_mv = malloc(inMesh.NumDesMat*sizeof(bool));
	            int num_reduced=0;
	            for(i=0;i<inMesh.NumDesMat;i++)
	            {
	                j = inMesh.mat_elems[i]; // elem number
	                if(alpha[j] > control.aMin){inc_mv[i] = 1; num_reduced++;}
	                else{inc_mv[i]=0;}
	            }

	            double *mat_sens_red = malloc(num_reduced*(1+lsprob.num)*sizeof(double));
	            double *mat_step_red = malloc(num_reduced*sizeof(double));
	            j=0;
	            for(i=0;i<inMesh.NumDesMat;i++)
	            {
	                // now copy sensitivity & limit data
	                if(inc_mv[i])
	                {
	                    mat_max[j] = mat_max[i];
	                    mat_min[j] = mat_min[i];

	                    for(k=0;k<=lsprob.num;k++)
	                    {
	                        j2 = k*inMesh.NumDesMat + i; // original
	                        i2 = k*num_reduced + j; // reduced
	                        mat_sens_red[i2] = mat_sens[j2]; // copy
	                    }
	                    j++; // next reduced entry
	                }
	            }

	            j = num_reduced * (lsprob.num+1);

	            //temp = SLPsubSol4(&inMesh, cfl, delcon, lsprob.num, sens_ptr, Lbound, num_Lbound, Lbound_nums, lambda, active, Bndry.AuxNodes, Vnorm, pred_temp, inMesh.NumDesMat, mat_sens, mat_min, mat_max, mat_step, control.pinfo);

	            if(inMesh.dm_sim)
	            {
	                // solve for simultaneous update
	                temp = clevelset.SLPsubSol4(&inMesh, &levelset, &Bndry, cfl, delcon, lsprob.num, sens_ptr, Lbound, num_Lbound, Lbound_nums, lambda, active, Vnorm, pred_temp, num_reduced, mat_sens_red, mat_min, mat_max, mat_step_red, control.pinfo);

	                // also ensure smooth transition for materials in elems just out from boundary
	                int n, m;
	                for(i=0;i<inMesh.NumDesMat;i++)
	                {
	                    // if not included
	                    if(!inc_mv[i])
	                    {
	                        // compute element indices
	                        ftemp = (double)inMesh.mat_elems[i] / (double)inMesh.elemX;
	                        m = floor(ftemp);
	                        n = inMesh.mat_elems[i] - m*inMesh.elemX;

	                        // average neighbours that are part of structure
	                        k=0;
	                        ftemp=0.0;
	                        for(j2=m-1;j2<m+1;j2++)
	                        {
	                            for(i2=n-1;i2<n+1;i2++)
	                            {
	                                if(j2<inMesh.elemY && j2 > 0 && i2<inMesh.elemX && i2 > 0)
	                                {
	                                    temp2 = j2*inMesh.elemX + i2; // neighbouring elem number
	                                    if(alpha[temp2] > control.aMin)
	                                    {
	                                        ftemp += inMesh.mat_vars[inMesh.mat_elems[temp2]];
	                                        k++;
	                                    }
	                                }
	                            }
	                        }
	                        if(k>0){ftemp /= (double)k; inMesh.mat_vars[i] = ftemp;}
	                    }
	                }
	            }
	            else
	            {
	                // just solve for material variable steps
	                for(i=0;i<=lsprob.num;i++)
	                {
	                    // shift constraint
	                    if(i>0) {

	                        // check max and min constraint change
	                        oMax=0.0; oMin=0.0;
	                        k=i*num_reduced;
	                        for(j=0;j<num_reduced;j++)
	                        {
	                            if(mat_sens_red[k] > 0.0){
	                                oMax += mat_sens_red[k]*mat_max[j];
	                                oMin += mat_sens_red[k]*mat_min[j];
	                            }
	                            else {
	                                oMax += mat_sens_red[k]*mat_min[j];
	                                oMin += mat_sens_red[k]*mat_max[j];
	                            }
	                            k++;
	                        }

	                        if(delcon[i-1] < 0.0) // if reduction in constraint
	                        {
	                            delcon[i-1] = (delcon[i-1] < oMin) ? oMin : delcon[i-1]; // set limit
	                        }
	                        else // if increase (or no change) in constraint
	                        {
	                            if(delcon[i-1] > oMax)
	                            {
	                                delcon[i-1] = oMax;
	                            }
	                            else
	                            {
	                            	// no-operation -- jeehanglee@gmail.com
	                                // delcon[i-1] = delcon[i-1];
	                            }
	                        }

	                        delcon[i-1] -= cblas_ddot(num_reduced, &mat_sens_red[num_reduced*i], 1, mat_min, 1);
	                    }

	                    // scale equations
	                    oMax = 0.0;
	                    for(k=0;k<num_reduced;k++)
	                    {
	                        j=k+(num_reduced*i);
	                        oMax = (fabs(mat_sens_red[j]) > oMax) ? fabs(mat_sens_red[j]) : oMax;
	                    }
	                    oMax = (oMax < 1.0e-20) ? 1.0 : oMax; // do not scale if max is 0

	                    for(k=0;k<num_reduced;k++)
	                    {
	                        j=k+(num_reduced*i);
	                        mat_sens_red[j] /= oMax; // scaled
	                    }

	                    if(i>0){ delcon[i-1] /= oMax; } // scale constraint value
	                }

	                for(i=0;i<num_reduced;i++)
	                {
	                    mat_max[i] -= mat_min[i]; // adjust so min is zero
	                }

	                // Can directly solve using LP sub-problem
	                csolver.LPsolve(num_reduced, lsprob.num, num_reduced, mat_step_red, mat_sens_red, &mat_sens_red[num_reduced], delcon, mat_max, 3);

	                for(i=0;i<num_reduced;i++)
	                {
	                    mat_step_red[i] += mat_min[i]; // re-adjust
	                }
	            }

	            // un-map reduced mat vars
	            j=0;
	            for(i=0;i<inMesh.NumDesMat;i++)
	            {
	                if(inc_mv[i]){ mat_step[i] = mat_step_red[j++]; }
	                else{ mat_step[i] = 0.0; }
	            }

	            free(mat_sens_red);
	            free(mat_step_red);
	            free(inc_mv);
	        }
	        else
	        {
	            // just solve for level-set Vnorm
	            temp = clevelset.SLPsubSol4(&inMesh, &levelset, &Bndry, cfl, delcon, lsprob.num, sens_ptr, Lbound, num_Lbound, Lbound_nums, lambda, active, Vnorm, pred_temp, 0, 0, 0, 0, 0, control.pinfo);
	        }
		}
		else
		{
			temp = clevelset.SLPsubSol4(&inMesh, &levelset, &Bndry, cfl, delcon, lsprob.num, sens_ptr, Lbound, num_Lbound, Lbound_nums, lambda, active, Vnorm, pred_temp, 0, 0, 0, 0, 0, control.pinfo);
		}

		// copy predicted obj & constraint changes
		k=(1+lsprob.num) * itt;
		for(i=0;i<=lsprob.num;i++)
		{
			pred[k++] = pred_temp[i];
		}

		// clear memory
		free(Lbound);
		free(Lbound_nums);
		free(Nsens);
		if(vol_sens){ free(vol_sens); }
		if(mass_sens){ free(mass_sens); }
		if(zero_sens){ free(zero_sens); }

		/*--------------------------------------------------------------------------------------/
		/																						/
		/		Calculate Extension Velocities, Gradients and Update the level set function		/
		/																						/
		/--------------------------------------------------------------------------------------*/

		// Fast marching method to determine Extension Velocities
		clevelset.Vext(&inMesh, &levelset, &Bndry, Vnorm);

		Grad = calloc(inMesh.NumNodes, sizeof(double)); // array for gradinet info

		// Calculate lsf gradients & update
		j2=inMesh.NodeY-1;
		i2=inMesh.NodeX-1;
		for(j=1;j<j2;j++)
		{
			for(i=1;i<i2;i++)
			{
				k = inMesh.Nodes2[i][j]; // read in node number
				if(levelset.active[k]) // If node is active (not fixed and within narrow band)
				{
					if(fabs(Vnorm[k]) > 1.0e-6) // If normal velocity not zero
					{
						// work out which sign to use for gradient calculation - upwind scheme
						temp = (Vnorm[k] < -0.0) ? -1 : 1;
						Grad[k] = clevelset.GradWENO(i, j, k, levelset.lsf, &inMesh, temp, Vnorm[k], 1.0);
					}
				}
			}
		}

		j = inMesh.NumNodes;
		for(k=0;k<j;k++)
		{
			if(levelset.active[k])
			{
				levelset.lsf[k] -=  Grad[k] * Vnorm[k]; // calculate updated lsf value
			}
		}

		if(control.pinfo == 3)	{
			// Write Node Gradient & Velocity information file
			coutput.OutHJVTK(&inMesh, Vnorm, Grad, itt, filename);
		}
		free(Grad);
		free(Vnorm);

		// update bar areas
		if(inMesh.bars)
		{
			for(j=0;j<inMesh.NumBars;j++)
			{
				inMesh.bar_areas[j] += bar_step[j];
			}
		}

		// update bc variable
		if(inMesh.des_bc)
		{
			for(j=0;j<inMesh.NumBC;j++)
			{
				inMesh.K_bc[j] += bc_step[j];
			}
		}

		// update material variables
		if(inMesh.des_mat && mat_opt)
		{
			for(j=0;j<inMesh.NumDesMat;j++)
			{
				inMesh.mat_vars[j] += mat_step[j];
			}
		}

		/*------------------------------------------------------------------/
		/																	/
		/		Re-initalise lsf if mine hit in Narrow Band					/
		/		or 20 iterations have passed since last re-initalization	/
		/																	/
		/------------------------------------------------------------------*/

		temp = 0; // flag to see if lsf need to be re-initalised
		if(ReCount < 19)
		{
			for(j=0;j<levelset.numMine;j++) // for all mine nodes
			{
				// if boundary is within 1 grid length of mine then re-initalise
				if( (fabs(levelset.lsf[levelset.mine[j]]) -  inMesh.h) < inMesh.tol )
				{
					temp = 1;
					printf("\n\nMine Hit! ..........Re-initalizing");
					break; // no need to keep hunting for mines
				}
			}
		}
		// set limit of 20 iterations before forcing a re-initialization
		if(temp == 0 && ReCount == 19)
		{
			printf("\n\n20 Iterations have passed! ..........Re-initalizing ");
			temp = 1;
		}
		else { ReCount++; }

		if(temp == 1)
		{
			ReCount = 0; // reset reinitialization count

			clevelset.ReInt(&inMesh, &levelset); // reinitialize the lsf
			clevelset.NarBand(&inMesh, &levelset, control.lband); // re-calculate the bandwidth

			printf("\nSigned distance function Re-initialized");
		}

		itt++; // next iteration

	} while (itt < control.maxItt); // automatically stop after max iterations

	if(itt == control.maxItt) {
		printf("\nSolution stopped after %i iterations\nObjective Value = %12.4e\n",itt,obj[itt-1]);
	}

	// Output the final ls function
	if(control.pinfo == 1) {
		coutput.OutPLotShapeVTK2(&inMesh, levelset.lsf, alpha, 1, itt, filename);
	}

	// Write Objetcive and Constraint convergence info file
	if(itt == control.maxItt){ itt--; }
	coutput.OutConv(itt, &lsprob, obj, cnstr, filename);

	// Write eigen frequency info file (is required)
	if(comp_eig==1){ coutput.OutFreq(itt, num_eig, save_freq, filename); }

	printf("\n\n------*---End of BLES Version 5.4 Program---*------\n\n");
	// return (0);
}
