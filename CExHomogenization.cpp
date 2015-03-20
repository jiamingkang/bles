/*
 * CExHomogenization.cpp
 *
 *  Created on: Dec 1, 2014
 *      Author: jeehang, Khalid Ismail
 */

#include <stdlib.h>
#include "CExHomogenization.h"

//
// Constructor / Destructor
//

CExHomogenization::CExHomogenization() {
	// TODO Auto-generated constructor stub
	m_pFileIn = NULL;
	memset(m_fileOut, 0x00, 128 * sizeof(char));
}

CExHomogenization::~CExHomogenization() {
	// TODO Auto-generated destructor stub
}

//
// Utilities
//

void CExHomogenization::_SetFilePath(char *path)
{
	int i = 0;
	int lenPath = strlen(path);

	if (path == NULL || lenPath < 1)
	{
		printf("File path is not specified. Process failed...\n\n");
		return;
	}

	// set input file path
	m_pFileIn = path;

	// set output file path
	for (i = 0; i < lenPath; i++)
    {
        if (path[i] == '.') 
		{
			// chop off file extension (if found)
			break;
		} 
        else if (path[i] == '\0')
		{
			break;
		}
        else 
		{
			m_fileOut[i] = path[i];
		}
    }
    m_fileOut[i] = '\0';
}

//
// implementation
//

//knai20@bath.ac.uk: OOD implementations
//void CExMinimiseCompliance::Initialise(char* arg, char* filename)
int CExHomogenization::InitialiseOOP(char *path)
{
	int res = -1;

	//Set file path - both input and output
	_SetFilePath(path);
	    
	// filename --> m_fileOut;
	// arg --> m_pFileIn;

	// read the input file
	res = m_input.read_input(m_pFileIn,
							&m_mesh,
							&numMat, 
							m_material,
							&m_levelset,
							&m_lsProblem,
							&m_control,
							&fixDof,
							&numCase, 
							&load, 
							&freeDof, 
							&m_lumpMass,
							&sw, 
							&acc);
	
	if (res == -1) // exit on error
	{
		printf("Failure in reading input file....\n");
		return 0;
	}
	
	m_mesh = m_input.GetMesh();

	if(m_control.m_outInfo==3)
	{
		m_output.OutNumber(m_mesh, filename);
		m_output.OutNodeCoord(m_mesh, filename);
	}

	// compute gauss point coords (used in sensitivity smoothing)
	gCoord = (Coord *) malloc(4 * m_mesh.m_numElem * sizeof(Coord));
	m_mathUtility.Gauss_Coord(m_mesh, gCoord);

    
	// calculate IN element stiffness matrix (& mass matrix)
	AreaElem = m_mesh.m_lenEdge * m_mesh.m_lenEdge; // Area of an element
	KE = (double **) malloc(numMat * sizeof(double*));
	ME = (double **) malloc(numMat * sizeof(double*));
	for (i = 0; i < numMat; i++)
	{
		KE[i] = (double *) malloc(KE_SIZE*sizeof(double));
		ME[i] = (double *) malloc(KE_SIZE*sizeof(double));
		
		// use 1.0 for thickness, as E & rho already multipled by thickness
		m_fem.KEMatrix(KE[i], m_material[i], m_mesh.m_thinkness);	// Element Stiffness Matrix for an IN element
		m_fem.MEMatrix(ME[i], m_material[i].m_rho, AreaElem, m_mesh.m_thinkness); // Element Mass Matrix for an IN element
	}

	/*------------------------------------------------------------------/
	/																	/
	/		Optimise structure boundary by the level set method			/
	/																	/
	/------------------------------------------------------------------*/

	alpha = (double *) malloc(m_mesh.m_numElem * sizeof(double)); // Array to store element areas
	
	sens_ptr = (double **) malloc( (1+m_lsProblem.m_numConstraint) * sizeof(double *) ); // pointer to sensitivity arrays for each fucntion (obj + constraints)

	active = (int *) malloc(m_lsProblem.m_numConstraint * sizeof(int)); // array for active constraints

	// Initialise length of AuxNodes & Boundary segment arrays - generous initial size
	m_boundary.m_pAuxNodes = (Coord *) malloc(m_mesh.m_numElem * sizeof(Coord));
	m_boundary.Bound = (Bseg *) malloc(m_mesh.m_numElem * sizeof(Bseg));
	m_boundary.m_pLenBseg = (double *)malloc(sizeof(double)); // initial memory allocation
	m_boundary.m_pWeightBseg = (double *) malloc(sizeof(double)); // initial memory allocation
	m_boundary.m_pConnAuxNode =(int *) malloc(sizeof(int)); // initial memory allocation
	m_boundary.m_pIndConnAuxNode = (int *) malloc( (m_mesh.m_numNodes+1) * sizeof(int)); // memory allocation
	lambda = (double *) calloc(1+m_lsProblem.m_numConstraint, sizeof(double)); // lambda values

	// stuff for Ku = f solve using MA57 solver (& eigenvalue solver)
	Kg.m_indRow=0; Kg.m_indCol=0; Kg.m_pMatEntry=0; Mg.m_indRow=0; Mg.m_indCol=0; Mg.m_pMatEntry=0;
    num_eig = (m_lsProblem.m_idObjType == 3) ? 1 : 0; // number of eigenvalues to compute
    comp_eig = (m_lsProblem.m_idObjType == 3) ? 1 : 0; // flag to decide if eigenvalues are to be computed
    order = (m_mesh.m_numNodes * NUM_DOF); // full matrix order
	numEnt = m_mesh.m_numElem * ((NUM_DOF*4)+1)*(NUM_DOF*2); // max number of index entries (symmetry!)

	// bar additonal design variables
	if(m_mesh.bars){numEnt += m_mesh.NumBars * 4;} // additional entries for bar elements
	if(m_mesh.bars){j=m_mesh.NumBars; bar_max= (double *) malloc(j*sizeof(double));
		bar_min= (double *) malloc(j*sizeof(double)); bar_step=(double *) malloc(j*sizeof(double));
		bar_sens = (double *) malloc(j * (1+m_lsProblem.m_numConstraint) * sizeof(double));}

	// designable bc addtional variables
	if(m_mesh.des_bc){numEnt += m_mesh.NumBC * 4 * NUM_DOF;} // additional entries for designable bcs
	if(m_mesh.des_bc){j=m_mesh.NumBC; bc_max= (double *) malloc(j*sizeof(double));
		bc_min=(double *) malloc(j*sizeof(double)); bc_step=(double *) malloc(j*sizeof(double));
		bc_sens =(double *) malloc(j * (1+m_lsProblem.m_numConstraint) * sizeof(double));}


	// designable material varibles
    mat_opt = m_mesh.dm_sim; // flag to switch between seq & sim opt
		if(m_mesh.des_mat){j=m_mesh.NumDesMat; mat_max=(double *) malloc(j*sizeof(double));
			mat_min=(double *) malloc(j*sizeof(double)); mat_step=(double *) malloc(j*sizeof(double));
			mat_sens = (double *) malloc(j * (1+m_lsProblem.m_numConstraint) * sizeof(double));}


	dispLen = order * numCase; // length of disp array (inc fixed dof)
    loadLen = freeDof * numCase; // length of load array (excld fixed dof)

	for(i=0;i<m_lsProblem.m_numConstraint;i++)
	{
		if(m_lsProblem.m_pConstraint[i].type==5){ num_adj++; } // displacement constraint
	}

	// compliance, 1st eigenvalue or mechanism objective
	if(m_lsProblem.m_idObjType==1 || m_lsProblem.m_idObjType==3 || m_lsProblem.m_idObjType==5){ num_sens++; }
	for(i=0;i<m_lsProblem.m_numConstraint;i++)
	{
		j = m_lsProblem.m_pConstraint[i].type;
		if(j==3 || j==5 || j==7 || j==8){ num_sens++; } // compliance or displacement constraint
		else if(j==6){num_sens++; num_eig++;} // eigenvalue ratio constraint
	}
	if(m_lsProblem.m_idObjType==5){num_sens++;}
	if(m_lsProblem.m_idObjType==3){numCase += num_eig;} // 1 primary state for each eigenvector


	if(m_lsProblem.m_idObjType==5){ prim_ptr = (double **) malloc(3 * sizeof(double*));
						dual_ptr = (double **) malloc(3 * sizeof(double*));
						wgt = (double *) malloc(sizeof(double));
						wgt[0] = 1.0; } // special cosideration for compliant mech
	else{ prim_ptr = (double **) malloc(numCase * sizeof(double*));
			dual_ptr = (double **) malloc(numCase * num_sens * sizeof(double*));
			wgt = (double *) malloc(numCase * sizeof(double));
			for(i=0;i<numCase;i++){wgt[i]=1.0;} } // set to 1.0 (for now)

	// compute load magnitudes for compliant mech design
	if(m_lsProblem.m_idObjType==5)
	{
		for(i=0;i<m_lsProblem.m_numConstraint;i++)
		{
			if(m_lsProblem.m_pConstraint[i].type==7)
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
	obj = (double *) malloc(m_control.m_maxIter*sizeof(double));
	cnstr = (double *) malloc(m_lsProblem.m_numConstraint*m_control.m_maxIter*sizeof(double));
	pred = (double *) malloc((1+m_lsProblem.m_numConstraint)*m_control.m_maxIter*sizeof(double)); // predicted obj & constraint change
	pred_temp = (double *) malloc((1+m_lsProblem.m_numConstraint)*sizeof(double)); // temp to feed SLPsubSol
	delcon = (double *) malloc(m_lsProblem.m_numConstraint*sizeof(double)); // array to store distance from constraint (constraint violation)
	relax = (double *) malloc(m_lsProblem.m_numConstraint*sizeof(double)); // relaxation factor
	if(comp_eig==1){save_freq = (double *) malloc(num_eig*2*m_control.m_maxIter*sizeof(double)); } // frequencies
	for(i=0;i<m_lsProblem.m_numConstraint;i++){relax[i] = 1.0;} // initialize to 1


    
    return m_control.m_maxIter;
    
}


int CExHomogenization::AnalyseOOP(int itt)
{


	//do {
		printf("\n\n------------\nStart iteration %i ...",itt);
		/*--------------------------------------------------------------/
		/																/
		/	   Analyse Structure, Solve Primary & Adjoint Equations		/
		/																/
		/--------------------------------------------------------------*/

		// Discretize the boundary & compute area ratios for AFG method
		m_mesh.Find_struct(&m_mesh, &m_levelset, &m_boundary, alpha, m_control.m_minArea);
		Ntot = m_mesh.m_numNodes + m_boundary.m_numAux; // total number of nodes

		if(m_control.m_outInfo > 1 || itt==0)	{
			// output lsf (& alpha) values
			m_output.OutPLotShapeVTK2(m_mesh, m_levelset.m_pNodalLsf, alpha, m_control.m_outInfo, itt, filename);
		}

		if(m_control.m_outInfo > 1) {
			printf("\nTotal nodes (inc. auxillary) = %i",Ntot);
		}

		// Assemble global matrices (in triplet form)

		// Initialize Index Arrays using dynamic memory allocation
		Kg.m_numEntry = numEnt; m_fem.set_sp_mat(&Kg);

		if(comp_eig==1)
		{
			Mg.m_numEntry = numEnt+m_lumpMass.m_numEntry; m_fem.set_sp_mat(&Mg); // mass matrix
			eig_vals = (double *) malloc(num_eig * 2 * sizeof(double)); // eigenvalues
			eig_vecs = (double *) malloc(num_eig * 2 * order * sizeof(double)); // eigenvectors
		}

		// calculate total structure volume (area) by summing element area array
		ftemp = 0.0;
		for(i=0;i<m_mesh.m_numElem;i++) {
			ftemp += alpha[i];
		}
		ftemp *= AreaElem;
		printf("\nTotal volume = %12.4e",ftemp);
		if(m_lsProblem.m_idObjType==2) // if volume is objective or constraint
		{
			obj[itt] = ftemp; // store value
		}
		else
		{
			for(i=0;i<m_lsProblem.m_numConstraint;i++){
				if(m_lsProblem.m_pConstraint[i].type==1){
					j = m_lsProblem.m_numConstraint * itt + i; // point to place in cnstr
					cnstr[j] = ftemp; // store volume as constraint value
				}
			}
		}

		// compute total mass (if required)
		for(i=0;i<m_lsProblem.m_numConstraint;i++)
		{
			if(m_lsProblem.m_pConstraint[i].type==2)
			{
				ftemp = 0.0;
				if(m_mesh.des_mat){temp=0; oMin = m_material[m_mesh.mat2].m_rho / m_material[m_mesh.mat1].m_rho;}
				for(j=0;j<m_mesh.m_numElem;j++) {
					oMax = alpha[j] * m_material[m_mesh.mat_type[j]].m_rho * m_mesh.m_thinkness; // area ratio x density
					if(m_mesh.des_mat && m_mesh.mat_elems[temp]==j) {
						//oMax *= m_mesh.mat_vars[temp] + (1.0-m_mesh.mat_vars[temp])*oMin; // modify for current variable
	                    oMax *= m_mesh.mat_vars[temp]*oMin + (1.0-m_mesh.mat_vars[temp]); // modify for current variable
						temp++; // next variable
					}
					ftemp += oMax;
				}
				ftemp *= m_mesh.m_lenEdge * m_mesh.m_lenEdge; // multiply by side length sqrd to get total mass
				printf("\nTotal mass (plane) = %12.4e",ftemp);

				// add bar mass (if there are some)
				if(m_mesh.bars)
				{
					oMax = 0.0;
					oMin = m_mesh.bar_mat->rho * m_mesh.m_lenEdge; // density x length
					for(j=0;j<m_mesh.NumBars;j++) {
						oMax += m_mesh.bar_areas[j] * oMin; // bar mass
					}
					printf("\nTotal mass (bar) = %12.4e",oMax);
					ftemp += oMax;
				}

				j = m_lsProblem.m_numConstraint * itt + i; // point to place in cnstr
				cnstr[j] = ftemp; // store mass as constraint value
				printf("\nTotal mass = %12.12e",ftemp);
			}
		}

	    // compute % of each material if using designable materials
	    if(m_mesh.des_mat)
	    {
	        temp = 0;
	        oMax = 0.0; // % material 2
	        oMin = 0.0; // sum alpha
	        for(j=0;j<m_mesh.m_numElem;j++)
	        {
	            if(m_mesh.des_mat && m_mesh.mat_elems[temp]==j)
	            {
	                oMin += alpha[j];
	                oMax += m_mesh.mat_vars[temp]*alpha[j]; // multiply by area ratio
	                temp++; // next variable
	            }
	        }
	        oMax /= oMin; // % material 2
	        printf("\nPcnt material 1 = %f",1.0-oMax);
	        printf("\nPcnt material 2 = %f",oMax);
	    }

		// asemble gloabl stiffness (& maybe mass) matrix
		m_fem.AFG_Matrix(comp_eig, KE, ME, &Kg, &Mg, &m_lumpMass, alpha, &m_mesh, m_material, m_control.m_minArea, m_control.m_minMass);

		// Solve the FE equation using the MA57 solver

		// solve Ku = f (for all cases)
		// remove fixed dof from global matrix
	    m_fem.rem_sp_mat(&Kg, fixDof, 1);

		if(m_lsProblem.m_idObjType != 3)
		{
	        disp = (double *) malloc(dispLen * sizeof(double));
	        // compute total load vector (in case of self-weight)
	        if(sw)
	        {
	            load_sw = (double *) malloc(dispLen * sizeof(double));
	            cmat.self_weight(&m_mesh, m_material, m_control.m_minArea, m_control.m_minMass, alpha, freeDof, fixDof, numCase, load, load_sw, acc);
	            cblas_dcopy(loadLen, load_sw, 1, disp, 1); // copy in load to send to FE_solve
	        }
	        else
	        {
	            cblas_dcopy(loadLen, load, 1, disp, 1); // copy in load to send to FE_solve
	        }
			m_solver.FE_Solve(&Kg, fixDof, disp, freeDof, order, numCase);
		}

		// compute eigenvalues (if required)
		if(comp_eig==1)
		{
			m_fem.rem_sp_mat(&Mg, fixDof, 1); // remove fixed dof from global matrix
			m_solver.eig_solve(num_eig, &Kg, &Mg, freeDof, order, fixDof, eig_vals, eig_vecs, m_control.m_outInfo);
			if(m_lsProblem.m_idObjType==3){ obj[itt] = sqrt(eig_vals[0]); } // 1st freq
			for(i=0;i<m_lsProblem.m_numConstraint;i++){
				if(m_lsProblem.m_pConstraint[i].type==6){
					j = m_lsProblem.m_numConstraint * itt + i; // point to place in cnstr
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
		if(m_lsProblem.m_idObjType==1){temp=1;}
		else {
			for(i=0;i<m_lsProblem.m_numConstraint;i++){ if(m_lsProblem.m_pConstraint[i].type==3){temp=1; break;} }
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

			if(m_lsProblem.m_idObjType==1){obj[itt]=obj_val;}
			for(i=0;i<m_lsProblem.m_numConstraint;i++){
				if(m_lsProblem.m_pConstraint[i].type==3){
					j = m_lsProblem.m_numConstraint * itt + i; // point to place in cnstr
					cnstr[j] = obj_val; // store compliance as constraint value
				}
			}
		}

	    if(sw){free(load_sw);} // no longer needed

		// calculate constraint values and adjoints for displacement constraints
		if(num_adj>0){adjont = (double *) calloc(order*num_adj, sizeof(double));}
		temp=0; // count adjoint load cases
		for(i=0;i<m_lsProblem.m_numConstraint;i++)
		{
			if(m_lsProblem.m_pConstraint[i].type==5)
			{
				j = (m_lsProblem.m_numConstraint * itt) + i; // point to place in cnstr
				k = (order * (int)m_lsProblem.m_pConstraint[i].data[2]) + (int)m_lsProblem.m_pConstraint[i].data[1]; // place in disp vector
				// compute constraint value
				ftemp = disp[k] - m_lsProblem.m_pConstraint[i].data[0];
				//printf("\nDisplacement constraint %i violated by %12.4e\n",i+1,ftemp);

				i2 = (int)(m_lsProblem.m_pConstraint[i].data[1]);
				k = (temp * freeDof) + fixDof[i2]; // place in adjoint load vector

				cnstr[j] = fabs(ftemp); // quicker than ^2 then sqrt!
				adjont[k] = ftemp/cnstr[j]; // adjoint load
				printf("\nAdjoint load %i = %12.4e\n",i+1,adjont[k]);
				temp++; // next adjoint
			}

			// compliant mechanism disp constraint
			else if(m_lsProblem.m_pConstraint[i].type==7)
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
				if(m_lsProblem.m_pConstraint[i].data[3] > 0.0){ ftemp += p2/m_lsProblem.m_pConstraint[i].data[3]; } // add spring stiffness
				ftemp = (p2*u21)/ftemp;

				// objective
				obj[itt] = ftemp / p1; // output / input force
				printf("\nMechanical advantage = %12.4e",obj[itt]);

				// input displacement
				j = (m_lsProblem.m_numConstraint * itt) + i; // point to place in cnstr
				if(m_lsProblem.m_pConstraint[i].data[3] > 0.0)
				{
					ftemp = u21 - (ftemp/m_lsProblem.m_pConstraint[i].data[3]);
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
			m_solver.FE_Solve(&Kg, fixDof, adjont, freeDof, order, num_adj);
			// this could be more efficient (by saving the factorization)
		}

		// clear memory for K (&M) matrices
		m_fem.free_sp_mat(&Kg); if(comp_eig==1){m_fem.free_sp_mat(&Mg);}

		if(m_control.m_outInfo == 3) {
			// write displacement & mode shape info file
			temp = (m_lsProblem.m_idObjType==3) ? 0 : numCase;
			m_output.OutDispVTK(m_mesh, temp, disp, 2*num_eig, eig_vecs, itt, filename);
		}
}

int CExHomogenization::SensitivityOOP(void)
{

		/*--------------------------------------------------------------/
		/																/
		/		Calculate shape sensitivities and velocity function		/
		/																/
		/--------------------------------------------------------------*/

		// Calculate boundary node sensitivities
		Nsens = (double *) calloc(Ntot*num_sens,sizeof(double));

		// point to correct senstivities etc ...

		// setup array of pointers to correct places in disp (for each case)
		if(m_lsProblem.m_idObjType != 3 && m_lsProblem.m_idObjType != 5)
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
		if(m_lsProblem.m_idObjType==1){
			sens_ptr[0] = &Nsens[temp++]; // compliance sens
			for(i=0;i<numCase;i++){dual_ptr[i*num_sens] = prim_ptr[i];} // compliance is self-adjoint
		}
		else if(m_lsProblem.m_idObjType==2){vol_sens = (double *) malloc(Ntot*sizeof(double));
								for(k=0;k<Ntot;k++){ vol_sens[k] = -1.0; }
								sens_ptr[0] = vol_sens;}   // volume sens
		else if(m_lsProblem.m_idObjType==3){
			sens_ptr[0] = &Nsens[temp++]; // 1st eigenvalue sensitivity
			prim_ptr[0] = eig_vecs; // eigenvector
			dual_ptr[0] = prim_ptr[0];
		}
		else if(m_lsProblem.m_idObjType==5) {
			sens_ptr[0] = Nsens; // first array of Nsens for objective
		}

		// constraints
		for(i=0;i<m_lsProblem.m_numConstraint;i++)
		{
			j = m_lsProblem.m_pConstraint[i].type;
			if(j==1)
			{
				vol_sens = (double *) malloc(Ntot*sizeof(double));
				for(k=0;k<Ntot;k++){ vol_sens[k] = -1.0; } // sensitvity for volume is -1 (geometric constraint)
				sens_ptr[i+1] = vol_sens;
			}
			else if(j==2)
			{
				mass_sens = (double *) calloc(Ntot,sizeof(double));

				if(m_mesh.des_mat)
				{
					i2=0; j2=0; // i2 = material element, j2 = boundary segment
					for(k=0;k<m_mesh.m_numElem;k++)
					{
						if(k==m_boundary.Bound[j2].e)
						{
							ftemp = m_mesh.mat_vars[i2];
							//ftemp = ftemp*m_material[m_mesh.mat1].rho + (1.0-ftemp)*m_material[m_mesh.mat2].rho;
	                        ftemp = (1.0-ftemp)*m_material[m_mesh.mat1].rho + ftemp*m_material[m_mesh.mat2].rho;
							ftemp *= -0.5 * m_mesh.m_thinkness;
							mass_sens[m_boundary.Bound[j2].n1] += ftemp;
							mass_sens[m_boundary.Bound[j2].n2] += ftemp;

							j2++; // update count

							if(j2 < m_boundary.m_numBound && k==m_boundary.Bound[j2].e) // 2 segments / elem
							{
								mass_sens[m_boundary.Bound[j2].n1] += ftemp;
								mass_sens[m_boundary.Bound[j2].n2] += ftemp;

								j2++; // update count
							}
						}

						if(k==m_mesh.mat_elems[temp]){i2++;} // update count of material element
					}
				}

				else
				{
					ftemp = -m_material[0].rho * m_mesh.m_thinkness; // mass senstivity (assumed) density x thickness
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
				i2 = (int)m_lsProblem.m_pConstraint[i].data[2]; // load case number
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
				zero_sens = (double *) calloc(Ntot,sizeof(double));
				sens_ptr[i+1] = zero_sens;
			}
		}

		if(m_lsProblem.m_idObjType == 3)
		{
			// compute boundary shape sensitivies (for eigenvalues)
			// MAT VAR EDIT
			m_sensitivity.AFG_Sens(&m_mesh, &m_boundary, alpha, m_material, Nsens, prim_ptr, dual_ptr, num_sens, num_eig, eig_vals, gCoord, m_control.m_minArea, 1,fact,sw,acc);
			for(i=0;i<m_lsProblem.m_numConstraint;i++)
			{
				if(m_lsProblem.m_pConstraint[i].type==6)
				{
					// process individual senstivities for constraint
					ftemp = cnstr[m_lsProblem.m_numConstraint * itt + i] / (2.0*eig_vals[0]);
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
		else if(m_lsProblem.m_idObjType == 5)
		{
			// combine senstivities together
			for(i=0;i<m_lsProblem.m_numConstraint;i++)
			{
				if(m_lsProblem.m_pConstraint[i].type==7){temp = i; break;} // displacement constraint
			}

			// compute reaction force
			ftemp = u22;
			if(m_lsProblem.m_pConstraint[temp].data[3] > 0.0){ ftemp += p2/m_lsProblem.m_pConstraint[temp].data[3];} // add spring stiffness

			fact[0] = (1.0/p1)*ftemp;
			fact[1] = -fact[0] * (u21 / ftemp);

			ftemp = (m_lsProblem.m_pConstraint[temp].data[3] > 0.0) ? p2/(ftemp*m_lsProblem.m_pConstraint[temp].data[3]) : 0.0;
			fact[2] = (u12/u22)*(ftemp-1.0)/p2;
			fact[2] += (u21/u22)*(ftemp-1.0)/p1;

			if(m_lsProblem.m_pConstraint[temp].data[3] > 0.0)
			{
				oMax = 1.0 - ftemp - (u22/( u22 + (p2/m_lsProblem.m_pConstraint[temp].data[3]) ))*ftemp;
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
			m_sensitivity.AFG_Sens(&m_mesh, &m_boundary, alpha, m_material, Nsens, prim_ptr, dual_ptr, 3, 1, wgt, gCoord, m_control.m_minArea, 2,fact,sw,acc);
		}
		else
		{
			// compute boundary shape sensitivies (for compliance or displacement functions)
			m_sensitivity.AFG_Sens(&m_mesh, &m_boundary, alpha, m_material, Nsens, prim_ptr, dual_ptr, num_sens, numCase, wgt, gCoord, m_control.m_minArea, 0,fact,sw,acc);
			free(disp);
		}

		// bar senstivities
		if(m_mesh.bars && m_lsProblem.m_idObjType == 1)
		{
			// call routine for computing bar compliance senstivites
			// inlcuding multi-load cases
			for(j=0;j<m_mesh.NumBars;j++) {
				bar_sens[j] = 0.0; // zero initial
			}

			m_sensitivity.barSens(&m_mesh, bar_sens, prim_ptr, numCase, wgt);

			// mass senstivites
			for(i=0;i<m_lsProblem.m_numConstraint;i++)
			{
				if(m_lsProblem.m_pConstraint[i].type==2)
				{
					ftemp = m_mesh.bar_mat->rho * m_mesh.m_lenEdge; // bar mass senstivity (wrt area)
					k = m_mesh.NumBars * (i+1);
					for(j=0;j<m_mesh.NumBars;j++) {
						bar_sens[k++] = ftemp; // mass constraint
					}
				}
			}

			// also set move limits on bar areas here
			ftemp = 0.05 * (m_mesh.bar_max - m_mesh.bar_min);
			for(j=0;j<m_mesh.NumBars;j++)
			{
				oMax = m_mesh.bar_areas[j] + ftemp;
				oMin = m_mesh.bar_areas[j] - ftemp;
				bar_max[j] = (oMax < m_mesh.bar_max) ? ftemp : m_mesh.bar_max - m_mesh.bar_areas[j];
				bar_min[j] = (oMin > m_mesh.bar_min) ? -ftemp : m_mesh.bar_min - m_mesh.bar_areas[j];
			}

			// output bar info
			m_output.OutBars(m_mesh, (m_lsProblem.m_numConstraint+1), bar_sens, m_control.m_outInfo, itt, filename);
		}

		// additional bc sensitivities
		if(m_mesh.des_bc && m_lsProblem.m_idObjType == 1)
		{
			for(j=0;j<m_mesh.NumBC;j++) {
				bc_sens[j] = 0.0; // zero initial
			}

			// CALL SENSITIVITY FUNCTION
			m_sensitivity.bcSens(&m_mesh, bc_sens, prim_ptr, numCase, wgt);

			// Sensitivity filtering
			if(itt==0)
			{
				// analyse connectivity of des bc elements
				int xMin, xMax, yMin, yMax, n, m;
				temp = m_mesh.NumBC;
				bc_con = (int *) malloc(8*temp*sizeof(int));
				bc_dist = (double *) malloc(8*temp*sizeof(double));
				bc_con_ind = (int *) malloc((temp+1)*sizeof(int));
				int count = 0;
				for(i=0;i<temp;i++)
				{
					// compute element indices
					ftemp = (double)m_mesh.BC_nums[i] / (double)m_mesh.m_elemX;
					m = floor(ftemp);
					n = m_mesh.BC_nums[i] - m*m_mesh.m_elemX;
					bc_con_ind[i] = count;

					// search for neighbouring elements with designable bc
					xMin = n-1; xMax=n+2;
					yMin = m-1; yMax=m+2; // search limits
					xMin = (xMin < 0) ? 0 : xMin;
					yMin = (yMin < 0) ? 0 : yMin;
					xMax = (xMax > m_mesh.m_elemX) ? m_mesh.m_elemX : xMax;
					yMax = (yMax > m_mesh.m_elemY) ? m_mesh.m_elemY : yMax; // adjust search limits

					for(j2=yMin;j2<yMax;j2++)
					{
						for(i2=xMin;i2<xMax;i2++)
						{
							temp2 = j2*m_mesh.m_elemX + i2; // elem number
							for(j=0;j<temp;j++)
							{
								if(temp2 == m_mesh.BC_nums[j] && j != i)
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
				bc_con = (int *) realloc(bc_con, count * sizeof(int));
				bc_dist = (double *) realloc(bc_dist, count * sizeof(double));
			}

			// Filter!
			{
				temp = m_mesh.NumBC;
				double *bc_temp = (double *) malloc(temp*sizeof(double));
				cblas_dcopy(temp, bc_sens, 1, bc_temp, 1); // copy original

				for(i=0;i<temp;i++) // for all des bcs
				{
					temp2 = bc_con_ind[i+1];
					oMax = 1.5*m_mesh.K_bc[i]*bc_temp[i];
					oMin = 1.5;
					for(j=bc_con_ind[i];j<temp2;j++) // for all connected elems
					{
						k = bc_con[j]; // connectd bc number
						oMax += bc_dist[j] * m_mesh.K_bc[k] * bc_temp[k];
						oMin += ftemp;
					}

					bc_sens[i] = oMax / (oMin * m_mesh.K_bc[i]); // filtered sensitivity
				}

				free(bc_temp);
			}

			// scale maximum sensitiviy to 1
			oMax = 0.0;
			for(j=0;j<m_mesh.NumBC;j++) {
				ftemp = fabs(bc_sens[j]);
				oMax = (ftemp > oMax) ? ftemp : oMax;
			}
			for(j=0;j<m_mesh.NumBC;j++) {
				bc_sens[j] /= oMax;
			}

			// other senstivities
			for(i=0;i<m_lsProblem.m_numConstraint;i++)
			{
				// constraint cost sensitivities are 1.0
				// also compute the constraint value
				if(m_lsProblem.m_pConstraint[i].type==10)
				{
					ftemp = 0.0;
					k = m_mesh.NumBC;
					for(j=0;j<m_mesh.NumBC;j++) {
						bc_sens[k++] = 1.0;
						ftemp += m_mesh.K_bc[j];
					}

					j = (m_lsProblem.m_numConstraint * itt) + i; // point to place in cnstr
					cnstr[j] = ftemp;
				}
			}

			// also set move limits on bc variables here
			for(j=0;j<m_mesh.NumBC;j++)
			{
				oMax = m_mesh.K_bc[j] + 0.05;
				oMin = m_mesh.K_bc[j] - 0.05;
				bc_max[j] = (oMax < 1.0) ? 0.05 : 1.0 - m_mesh.K_bc[j];
				bc_min[j] = (oMin > 1.0e-4) ? -0.05 : 1.0e-4 - m_mesh.K_bc[j];
			}

			// output designable BC info
			m_output.OutDesBC(m_mesh, bc_sens, m_control.m_outInfo, itt, filename);
		}

		// material sensitvities
		if(m_mesh.des_mat)
		{
	        double *mat_sens_temp = (double *) calloc(numCase * m_mesh.NumDesMat, sizeof(double));

	        // call routine for computing material complinace senstivites
	        if(m_lsProblem.m_idObjType == 1)
	        {
	            // compliance objective sensitivity
	            if(m_mesh.mat_lin){ m_sensitivity.matSens_comp(&m_mesh, m_material, KE[0], mat_sens_temp, numCase, wgt, prim_ptr, dual_ptr, alpha, m_control.m_minArea, sw, acc); }
	            //else{ HS_Sens_eig(&m_mesh, m_material, KE[0], ME[0], mat_sens_temp, num_eig, eig_vals, prim_ptr, alpha, m_control.m_minArea); }
	            cblas_dcopy(m_mesh.NumDesMat, mat_sens_temp, 1, mat_sens, 1); // copy objective sensitivities
	        }

			// call routine for computing material eig senstivites
	        if(m_lsProblem.m_idObjType == 3)
	        {
	            // eigenvalue objective sensitivity
	            if(m_mesh.mat_lin){ m_sensitivity.matSens_eig(&m_mesh, m_material, KE[0], ME[0], mat_sens_temp, num_eig, eig_vals, prim_ptr, alpha, m_control.m_minArea); }
	            else{ m_sensitivity.HS_Sens_eig(&m_mesh, m_material, KE[0], ME[0], mat_sens_temp, num_eig, eig_vals, prim_ptr, alpha, m_control.m_minArea); }
	            cblas_dcopy(m_mesh.NumDesMat, mat_sens_temp, 1, mat_sens, 1); // copy objective sensitivities
	        }

	        // constraint senstivites
	        for(i=0;i<m_lsProblem.m_numConstraint;i++)
	        {
	            // volume constraint sensitivities
	            if(m_lsProblem.m_pConstraint[i].type==1)
	            {
	                k = m_mesh.NumDesMat * (i+1);
	                for(j=0;j<m_mesh.NumDesMat;j++) {
	                    mat_sens[k++] = 0.0; // material choice does not effect volume
	                }
	            }

	            else if(m_lsProblem.m_pConstraint[i].type==2)
	            {
	                // material mass senstivity
	                ftemp = (m_material[m_mesh.mat2].rho - m_material[m_mesh.mat1].rho) * m_mesh.m_lenEdge * m_mesh.m_lenEdge * m_mesh.m_thinkness;
	                k = m_mesh.NumDesMat * (i+1);
	                for(j=0;j<m_mesh.NumDesMat;j++) {
	                    mat_sens[k++] = alpha[m_mesh.mat_elems[j]] * ftemp; // mass sensitivity
	                }
	            }

	            // freq ratio sensitivites
	            else if(m_lsProblem.m_pConstraint[i].type==6)
	            {
	                // process individual senstivities for constraint
	                ftemp = cnstr[m_lsProblem.m_numConstraint * itt + i] / (2.0*eig_vals[0]);
	                oMax = 1.0 / (2.0*sqrt(eig_vals[0]*eig_vals[1]));
	                //oMin = 1.0 / (2.0*sqrt(eig_vals[0]));  // factors
	                k = m_mesh.NumDesMat * (i+1); // place in mat_sens for freq ratio sens
	                temp = m_mesh.NumDesMat; // start of 2nd eigenvalue sens
	                for(j=0;j<m_mesh.NumDesMat;j++)
	                {
	                    mat_sens[k++] = oMax*mat_sens_temp[temp++] - ftemp*mat_sens_temp[j]; // swap sign for reduction
	                    //Nsens[j] *= oMin; // do not scale 1st eigenvalue senstivity
	                }
	            }
	        }

	        // find max for FD check
	        /*oMax=0.0;
	        j=-1;
	        for(i=0;i<m_mesh.NumDesMat;i++)
	        {
	            ftemp = fabs(mat_sens[i]);
	            if(ftemp > oMax){ oMax=ftemp; j=i;}
	        }

	        printf("\nMax obj sens %i = %12.12e",j,mat_sens[j]);
	        for(i=0;i<m_lsProblem.m_numConstraint;i++)
	        {
	            k = m_mesh.NumDesMat * (i+1);
	            printf("\nConst sens %i = %12.12e",i,mat_sens[k+j]);
	        }*/

	        free(mat_sens_temp);

			// also set move limits on material variables
			ftemp = 0.02;
			for(j=0;j<m_mesh.NumDesMat;j++)
			{
				oMax = m_mesh.mat_vars[j] + ftemp;
				oMin = m_mesh.mat_vars[j] - ftemp;
				mat_max[j] = (oMax < 1.0) ? ftemp : 1.0 - m_mesh.mat_vars[j];
				mat_min[j] = (oMin > 0.0) ? -ftemp : 0.0 - m_mesh.mat_vars[j];
			}

			// output material info
			m_output.OutDesMat(m_mesh, alpha, m_control.m_minArea, (m_lsProblem.m_numConstraint+1), mat_sens, m_control.m_outInfo, itt, filename);
		}

		// free required memory
		if(comp_eig==1) { free(eig_vals); free(eig_vecs); }
		if(num_adj>0){ free(adjont); }

		if(m_control.m_outInfo == 3)	{
			// Write node sensitvity information file (as a boundary mesh)
			m_output.OutBoundVTK(m_mesh, &m_boundary, (m_lsProblem.m_numConstraint+1), sens_ptr, itt, filename);
		}

		// Check convergence criteria
		// calculate difference between current and target constraint
		temp = 1; // check for violated constraints
		j = m_lsProblem.m_numConstraint;
		for(i=0;i<j;i++)
		{
			k = j*itt + i; // pointer for cnstr
			ftemp = cnstr[k]; // current constraint value
			if(m_lsProblem.m_pConstraint[i].type==1)
			{
				delcon[i] = m_lsProblem.m_pConstraint[i].data[0] - ftemp;
				printf("\nVolume constraint violated by %12.4e\n",delcon[i]);
				if(delcon[i]/m_lsProblem.m_pConstraint[i].data[0] < -m_control.gm){temp=0;}
			}
			else if(m_lsProblem.m_pConstraint[i].type==2)
			{
				delcon[i] = m_lsProblem.m_pConstraint[i].data[0] - ftemp;
				printf("\nMass constraint violated by %12.4e\n",delcon[i]);
				if(delcon[i]/m_lsProblem.m_pConstraint[i].data[0] < -m_control.gm){temp=0;}
			}
			else if(m_lsProblem.m_pConstraint[i].type==3)
			{
				delcon[i] = m_lsProblem.m_pConstraint[i].data[0] - ftemp;
				printf("\nCompliance constraint violated by %12.4e\n",delcon[i]);
				if(delcon[i] < 0.0){temp=0;}
			}
			else if(m_lsProblem.m_pConstraint[i].type==5)
			{
				delcon[i] = -ftemp;
				printf("\nDisplacement constraint %i violated by %12.4e\n",i+1,delcon[i]);
				if(fabs(delcon[i]/m_lsProblem.m_pConstraint[i].data[0]) > 0.01){temp=0;} // within 1%
			}
			else if(m_lsProblem.m_pConstraint[i].type==6)
			{
				delcon[i] = ftemp - m_lsProblem.m_pConstraint[i].data[0]; // reverse sign (as greater than constraint)
				printf("\nEig freq ratio constraint %i violated by %12.4e\n",i+1,-delcon[i]);
				if(delcon[i] < 0.0){temp=0;}
			}
			else if(m_lsProblem.m_pConstraint[i].type==7)
			{
				delcon[i] = m_lsProblem.m_pConstraint[i].data[0] - ftemp;
				printf("\nInput disp constraint %i violated by %12.4e\n",i+1,delcon[i]);
				if(delcon[i] < 0.0){temp=0;}
			}
			else if(m_lsProblem.m_pConstraint[i].type==8)
			{
				delcon[i] = m_lsProblem.m_pConstraint[i].data[0] - ftemp;
				printf("\nOutput comp constraint %i violated by %12.4e\n",i+1,delcon[i]);
				if(delcon[i] < 0.0){temp=0;}
			}
			else if(m_lsProblem.m_pConstraint[i].type==10)
			{
				delcon[i] = m_lsProblem.m_pConstraint[i].data[0] - ftemp;
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
			if(ftemp < m_control.gm) {
	            if(m_mesh.des_mat && !m_mesh.dm_sim && !mat_opt)
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
			for(i=0;i<m_lsProblem.m_numConstraint;i++)
			{
				if(active[i] != 0)
				{
					// multiply target constraint change by: actual / predicted from last itt
					j=(1+m_lsProblem.m_numConstraint)*(itt-1) + i+1; // place in predicted (last itt)

					ftemp = cnstr[itt*m_lsProblem.m_numConstraint + i];
					ftemp -= cnstr[(itt-1)*m_lsProblem.m_numConstraint + i];  // actual constraint change
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
    return 0;
}

int CExHomogenization::OptimiseOOP(void)
{

		// Compute boundary intergral of shape senstivities
		// then use LP sub-problem to obtain the velocity function
		temp = m_lsProblem.m_numConstraint + 1; // number of "functions"
		Vnorm = (double *) calloc(Ntot,sizeof(double));		// normal velocity array
		Lbound = (double *) malloc(Ntot*temp*sizeof(double));	// boundary intergral data array
		Lbound_nums = (int *) malloc(Ntot*sizeof(int)); // global boundary point numbers

		// compute variabe weights for boundary intergral
		cboundary.BsegWgt(&m_boundary, &m_mesh);

		// compute coefficients for discretized boundary integral
		cboundary.BoundInt(&m_mesh, &m_levelset, &m_boundary, temp, sens_ptr, Lbound_nums, &num_Lbound, Lbound);

		if(m_control.m_outInfo == 3)	{
			// Write boundary integration information file
			m_output.OutBoundInt((m_lsProblem.m_numConstraint+1), num_Lbound, Lbound_nums, Lbound, itt, filename);
		}

		// realloc & free memory
		Lbound = (double *) realloc(Lbound, num_Lbound*temp*sizeof(double));
		Lbound_nums = (int *) realloc(Lbound_nums, num_Lbound*sizeof(int));

		// obtian Vnorm using SLP subproblem
		if(m_mesh.bars)
		{
			temp = m_levelset.SLPsubSol4(&m_mesh, &m_levelset, &m_boundary, cfl, delcon, m_lsProblem.m_numConstraint, sens_ptr, Lbound, num_Lbound, Lbound_nums, lambda, active, Vnorm, pred_temp, m_mesh.NumBars, bar_sens, bar_min, bar_max, bar_step, m_control.m_outInfo);
		}
		else if(m_mesh.des_bc)
		{
			// solve for level set velocity
			temp = m_levelset.SLPsubSol4(&m_mesh, &m_levelset, &m_boundary, cfl, delcon, m_lsProblem.m_numConstraint, sens_ptr, Lbound, num_Lbound, Lbound_nums, lambda, active, Vnorm, pred_temp, 0, 0, 0, 0, 0, m_control.m_outInfo);

			// designable boundary condition cost constraint
			// shift for LP sub-solve
			for(i=0;i<m_lsProblem.m_numConstraint;i++)
			{
				if(m_lsProblem.m_pConstraint[i].type==10)
				{
					j = (m_lsProblem.m_numConstraint * itt) + i; // point to place in cnstr
					ftemp = m_lsProblem.m_pConstraint[i].data[0] - cnstr[j];
					ftemp -= cblas_ddot(m_mesh.NumBC, &bc_sens[m_mesh.NumBC], 1, bc_min, 1);
					break;
				}
			}

			for(i=0;i<m_mesh.NumBC;i++)
			{
				bc_max[i] -= bc_min[i]; // adjust so min is zero
			}

			// Can directly solve using LP sub-problem
			m_solver.LPsolve(m_mesh.NumBC, 1, m_mesh.NumBC, bc_step, bc_sens, &bc_sens[m_mesh.NumBC], &ftemp, bc_max, m_control.m_outInfo);

			for(i=0;i<m_mesh.NumBC;i++)
			{
				bc_step[i] += bc_min[i]; // re-adjust
			}
		}
		else if(m_mesh.des_mat)
		{
	        if(m_mesh.dm_sim || mat_opt)
	        {
	            // reduce number of mat vars using a map
	            bool *inc_mv = (bool *) malloc(m_mesh.NumDesMat*sizeof(bool));
	            int num_reduced=0;
	            for(i=0;i<m_mesh.NumDesMat;i++)
	            {
	                j = m_mesh.mat_elems[i]; // elem number
	                if(alpha[j] > m_control.m_minArea){inc_mv[i] = 1; num_reduced++;}
	                else{inc_mv[i]=0;}
	            }

	            double *mat_sens_red = (double *) malloc(num_reduced*(1+m_lsProblem.m_numConstraint)*sizeof(double));
	            double *mat_step_red = (double *) malloc(num_reduced*sizeof(double));
	            j=0;
	            for(i=0;i<m_mesh.NumDesMat;i++)
	            {
	                // now copy sensitivity & limit data
	                if(inc_mv[i])
	                {
	                    mat_max[j] = mat_max[i];
	                    mat_min[j] = mat_min[i];

	                    for(k=0;k<=m_lsProblem.m_numConstraint;k++)
	                    {
	                        j2 = k*m_mesh.NumDesMat + i; // original
	                        i2 = k*num_reduced + j; // reduced
	                        mat_sens_red[i2] = mat_sens[j2]; // copy
	                    }
	                    j++; // next reduced entry
	                }
	            }

	            j = num_reduced * (m_lsProblem.m_numConstraint+1);

	            //temp = SLPsubSol4(&m_mesh, cfl, delcon, m_lsProblem.m_numConstraint, sens_ptr, Lbound, num_Lbound, Lbound_nums, lambda, active, m_boundary.AuxNodes, Vnorm, pred_temp, m_mesh.NumDesMat, mat_sens, mat_min, mat_max, mat_step, m_control.m_outInfo);

	            if(m_mesh.dm_sim)
	            {
	                // solve for simultaneous update
	                temp = m_levelset.SLPsubSol4(&m_mesh, &m_levelset, &m_boundary, cfl, delcon, m_lsProblem.m_numConstraint, sens_ptr, Lbound, num_Lbound, Lbound_nums, lambda, active, Vnorm, pred_temp, num_reduced, mat_sens_red, mat_min, mat_max, mat_step_red, m_control.m_outInfo);

	                // also ensure smooth transition for materials in elems just out from boundary
	                int n, m;
	                for(i=0;i<m_mesh.NumDesMat;i++)
	                {
	                    // if not included
	                    if(!inc_mv[i])
	                    {
	                        // compute element indices
	                        ftemp = (double)m_mesh.mat_elems[i] / (double)m_mesh.m_elemX;
	                        m = floor(ftemp);
	                        n = m_mesh.mat_elems[i] - m*m_mesh.m_elemX;

	                        // average neighbours that are part of structure
	                        k=0;
	                        ftemp=0.0;
	                        for(j2=m-1;j2<m+1;j2++)
	                        {
	                            for(i2=n-1;i2<n+1;i2++)
	                            {
	                                if(j2<m_mesh.m_elemY && j2 > 0 && i2<m_mesh.m_elemX && i2 > 0)
	                                {
	                                    temp2 = j2*m_mesh.m_elemX + i2; // neighbouring elem number
	                                    if(alpha[temp2] > m_control.m_minArea)
	                                    {
	                                        ftemp += m_mesh.mat_vars[m_mesh.mat_elems[temp2]];
	                                        k++;
	                                    }
	                                }
	                            }
	                        }
	                        if(k>0){ftemp /= (double)k; m_mesh.mat_vars[i] = ftemp;}
	                    }
	                }
	            }
	            else
	            {
	                // just solve for material variable steps
	                for(i=0;i<=m_lsProblem.m_numConstraint;i++)
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
	                m_solver.LPsolve(num_reduced, m_lsProblem.m_numConstraint, num_reduced, mat_step_red, mat_sens_red, &mat_sens_red[num_reduced], delcon, mat_max, 3);

	                for(i=0;i<num_reduced;i++)
	                {
	                    mat_step_red[i] += mat_min[i]; // re-adjust
	                }
	            }

	            // un-map reduced mat vars
	            j=0;
	            for(i=0;i<m_mesh.NumDesMat;i++)
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
	            temp = m_levelset.SLPsubSol4(&m_mesh, &m_levelset, &m_boundary, cfl, delcon, m_lsProblem.m_numConstraint, sens_ptr, Lbound, num_Lbound, Lbound_nums, lambda, active, Vnorm, pred_temp, 0, 0, 0, 0, 0, m_control.m_outInfo);
	        }
		}
		else
		{
			temp = m_levelset.SLPsubSol4(&m_mesh, &m_levelset, &m_boundary, cfl, delcon, m_lsProblem.m_numConstraint, sens_ptr, Lbound, num_Lbound, Lbound_nums, lambda, active, Vnorm, pred_temp, 0, 0, 0, 0, 0, m_control.m_outInfo);
		}

		// copy predicted obj & constraint changes
		k=(1+m_lsProblem.m_numConstraint) * itt;
		for(i=0;i<=m_lsProblem.m_numConstraint;i++)
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
		m_levelset.Vext(&m_mesh, &m_levelset, &m_boundary, Vnorm);

		Grad = (double *) calloc(m_mesh.m_numNodes, sizeof(double)); // array for gradinet info

		// Calculate lsf gradients & update
		j2=m_mesh.NodeY-1;
		i2=m_mesh.NodeX-1;
		for(j=1;j<j2;j++)
		{
			for(i=1;i<i2;i++)
			{
				k = m_mesh.Nodes2[i][j]; // read in node number
				if(m_levelset.active[k]) // If node is active (not fixed and within narrow band)
				{
					if(fabs(Vnorm[k]) > 1.0e-6) // If normal velocity not zero
					{
						// work out which sign to use for gradient calculation - upwind scheme
						temp = (Vnorm[k] < -0.0) ? -1 : 1;
						Grad[k] = m_levelset.GradWENO(i, j, k, m_levelset.m_pNodalLsf, &m_mesh, temp, Vnorm[k], 1.0);
					}
				}
			}
		}

		j = m_mesh.m_numNodes;
		for(k=0;k<j;k++)
		{
			if(m_levelset.active[k])
			{
				m_levelset.m_pNodalLsf[k] -=  Grad[k] * Vnorm[k]; // calculate updated lsf value
			}
		}

		if(m_control.m_outInfo == 3)	{
			// Write Node Gradient & Velocity information file
			m_output.OutHJVTK(m_mesh, Vnorm, Grad, itt, filename);
		}
		free(Grad);
		free(Vnorm);

		// update bar areas
		if(m_mesh.bars)
		{
			for(j=0;j<m_mesh.NumBars;j++)
			{
				m_mesh.bar_areas[j] += bar_step[j];
			}
		}

		// update bc variable
		if(m_mesh.des_bc)
		{
			for(j=0;j<m_mesh.NumBC;j++)
			{
				m_mesh.K_bc[j] += bc_step[j];
			}
		}

		// update material variables
		if(m_mesh.des_mat && mat_opt)
		{
			for(j=0;j<m_mesh.NumDesMat;j++)
			{
				m_mesh.mat_vars[j] += mat_step[j];
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
			for(j=0;j<m_levelset.numMine;j++) // for all mine nodes
			{
				// if boundary is within 1 grid length of mine then re-initalise
				if( (fabs(m_levelset.m_pNodalLsf[m_levelset.mine[j]]) -  m_mesh.m_lenEdge) < m_mesh.m_tolerance )
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

			m_levelset.ReInt(&m_mesh, &m_levelset); // reinitialize the lsf
			m_levelset.NarBand(&m_mesh, &m_levelset, m_control.m_lband); // re-calculate the bandwidth

			printf("\nSigned distance function Re-initialized");
		}

		itt++; // next iteration

	} while (itt < m_control.m_maxIter); // automatically stop after max iterations

void CExHomogenization::OutputOOP(void)
{
	if(itt == m_control.m_maxIter) {
		printf("\nSolution stopped after %i iterations\nObjective Value = %12.4e\n",itt,obj[itt-1]);
	}

	// Output the final ls function
	if(m_control.m_outInfo == 1) {
		m_output.OutPLotShapeVTK2(m_mesh, m_levelset.m_pNodalLsf, alpha, 1, itt, filename);
	}

	// Write Objective and Constraint convergence info file
	if(itt == m_control.m_maxIter){ itt--; }
	m_output.OutConv(itt, &m_lsProblem, obj, cnstr, filename);

	// Write eigen frequency info file (is required)
	if(comp_eig==1){ m_output.OutFreq(itt, num_eig, save_freq, filename); }

	printf("\n\n------*---End of BLES Version 5.4 Program---*------\n\n");
	// return (0);
}



//knai20@bath.ac.uk: Non OOD implementations


