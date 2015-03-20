/*
 * CEx3DMinimiseCompliance.cpp
 *
 *  Created on: 20 Mar 2015
 *      Author: knai20-admin
 */

#include "CEx3DMinimiseCompliance.h"

CEx3DMinimiseCompliance::CEx3DMinimiseCompliance() {
	// TODO Auto-generated constructor stub

}

CEx3DMinimiseCompliance::~CEx3DMinimiseCompliance() {
	// TODO Auto-generated destructor stub
}


void CEx3DMinimiseCompliance::Execute()
{

	// timing variables and file
	struct  timeval clck1, clck2; // timing variables
    gettimeofday(&clck1,NULL);
    cout << "\n\nStart of Level-set 3D v2.3";
    //int numt = omp_get_max_threads(); // get number of thread
    //cout << "\nMax threads = " << numt;
    //omp_set_num_threads(28);

	int i,j,k,ind,ind2,err,end; // incrementors etc.
    char foname[200];
    int itt;
    int itt_in=0; // start iteration number

    // input files
    char *FE_file = 0;
    char *Opt_file = 0;
    char *restart = 0;

    // read input file names
    err = read_main_file(argv[1], &FE_file, 0, &Opt_file, 0, &restart, 0, &itt_in);
    if(err==1){return 0;}

    if(FE_file==0){ std::cout << "\nERROR! No FE mesh input file - abort"; return 0;}
    if(Opt_file==0){ std::cout << "\nERROR! No Optimization input file - abort"; return 0;}

    // prefix name for output files
    char infile[100];
    for(i=0;i<99;i++)
    {
        if(argv[1][i] == '.') {break;}
        else if(argv[1][i] == '\0') {break;}
        else{infile[i] = argv[1][i];}
    }
    infile[i] = '\0';

    // read controls & optimization problem definition
    cntrl control;
    prob lsprob; lsprob.obj = -1;
    int numMission = 1; // generalize to multiple load cases
    err = read_control(Opt_file, &control, &lsprob, &numMission);
    if(err==1) {return 0;}

    double *wgt_case = new double [numMission]; // array of mission weight numbers
	for(i=0;i<numMission;i++){ wgt_case[i]=1.0;}

    // check to see if a mass matrix is required (for eigenvalue freq analysis)
    if(lsprob.obj == 3 || lsprob.obj == 4){ control.comp_mass = true; }
    else
    {
        for(i=0;i<lsprob.num;i++)
        {
            if(lsprob.con[i].type==4){ control.comp_mass = true; break; }
        }
    }

    // check start iteration (if specified a restart)
    if(itt_in >= control.maxItt)
    {
        cout << "\nError! Specified start iteration (" << itt_in << ") > specified max (" << control.maxItt << ") Abort";
        return 0;
    }

    // variables used to read in mesh data etc...
	int numElem, numNode;
    int numElset, numNset;
    int numMat;
	int *node_nums_temp, *elem_def_temp, *elem_nums_temp, *elem_type;
	Coord3D *node_crds_temp;

	node_nums_temp = new int [MAX_NODES];
	node_crds_temp = new Coord3D [MAX_NODES];
    elem_nums_temp = new int [MAX_ELEMS];
    elem_type      = (int*)malloc(MAX_ELEMS * sizeof(int));
	elem_def_temp  = new int [MAX_ELEMS*8];

    // -------------------------------------------------------------------- //
    //                                                                      //
    //      Step 1: Read in structural input file                           //
    //              Make element stiffness (& mass) matrices                //
    //                                                                      //
    // -------------------------------------------------------------------- //

    // function to read in an Abaqus style input file of nodes and elements
    err = read_abq_mesh(FE_file, &numNode, &numElem, node_nums_temp, node_crds_temp, elem_nums_temp, elem_def_temp, elem_type);
    if(err==1) {return 0;}

	// create mesh (elements and nodes)
	mesh WingMesh(numNode,numElem);
    elem_type = (int*)realloc(elem_type, numElem*sizeof(int)); // resize elememt type array
    WingMesh.el_type = elem_type; // add element type array to mesh

	// copy node coords
	for(i=0;i<numNode;i++) { WingMesh.Nodes[i] = node_crds_temp[i]; }

    // count number of each element type
    int numC8=0; int numShell4=0; int numShell3=0; int numMass=0;
    for(i=0;i<numElem;i++)
    {
        if(elem_type[i]==24){ numC8++; }
        else if(elem_type[i]==20){ numShell4++; }
        else if(elem_type[i]==18){ numShell3++; }
        else if(elem_type[i]==3){ numMass++; }
    }

	// create arrays to store each element type
	C8M *WingElems; // use elements with incompatible modes
	WingElems = new C8M[numC8](); // array to store C8M elem objects
    //C8 *WingElems; // use elements without incompatible modes
	//WingElems = new C8[numC8](); // array to store C8 elem objects
    Shell3 *WingShells3;
    Shell4 *WingShells4;
    Mass1 *WingMass;
    if(numShell3 > 0)
    {
        control.shells = true;
        WingShells3 = new Shell3[numShell3](); // array to store Shell3 objects
    }
    if(numShell4 > 0)
    {
        control.shells = true;
        WingShells4 = new Shell4[numShell4](); // array to store Shell3 objects
    }
    if(numMass > 0)
    {
        WingMass = new Mass1[numMass](); // array to store Mass1 objects
    }

    int *dof_ind = new int [numNode+1];
    WingMesh.dof_ind = dof_ind;

    // Add elements to the mesh
    j=0; k=0; ind=0; ind2=0;
    for(i=0;i<numElem;i++)
    {
        if(elem_type[i]==24){ WingMesh.addElems(1, i, &WingElems[j++]); }
        else if(elem_type[i]==20){ WingMesh.addElems(1, i, &WingShells4[k++]); }
        else if(elem_type[i]==18){ WingMesh.addElems(1, i, &WingShells3[ind2++]); }
        if(elem_type[i]==3){ WingMesh.addElems(1, i, &WingMass[ind++]); }
    }

    // read element nodes & consolodate element numbers
    int *node_dof_temp = new int [numNode];
    for(i=0;i<numNode;i++){ node_dof_temp[i] = 0; } // assume 0 dof/node to start with
    for(i=0;i<numElem;i++)
	{
        if(elem_type[i]==24){ end = 8; } // 8 nodes
        else if(elem_type[i]==18){ end = 3; } // 3 nodes
        else if(elem_type[i]==20){ end = 4; } // 4 nodes
        else if(elem_type[i]==3){ end = 1; }  // 1 node
        for(j=0;j<end;j++) // for all elem nodes
        {
            ind = elem_def_temp[index2(i,j,8)];
            // need to select consolodated node number
            k = findPos2(ind, numNode, node_nums_temp);

            if(k==-1) {
                cout << "\nERROR! Could not find number for node " << j+1 << " Element " << i+1 << " aborting";
                return 0;
            }
            else { WingMesh.Elems[i]->addNode(j,k);
                node_dof_temp[k] = 3; // now increase to 3 tranlational dof
                if(end==3 || end==4){node_dof_temp[k] = 6;} // if shell element, increase dof to 6 (extra rotational dof)
                else if(end==1){node_dof_temp[k] = 0;} } // if mass element, decrease to zero (dof eliminated in constraint eq.)
        }
    }

    // create dof indicator array
    k = 0;
    for(i=0;i<numNode;i++)
    {
        dof_ind[i] = k;
        k += node_dof_temp[i];
    }
    dof_ind[numNode] = k; // total dof (order of mesh)
    int numDof  = k; // total dof in the mesh (full order of the sturtcural system)

	// delete temp data
	delete [] node_crds_temp;
    delete [] node_dof_temp;

    // set up array of pointer to (as yet) undefined nest and elsets
    mset *Nsets[MAX_SETS];
    mset *Elsets[MAX_SETS];

    // read in sets
    err = read_abq_sets(FE_file, &WingMesh, &numNset, &numElset, Nsets, Elsets, node_nums_temp, elem_nums_temp,MAX_SETS,MAX_SETS);
    if(err==1) {return 0;}

    // read in boundary and loading conditions from input file
    mset fixed(numDof); // set for fixed dof
    double *static_load = new double[numMission*numDof](); // static load vector
    err = read_abq_bc(FE_file, numDof, numNode, node_nums_temp, numNset, Nsets, &fixed, static_load, WingMesh.dof_ind);
    if(err==1) {return 0;}

    // read in material properties from input file
    mat *Materials[10]; // allow up to 10 materials
    err = read_abq_mat(FE_file, &numMat, Materials);
    if(err==1) {return 0;}

    // read and apply section data to elements
    err = read_abq_sec(FE_file, &WingMesh, numMat, Materials, numNset, numElset, Nsets, Elsets, node_nums_temp, elem_nums_temp);
    if(err==1) {return 0;}

    // if there are shell thk design variables, read them in now
    // can put back later

    delete [] node_nums_temp;
    delete [] elem_nums_temp;
	delete [] elem_def_temp;

    cout << "\n\nFE Mesh file read";
	cout << "\nNumber of Nodes = " << numNode;
    cout << "\nNumber of Brick Elements = " << numC8;
    cout << "\nNumber of Shell3 Elements = " << numShell3;
    cout << "\nNumber of Shell4 Elements = " << numShell4;
    cout << "\nNumber of Mass Elements = " << numMass;
	cout << "\nTotal Number of Elements = " << numElem;
    cout << "\nNumber of Node sets = " << numNset;
	cout << "\nNumber of Element sets = " << numElset;
    cout << "\nNumber of Materials = " << numMat;
    //cout << "\nNumber of Shell Thickness variables = " << numShell_var;

	// create element stiffness matrices
    //omp_set_num_threads(12); // stop multi-threading
    #pragma omp parallel for private(i)
	for(i=0;i<numElem;i++)
	{
		WingMesh.Elems[i]->makeKe(); // make the stiffness matrix
	}
    cout << "\nElement stiffness matrices made";

    if(control.comp_mass) // if mass matrix required
    {
        #pragma omp parallel for private(i)
        for(i=0;i<numElem;i++)
        {
            WingMesh.Elems[i]->makeMe(); // make the mass matrix
        }
        cout << "\nElement mass matrices made";
    }

    // -------------------------------------------------------------------- //
    //                                                                      //
    //      Step 2: Create the level set object                             //
    //              Initialize the implicit function                        //
    //                                                                      //
    // -------------------------------------------------------------------- //

	// MMA variables for SIMP
	// ITER = iteration
	// M = no. constraints
	// N = no. variables
    /*
	double GEPS = 1.0e-7; // convergence tolerance for constraints
	int *IYFREE; // WS - need to initialise to size N
	// XVAL = design variables
	double *XMMA; // optimal values from mmasub
	double *XMIN; // minimum values
	double *XMAX; // maximum values
	double *XOLD1; // Previous design variables
	double *XOLD2; // Previous previous design variables
	double *XLOW; // lower asymptote
	double *XUPP; // upper asymptote
	double *ALFA; // WS - need to initialise to size N
	double *BETA; // WS - need to initialise to size N
	double *A; // some coefficient - length M
	double *B; // WS - length M
	double *C; // some coefficient - length M
	double *Y; // Artificial value - length M
	double Z; // min-max variable Z
	double *ULAM; // dual variable - length M
	// F0VAL = objective value
	// FVAL = constraint values - length M
	double *FMAX; // rhs of constraints - length M
	// DF0DX = derivative of objective - length N
	// DFDX = derivative of constraints - length N*M : indexing = j*M + i
	double *P; // WS - length N*M
	double *Q; // WS - length N*M
	double *P0; // WS - length N
	double *Q0; // WS - length N
	double *UU; // WS - length M
	double *GRADF; // WS - length M
	double *DSRCH; // WS - length M
	double *HESSF; // WS - length M*N

	// Filtering variables
	double *filt_wgt;
	int *filt_ind;
	double filt_rad = 0.151; // Hard coded for now TODO
     */

    // compute element volumes (and total volume)
    double wing_mass_full = 0.0;
    double wing_vol_full = 0.0;
    double sum;
    for(i=0;i<numElem;i++)
	{
        WingMesh.Elems[i]->makeVol(); // first compute (& store) the element volume
        sum = WingMesh.Elems[i]->getVol();
        if(elem_type[i]==24){wing_vol_full += sum;} // only used for 3D level-set constraint
        if(elem_type[i]!=3){sum *= WingMesh.Elems[i]->getMat()->getDens();} // if not a mass element, x density
        wing_mass_full += sum; // add element mass
    }

    cout << "\n\nTotal wing mass = " << wing_mass_full;
    cout << "\nTotal wing full volume = " << wing_vol_full;

    // multiply volume constraint by wing volume
    for(i=0;i<lsprob.num;i++)
    {
        if(lsprob.con[i].type == 1)
        {
            lsprob.con[i].data[0] *= wing_vol_full;
        }
    }

    // create a level set object
    unLevelSet3D levelset(&WingMesh);

    if(control.method==1)
    {
		// find the zero lsf node set
		j = findSet("zero_lsf", numNset, Nsets);
		if(j==-1)     // if the set does not exist - abort
		{
			cout << "\nERROR in input, cannot find Nset for the zero lsf (should contain zero_lsf in the name) - abort";
			return 0;
		}

		levelset.initialize(Nsets[j], control.shells); // initialize the implicit lsf function (from outer faces of domain)
		// levelset.initialize(Nsets[j], 0); // SHELL EDIT

		if(restart != 0)
		{
			err = levelset.read_init(restart); // initialize from file (restart)
			if(err==1) {return 0;}
		}

		// either way, outer face nodes cannot become +ve
		mset *lsf_fix=0;

		// add nodes associated with mass elements to fixed set (as they are not included in level-set mesh)
		j = findSet("fix_lsf", numNset, Nsets);
		if(j==-1)     // if the set does not exist - abort
		{
			lsf_fix = new mset(1);
		}
		else
		{
			lsf_fix = Nsets[j];
		}

		// fixed lsf set
		levelset.setFixed(lsf_fix,0,0,0); // SHELL EDIT

		//if(control.shells){levelset.setFixed(lsf_fix,0,0,control.shells);}

		/*
		//k = (control.shells) ? 0 : lsf_fix->getPos(); //EDIT
		k = lsf_fix->pos; // EDIT
		k += numMass;
		lsf_fix->resize(k); // set storage space
		if(control.shells){lsf_fix->setPos(0);} // reset

		if(numMass>0)
		{
			int nrig[100];
			for(i=0;i<numMass;i++)
			{
				k = WingMass[i].getGlobal(0); // node number for mass elem
				lsf_fix->add(k); // add to set (referenced by level-set as fixed)

				// also check for rigid attachments
				ind =  WingMass[i].getRigid(nrig);
				if(ind>0){
					lsf_fix->resize(lsf_fix->getLen()+ind); // increase storage space
					for(j=0;j<ind;j++){ lsf_fix->add(nrig[j]);}
				}
			}
		}
		*/

		//levelset.setNeg(lsf_fix);
		// use the node numbering structure to define initial holes
		// Warning!! This will only work for a specific node number structure!
		/*int numH[3] = {62,6,3}; // number of holes in y,x,z
		int numNd[3] = {370,31,6}; // number of nodes in y,x,z
		int numHole = numH[0]*numH[1]*numH[2]; */
		int numHole;
		Coord3D *holeCentre = new Coord3D [MAX_HOLE];
		double *holeRadius = new double [MAX_HOLE];
		int *hole_type = new int [MAX_HOLE]();

		if(restart == 0)
		{
			//defineHoles(holeCentre, holeRadius, numH, numNd, &WingMesh, numHole);
			read_holes(FE_file, &numHole, holeCentre, holeRadius, hole_type); // read hole data from file

			for(i=0;i<numHole;i++) {
				if(hole_type[i]==0){levelset.addHoleSphere(holeCentre[i], holeRadius[i]);}
				else{levelset.addHoleCylinder(hole_type[i], holeCentre[i], holeRadius[i]);}
			}
			cout << "\n" << numHole << " initial holes added";

			// re-set initial narrow band
			levelset.NarBand();
		}

		delete [] holeRadius;
		delete [] holeCentre;
		delete [] hole_type;
    }
    // otherwise initialize densities
    /*
    else
    {
    	sum=1.0;
    	// find volume constraint value
    	for(i=0;i<lsprob.num;i++) { if(lsprob.con[i].type==1){ sum=lsprob.con[i].data[0]/wing_vol_full; break; } }

    	// initialise values
    	end = WingMesh.getNumEl();
    	for(i=0;i<end;i++){ levelset.setVR(i,sum); }

    	// Initialise MMA variables

    	k = lsprob.num;
    	IYFREE = new int [end]; // WS - need to initialise to size N
    	XMMA = new double [end]; // optimal values from mmasub
    	XMIN = new double [end]; // minimum values
    	XMAX = new double [end]; // maximum values
    	XOLD1 = new double [end]; // Previous design variables
    	XOLD2 = new double [end]; // Previous previous design variables
    	XLOW = new double [end]; // lower asymptote
    	XUPP = new double [end]; // upper asymptote
    	ALFA = new double [end]; // WS - length N
    	BETA = new double [end]; // WS - length N
    	A = new double [k]; // some coefficient - length M
    	B = new double [k]; // WS - length M
    	C = new double [k]; // some coefficient - length M
    	Y = new double [k]; // Artificial value - length M
    	ULAM = new double [k]; // dual variable - length M
    	// F0VAL = objective value
    	// FVAL = constraint values - length M
    	FMAX = new double [k](); // rhs of constraints - length M
    	// DF0DX = derivative of objective - length N
    	// DFDX = derivative of constraints - length N*M : indexing = j*M + i
    	P = new double [k*end]; // WS - length N*M
    	Q = new double [k*end]; // WS - length N*M
    	P0 = new double [end]; // WS - length N
    	Q0 = new double [end]; // WS - length N
    	UU = new double [k]; // WS - length M
    	GRADF = new double [k]; // WS - length M
    	DSRCH = new double [k]; // WS - length M
    	HESSF = new double [k*end]; // WS - length M*N

		#pragma omp parallel for private (i)
    	for(i=0;i<end;i++)
    	{
    		XMIN[i] = 0.001;
    		XMAX[i] = 1.0;
    		XLOW[i] = levelset.getVolrat(i);
    		XUPP[i] = levelset.getVolrat(i);
    		XOLD1[i] = levelset.getVolrat(i);
			XOLD2[i] = levelset.getVolrat(i);
    	}
    	for(i=0;i<k;i++)
    	{
    		A[i] = 0.0;
    		C[i] = 100.0 * Materials[0]->getMod() / wing_vol_full;
    		//FMAX[i] = lsprob.con[i].data[0]; // upper limit on constraint
    	}

    	// Initialise filtering variables

    	// compute element centres
    	end = WingMesh.getNumEl();
    	Coord3D cent;
    	Coord3D *el_cent = new Coord3D [end];
		#pragma omp parallel for private (i,j,k,cent)
    	for(i=0;i<end;i++)
    	{
    		// average node coords
    		cent.x=0.0; cent.y=0.0; cent.z=0.0; // reset
    		for(j=0;j<8;j++)
    		{
    			k = WingMesh.getElem(i)->getGlobal(j); // global node num
    			cent = cent + WingMesh.Nodes[k];
    		}
    		el_cent[i] = cent * 0.125;
    	}

    	// compute distances between centres
    	filt_wgt = new double [100*end]; // assume maximum 50 per element
    	filt_ind = new int [100*end];

		#pragma omp parallel for private (i,j,k,sum,ind,cent)
    	for(i=0;i<end;i++)
    	{
    		k=0; // count total in filter for each element
    		for(j=0;j<end;j++)
    		{
    			sum = filt_rad - veclen(el_cent[i] - el_cent[j]);
    			if(sum > -1.0e-6)
    			{
    				ind = index2(i,k,100);
    				filt_wgt[ind] = sum;
    				filt_ind[ind] = j;
    				k++;
    			}
    		}
    		if(k>99){cout << "\nError, more than 100 local in filter! abort";}
    		ind = index2(i,k,100);
    		filt_ind[ind] = -1; // end marker
    	}
    	delete [] el_cent;
    }
     */

    // -------------------------------------------------------------------- //
    //                                                                      //
    //     Step 5: Compute constant values, vectors and matrices for solve  //
    //                                                                      //
    // -------------------------------------------------------------------- //

    // variables required for aero-structural solve (primal & dual)
    double wing_weight, wing_vol;
    double *fg, *fgM, *f_prim;
    double *disp_prim, *disp_dual;
    double **ptr_dp, **ptr_dd, **ptr_eig;

    // for eigenvalue problems
    int eig_opt, nmodes;
    int *eig_order;
    double *eig_vals, *eig_vecs, *eig_vecs_old;

    // for buckling constraint
    int numBuck; // max number of buckling modes for each buckling zone
    int numBzone; // number of buckling zones (i.e. single panels, or upper skin, etc.)
    int num_buck_con; // total number of buckling KS functions
    mset **Bzones; // pointer to elsets that contain definition of buckling zones
    double *KS_buck_con; // buckling K-S constraint values

    // for stress constraint
    int numSzone; // number of stress zones (each has own K-S function)
    int num_stress_con; // total number of stress K-S constraints
    mset **Szones; // pointer to elsets defining the stress zones
    double *KS_stress_con; // stress K-S constraint values

    k=0;
    if(lsprob.obj == 1){ k=1; }
    else {
        for(i=0;i<lsprob.num;i++) { if(lsprob.con[i].type>4){ k=1; break; } } // buckling or stress constraint
    }
    // obj or constraint requires solution of 1 or more static problems
    if(k==1)
    {
        // displacement vector
        disp_prim = new double [numDof*numMission](); // one displacement vector for each mission

        // disp vector pointers
        ptr_dp = new double* [numMission]; // points to primary displacement vector for each mission
        ptr_dd = new double* [numMission]; // points to dual displacement vector for each mission

        // set-up data to feed into the sensitivity calculation
        for(i=0;i<numMission;i++)
        {
        	j = i*numDof;
        	ptr_dp[i] = &disp_prim[j];
        	ptr_dd[i] = ptr_dp[i]; // self-adjoint - EDIT
        	//if(lsprob.obj==1){ptr_dd[i] = (control.AE_feedback) ? &disp_dual[j] : ptr_dp[i];} // if no aeroelastic feedback, no adjoint
        }
    }

    // mass objective
    if(lsprob.obj==2)
    {
        for(i=0;i<lsprob.num;i++)
        {
            // buckling constraint
            if(lsprob.con[i].type==5)
            {
                numBuck = (int)lsprob.con[i].data[1]; // number of buckling modes per zone

                // find all elsets with "buckle" in name
                int *btemp = new int [numElset];

                ind=0;
                for(k=0;k<numElset;k++)
                {
                    j = findSet("buckle", 1, &Elsets[k]);
                    if(j>-1){ btemp[ind++] = k; } // record set position in Elsets
                }

                if(ind==0) {
                    cout << "\nERROR! Could not a buckling Elset when buckling is a constraint, set name must contain 'buckle' - aborting";
                    return 0;
                }

                numBzone = ind; // record number of zones to compute buckling
                Bzones = new mset* [ind];

                for(j=0;j<ind;j++) { Bzones[j] = Elsets[btemp[j]]; }
                delete [] btemp;

                // space for maximum number of active buckling constraints (using K-S agregation)
                num_buck_con = numBzone*numMission; // total number of buckling KS constraints
                KS_buck_con = new double [num_buck_con];
            }

            // stress constraint
            if(lsprob.con[i].type==6)
            {
                // find all elsets with "stress" in name
                int *btemp = new int [numElset];

                ind=0;
                for(k=0;k<numElset;k++)
                {
                    j = findSet("Stress", 1, &Elsets[k]);
                    if(j>-1){ btemp[ind++] = k; } // record set position in Elsets
                }

                if(ind==0) {
                    cout << "\nERROR! Could not a stress Elset when buckling is a constraint, set name must contain 'Stress' - aborting";
                    return 0;
                }

                numSzone = ind; // record number of zones to compute buckling
                Szones = new mset* [ind];

                for(j=0;j<ind;j++) { Szones[j] = Elsets[btemp[j]]; }
                delete [] btemp;

                // space for maximum number of active buckling constraints (using K-S agregation)
                num_stress_con = numSzone*numMission; // total number of buckling KS constraints
                KS_stress_con = new double [num_stress_con];
            }
        }
    }

    // eigenvalue objective
    if(lsprob.obj==3)
    {
        // for eigenvalue problems
        eig_opt = (int)lsprob.data[0];
        nmodes = eig_opt + 3; // in-case modes switch
        eig_opt--; // -1 for tracking in C

        ptr_eig = new double* [1];
        eig_vals = new double [nmodes];
        eig_order = new int [nmodes];
        eig_vecs  = new double [nmodes*numDof];
        eig_vecs_old  = new double [nmodes*numDof];
        for(i=0;i<nmodes;i++){eig_order[i]=i;} // order at start of optimization
    }

    // variables to pass data to velocity sub-solve function
    int numTot; // number of node + aux
    int numCon; // total number of constraints
    double *comp_sens = 0;
    double *mass_sens = 0;
    double *eig_max_sens = 0;
    double *vol_sens = 0;
    double *eig_ratio_sens = 0;
    double *buck_sens = 0;
    double *stress_sens = 0;
    double *con_violate;
    double *con_relax; // relaxtion factors
    int *con_active;   // active constraints
    double *con_pred;  // predicted constraint change
    double *con_old;   // store old constraint values (used to compute actual constraint changes)
    double *con_new;   // new constraint values

    // for shell thk design vars
    double **shp_sens_ptr; // pointers to shape sens arrays
    double **shell_sens_ptr;  // pointers to shell sens arrays
    double *shell_comp_sens = 0;
    double *shell_mass_sens = 0;
    double *shell_mass_sens_temp = 0;
    double *shell_eig_max_sens = 0;
    double *shell_vol_sens = 0;
    double *shell_eig_ratio_sens = 0;
    double *shell_buck_sens = 0;
    double *shell_stress_sens = 0;

    // pre-compute geometric shell thk sensitivities
    // can put back later - EDIT

    // compute maximum number of constraints
    k = lsprob.num;
    for(i=0;i<lsprob.num;i++)
    {
        // buckling constriants
        if(lsprob.con[i].type == 5)
        {
            k += num_buck_con - 1; // add extra constraints (using K-S agregation)
        }

        //  stress constriants
        if(lsprob.con[i].type == 6)
        {
            k += num_stress_con - 1; // add extra constraints (using K-S agregation)
        }
    }
    if(k > 0)
    {
        con_violate = new double [k];
        con_active = new int [k];
        con_relax = new double[k];
        con_pred = new double [k];
        con_old = new double [k];
        con_new = new double [k];
        for(i=0;i<k;i++){con_relax[i] = 1.0;} // initial relaxation factor
    }
    shp_sens_ptr = new double* [k+1]; // +1 for objective
    //if(numShell_var>0){shell_sens_ptr = new double * [k+1];} - EDIT

    // -------------------------------------------------------------------- //
    //                                                                      //
    //      Step 6: Start the level set optimization                        //
    //                                                                      //
    // -------------------------------------------------------------------- //

    int rewind_count = 0;
    int reint_cnt=30;
    double *obj = new double [control.maxItt];
    double *cnstr; if(lsprob.num > 0){ cnstr = new double [lsprob.num * control.maxItt]; }
    double c_max, c_min; // used to check convergence criterion

    int numEnt = 300*numC8 + 171*numShell3 + 300*numShell4; // maximum num entries in K
    int numEntM = numEnt; // maximum entries in M
    if(numMass > 0)
    {
        // add mass element entries
        for(i=0;i<numMass;i++)
        {
            j = WingMass[i].getNnodes()*3;
            c_max = (double)j;
            j = (c_max+1.0)*(c_max*0.5); // triangular storage space
            numEntM += j;
        }
    }

    matrix_spsy *K; // global stiffness matrix
    matrix_spsy *M; // global mass matrix
    hsl_matrix_type mtype = HSL_MATRIX_REAL_SYM_PSDEF;
    matrix_csc *Kcsc; // csc version of K
    matrix_csc *Mcsc; // csc version of M
    MA87d *Kcsc_fact;
    int *Kcsc_order = NULL;
    double elK_temp[300]; // array to store C8 element k matrix multiplied by vol ratio
    double *kp; // pointer for element stiffness matrix
    int *mloc;  // pointer for connected nodes of a mass element
    Mass1 *maptr; // mass element pointer
    elem *ep;  // generic element pointer
    Coord3D minus_z = {0.0,0.0,-1.0}; // used to set gravity vector

    double del_time;
    time_t ts, te; // timing variables

    double *NSW, *Nfact, *wgt_sens_fact;
    if(!control.SW_flag)
    { NSW = new double [numMission](); } // dummy weight for NSW problems
    else
    { Nfact = new double [numMission];
    	for(i=0;i<numMission;i++){Nfact[i]=1.0;} }
    wgt_sens_fact = new double [numMission]; // factor for weight sensitivity, for each mission

    // create an output file (for compliance, volume, etc. / iteration)
    string fname(infile);
    fname.append("_convergence.txt");
    char *cstr = new char [fname.length()+1];
    strlcpy(cstr, fname.c_str(), fname.length()+1);

    // open the output file
    ofstream outFile(cstr);
    if (!outFile.good()) {
		cout << "\nCould not open the results file: exit program...";
		return 0; // exit if file not opened
	}

    // write some column headings
    outFile << "itt\tWing mass\tVol frac";

    ind=0;
    if(lsprob.obj==1){ind=1;} // if objective is compliance
    else{  for(k=0;k<lsprob.num;k++){ if(lsprob.con[k].type>4){ind=1; break;} } } // buckling or stress constraint

    if(lsprob.obj==1)
    {
        for(i=0;i<numMission;i++)
        {
            outFile << "\tCompliance m" << i+1;
        }
        outFile << "\tTotal Compliance";
    }
    else if(control.comp_mass)
    {
        for(i=0;i<nmodes;i++)
        {
            outFile << "\tFreq" << i+1 << " (Hz)";
        }
    }
    for(i=0;i<lsprob.num;i++)
    {
        // Buckling constraint K-S values
        if(lsprob.con[i].type == 5)
        {
            for(j=0;j<num_buck_con;j++)
            {
                outFile << "\tBuckle K-S " << j+1;
            }
            break;
        }
    }
    for(i=0;i<lsprob.num;i++)
    {
        // Stress constraints K-S values
        if(lsprob.con[i].type == 6)
        {
            for(j=0;j<num_stress_con;j++)
            {
                outFile << "\tStress K-S " << j+1;
            }
            break;
        }
    }
    outFile << "\n";

    // open the timing output file
       ofstream timeFile("timing.txt");
       timeFile.precision(12);
       if (!timeFile.good()) {
   		cout << "\nCould not open the timing file: exit program...";
   		return 0; // exit if file not opened
   	}
    // write pre-proc time
    gettimeofday(&clck2,NULL);
    sum = (clck1.tv_sec+(clck1.tv_usec/1000000.0)) - (clck2.tv_sec+(clck2.tv_usec/1000000.0));
    timeFile << "Pre-proc took " << -sum;
    gettimeofday(&clck1,NULL);

    // write some column headings
    timeFile << "\nitt\tBoundary disc\tFE solve\tShape sens\tLevel-set opt";

    // the main optimization loop ...
    if(itt_in > 0){
        if(control.method == 1){levelset.NarBand();}
        for(i=0;i<itt_in;i++){obj[i]=0.0;}
    }
    itt=itt_in;
    if(itt==0){reint_cnt=0;} // re-int at first iteration

    do {
        ts = time(0);
        cout << "\n\n___Start iteration " << itt << "\n";
        outFile << itt;
        timeFile << "\n" << itt;
        // -------------------------------------------------------------------- //
        //                                                                      //
        //      Step 6.1:  extract volume ratios and boundary from lsf          //
        //                                                                      //
        // -------------------------------------------------------------------- //

        if(control.method == 1)
        {
        	levelset.getBound_para(); // get boundary discretization and volume ratios
        	if(reint_cnt == 0){levelset.fast_march2(0); reint_cnt=30; levelset.getBound_para();} // reinitialize
        }
        //levelset.outBound(infile, itt, 0, 0);
        //levelset.setVR(151215, 4.00E-06); // FD check

        cout << "\nStructure boundary and volume ratios computed\n";
        if(control.pinfo>0) {
            levelset.outStruct(infile, itt);// output the element volume ratios & lsf
            //if(numShell_var>0){ sprintf(foname,"%s_shell-thk_it%i",infile,itt);
            //   outShellMesh(&WingMesh, numShell3, numShell4, 0, 0, 0, 0, foname); } // shell mesh & thk output
        }

        gettimeofday(&clck2,NULL);
        sum = (clck1.tv_sec+(clck1.tv_usec/1000000.0)) - (clck2.tv_sec+(clck2.tv_usec/1000000.0));
        timeFile << "\t" << -sum;
        gettimeofday(&clck1,NULL);

        // -------------------------------------------------------------------- //
        //                                                                      //
        //      Step 6.2:  assemble the global stiffness / mass matrices        //
        //                 compute current wing weight                          //
        //                                                                      //
        // -------------------------------------------------------------------- //

        // create global stiffness (maybe mass) matrix
        K = new matrix_spsy(numDof,numEnt);
        if(control.comp_mass){ M = new matrix_spsy(numDof,numEntM); }

        ind=0; ind2=0;
        wing_weight=0.0; wing_vol=0.0;
        for(i=0;i<numElem;i++)
        {
            ep = WingMesh.Elems[i]; // pinter to current element

            if(WingMesh.el_type[i]==24)
            {
                // need to multiply element stiffness matrix by volume ratio
                sum = levelset.getVolrat(i); // element vol ratio

                // only multiply if vol ratio < 1
                if(sum > 0.0)
                {
                	if(control.method==2){sum = pow(sum,3.0);} // penalization for SIMP
                    if(sum < 0.999999) {
                        for(j=0;j<300;j++) {
                            elK_temp[j] = sum*ep->getKe()[j];
                        }
                        kp = elK_temp; // point to multiplied matrix
                    }
                    // otherwise use original matrix
                    else {kp = ep->getKe();}

                    // assemble element into Global matrices - sparse storage
                    addElem(8, kp, K, 3, WingMesh.dof_ind, ep->getNodes(), &ind);

                    // mass matrix (used for frequency problems)
                    if(control.comp_mass)
                    {
                        if(sum < 1.01E-6){sum = 1.0E-8;} // extra penalization on void mass region
                        if(sum < 0.999999) {
                            for(j=0;j<300;j++) {
                                elK_temp[j] = sum*ep->getMe()[j];
                            }
                            kp = elK_temp; // point to multiplied matrix
                        }
                        // otherwise use original matrix
                        else {kp = ep->getMe();}

                        // assemble element into Global matrices - sparse storage
                        addElem(8, kp, M, 3, WingMesh.dof_ind, ep->getNodes(), &ind2);
                    }

                    // also compute total wing weight
                    sum = levelset.getVolrat(i); // element vol ratio
                    sum *= ep->getVol(); // volume x volume ratio
                    wing_weight += sum * ep->getMat()->getDens(); // weight = g x volume x density
                    wing_vol += sum; // sum volume
                }
            }

            // assemble shell element into Global matrix
            else if(WingMesh.el_type[i]==18 || WingMesh.el_type[i]==20)
            {
                // can add thickness change here later!!
                if(WingMesh.el_type[i]==18){addElem(3, ep->getKe(), K, 6, WingMesh.dof_ind, ep->getNodes(), &ind);}
                else{addElem(4, ep->getKe(), K, 6, WingMesh.dof_ind, ep->getNodes(), &ind);}
                if(control.comp_mass) {
                    if(WingMesh.el_type[i]==18){addElem(3, ep->getMe(), M, 6, WingMesh.dof_ind, ep->getNodes(), &ind2);}
                    else{addElem(4, ep->getMe(), M, 6, WingMesh.dof_ind, ep->getNodes(), &ind2);}
                }
                wing_weight += ep->getVol() * ep->getMat()->getDens(); // weight = g x volume x density
                // do not update wing volume (only used in level-set constraint)
            }

            // need to add lumped masses to total mass (and mass matrix if required)
            else if(WingMesh.el_type[i]==3)
            {
                if(control.comp_mass)
                {
                    maptr = (Mass1*)ep; // convert to mass elem pointer
                    j = maptr->getNnodes(); // number of connected nodes
                    mloc = maptr->getCNodes(); // connected node numbers
                    addElem(j, maptr->getMe(), M, 3, WingMesh.dof_ind, mloc, &ind2);
                }
                wing_weight += ep->getVol(); // mass of mass element
            }
        }
        K->setLen(ind); // update number of entries

        // remove fixed dof
        K->removeDof(fixed.getPos(),fixed.nums,1);

        // convert matrix storage
        Kcsc = new matrix_csc(mtype, K);
        for(i=0;i<lsprob.num;i++){if(lsprob.con[i].type==5){break;}}
        if(i==lsprob.num){delete K;} // do not delete if buckling analysis required

        cout << "\nStiffness Matrix made\n";

        //Kcsc->mOut(20);

        // print K to file
        /*sprintf(foname,"%s_Kmatrix_it%i",infile,itt);
        Kcsc->mOut_file(1e20,foname);

        /*{
        	ifstream mainFile;
        	mainFile.open("jCRM_06d_Buckling_Kmatrix_it0.txt");

        	int ent = 23232804;
        	matrix_spsy Ka(597954,ent,1);

        	int *irn = Ka.get_irn();
        	int *jcn = Ka.get_jcn();
        	double *Aa = Ka.getA();

        	for(i=0;i<ent;i++)
        	{
        		mainFile >> irn[i];
        		mainFile >> jcn[i];
        		mainFile >> Aa[i];
        	}

        	 // convert matrix storage
        	 matrix_csc *KAcsc = new matrix_csc(mtype, &Ka);

        	 // factor KAcsc
        	         MA87d *KAcsc_fact = new MA87d(KAcsc->getOrd(),0,1);
        	         KAcsc_fact->call_analyse(KAcsc->get_ptr(), KAcsc->get_row(), &Kcsc_order);
        	         KAcsc_fact->call_factor(KAcsc->get_ptr(), KAcsc->get_row(), KAcsc->getA(), Kcsc_order);

        	double eval[10];
        	double *evec = new double [10*numDof];
        	i = eig_solve_simple3(5, KAcsc_fact, Kcsc_order, &fixed, eval, evec, 2);
        	delete [] evec;
        }*/

        // repeat for mass matrix
        if(control.comp_mass)
        {
            M->setLen(ind2);
            M->removeDof(fixed.getPos(),fixed.nums,1);
            Mcsc = new matrix_csc(mtype, M);
            delete M;
            cout << "\nMass Matrix made\n";
        }

        // record volume / mass as obj or constraint values
        cout << "\nWing mass = " << wing_weight;
        sum = wing_vol / wing_vol_full; // vol fraction
        cout << "\nVolume fraction = " << sum;
        outFile << "\t" << wing_weight << "\t" << sum;

        if(lsprob.obj==2){ obj[itt] = wing_weight; } // if mass is objective or constraint
        else
        {
            for(i=0;i<lsprob.num;i++)
            {
                if(lsprob.con[i].type==2){
                    j = lsprob.num * itt + i; // point to place in cnstr
                    cnstr[j] = wing_weight; // store mass as constraint value
                    break;
                }
                else if(lsprob.con[i].type==1){
                    j = lsprob.num * itt + i; // point to place in cnstr
                    cnstr[j] = wing_vol; // store volume as constraint value
                    break;
                }
            }
        }

        wing_weight *= 9.81; // multiply by grav const.

        // -------------------------------------------------------------------- //
        //                                                                      //
        //      Step 6.3:  compute the self-weight load vector                  //
        //                                                                      //
        // -------------------------------------------------------------------- //

        ind=0;
        if(lsprob.obj==1){ind=1;} // if objective is compliance
        else{  for(k=0;k<lsprob.num;k++){ if(lsprob.con[k].type>4){ind=1; break;} } } // buckling or stress constraint
        if(ind==1)
        {
            // need to compute the self weight load for N=1
            // assume self weight acts in -ve Z direction (for steady level flight)
            fg = new double [numDof*numMission]; // self-weight load for 1g (-z direction)
            //fgM = new double [numDof*numMission]; //  self-weight loading vector (for each mission)
            //f_prim = new double [numDof*numMission]; // converged total force vector from primary solution (for each mission)
            //cblas_dcopy(numDof, static_load, 1, fg, 1); // initialize as the static load (assume this is from other self-weight loads)

            // now add the self weight loads from the current wing-box structure
            if(control.SW_flag)
            {
                ind=0; ind2=0; k=0; j=0;
                for(i=0;i<numElem;i++)
                {
                    if(WingMesh.el_type[i] == 24)
                    {
                        sum = levelset.getVolrat(i);
                        //if(sum > 1.0E-3) // only include elements that are above small volume ratio
                        if(sum > 0.0)
                        {
                            sum *= 9.81; // multiply vol ratio by gravity constant
                            WingElems[ind++].addWeight(minus_z, sum, WingMesh.dof_ind, fg);
                        }
                    }
                    // shell element
                    else if(WingMesh.el_type[i] == 18)
                    {
                        WingShells3[ind2++].addWeight(minus_z, 9.81, WingMesh.dof_ind, fg);
                    }
                    else if(WingMesh.el_type[i] == 20)
                    {
                        WingShells4[j++].addWeight(minus_z, 9.81, WingMesh.dof_ind, fg);
                    }
                    // mass element
                    else if(WingMesh.el_type[i] == 3)
                    {
                        WingMass[k++].addWeight(minus_z, 9.81, WingMesh.dof_ind, fg);
                    }
                }
                cout << "\nSelf-weight load vector computed...";
            }

            /*
            // self weight load vector for each mission (fg x load factor)
			for(i=0;i<numMission;i++)
			{
				sum = (control.SW_flag) ? Missions.Nfact[i] : 0.0; // load factor
				for(j=0;j<numDof;j++)
				{
					fgM[index2(i,j,numDof)] = fg[j] * sum; // self-weight x load factor
				}
			}*/
        }

        // -------------------------------------------------------------------- //
        //                                                                      //
        //      Step 6.4:  solve the governing equations                        //
        //                  aero-structural system & adjoint                    //
        //                  or an eigenvalue problem                            //
        //                  or a flutter problem                                //
        //                                                                      //
        // -------------------------------------------------------------------- //

        // factor Kcsc
        Kcsc_fact = new MA87d(Kcsc->getOrd(),0,1);
        Kcsc_fact->call_analyse(Kcsc->get_ptr(), Kcsc->get_row(), &Kcsc_order);
        Kcsc_fact->call_factor(Kcsc->get_ptr(), Kcsc->get_row(), Kcsc->getA(), Kcsc_order);

        ind=0;
        if(lsprob.obj==1){ind=1;} // if objective is compliance
        else{  for(k=0;k<lsprob.num;k++){ if(lsprob.con[k].type>4){ind=1; break;} } } // buckling or stress constraint
        if(ind==1)
        {
            // removed fixed dof from input force vector (i.e. currently static_load)
            ind = 1;
            k = fixed.nums[ind];
            j=0;
            end = dof_ind[numNode]; // total dof
            for(i=0;i<end;i++)
            {
                if(i != k) {		// If dof not fixed copy across the value
                	for(ind2=0;ind2<numMission;ind2++) // for each load case
                	{
                		fg[(ind2*Kcsc->getOrd())+j] = static_load[(ind2*numDof)+i];
                	}
            		j++;
                    //if(fabs(static_load[i])>0.0){cout << "\n" << i << ": "<< static_load[i]; }
                }
                else {				// Else start looking at next fixed dof
                    k = fixed.nums[++ind];
                }
            }

            // solve primary static problem
            Kcsc_fact->call_solve(numMission,fg,Kcsc_order);
        	disp_dual = disp_prim; // self-adjoint

        	 // Re-insert zero displacements into the displacement array
			ind = 1;
			j=0;
			k = fixed.nums[ind];

			for(i=0;i<end;i++)
			{
				if(i != k) {		// If dof not fixed copy accross the value
					for(ind2=0;ind2<numMission;ind2++) // for each load case
					{
						disp_prim[(ind2*numDof)+i] = fg[(ind2*Kcsc->getOrd())+j];
					}
					j++;
				}
				else {				// Else dof is fixed! i.e zero displacement -> then look for next fixed dof
					for(ind2=0;ind2<numMission;ind2++) // for each load case
					{
						disp_prim[(ind2*numDof)+i] = 0.0;
					}
					k = fixed.nums[++ind];
				}
			}

            // output nodal displacment
            if(control.pinfo==2) {
                sprintf(foname,"%s_disp_it%i",infile,itt);
                outNodalVTK(&WingMesh, numC8, numShell3, numShell4, numMission, disp_prim, 0, foname);
                if(lsprob.obj==1){ sprintf(foname,"%s_adjoint_it%i",infile,itt);
                    outNodalVTK(&WingMesh, numC8, numShell3, numShell4, numMission, disp_dual, 0, foname); }
            }

            // compute the compliance (f_prim * disp_prim)
            // sum for each load case
            if(lsprob.obj==1)
            {
                obj[itt] = 0.0;
                for(i=0;i<numMission;i++)
                {
                    sum = cblas_ddot(numDof, &static_load[i*numDof], 1, &disp_prim[i*numDof], 1); // EDIT
                    obj[itt] += sum; //*wgt_case[i]; - EDIT
                    cout << "\nCompliance for case " << i+1 << " = " << sum;
                    outFile << "\t" << sum;
                }
                cout << "\nTotal compliance = " << obj[itt];
                outFile << "\t" << obj[itt];
            }

            // clear load vector memory
            delete [] fg;
            //delete [] fgM;
            //delete [] f_prim;
        }
        if(control.comp_mass)
        {
            // funciton to use ARPACK reverse communication to solve generalized eigenvalue problem (mode 3)
            err = eig_solve(nmodes, Kcsc_fact, Kcsc_order, Mcsc, &fixed, eig_vals, eig_vecs, 2);

            // mode tracking
            if(itt > itt_in)
            {
                // correlate modes from last iteration
                mode_corr_real(nmodes, nmodes, numDof, numDof, eig_vecs_old, eig_vecs, eig_order);
            }

            // save new eigenvectors for correlation next iteration (with tracking)
            for(i=0;i<nmodes;i++)
            {
                cblas_dcopy(numDof, &eig_vecs[eig_order[i]*numDof], 1, &eig_vecs_old[i*numDof], 1);
                outFile << "\t" << sqrt(eig_vals[eig_order[i]])*0.1591549430919;
            }

            if(lsprob.obj==3){ obj[itt] = sqrt(eig_vals[eig_order[eig_opt]]); } // opt eigenvalue (tracked)

            // output eigenvector file
            if(control.pinfo==2) {
                sprintf(foname,"%s_EigVecs_it%i",infile,itt);
                j = (nmodes < 7) ? nmodes : 6;
                outNodalVTK(&WingMesh, numC8, numShell3, numShell4, j, eig_vecs, 0, foname);
            }
        }

        // buckling constraint - perform series of linear buckling analyses + adjoint + sensitivity calculation
        for(k=0;k<lsprob.num;k++)
        {
            if(lsprob.con[k].type==5)
            {
                buck_sens = new double [num_buck_con*levelset.getTot()](); // storage for buckling sensitivities
                //if(numShell_var>0){ shell_buck_sens = new double [num_buck_con*numShell_var](); } // shell variable senstivities

                // call buckling optimization sub-function
                /*Buckle_opt(numDof, numBuck, numBzone, Bzones, numMission, &Missions, &WingMesh, Transfer,
                           numC8, numShell3, numShell4, WingShells3, WingShells4, ptr_dp, &fixed, Kcsc, KS_buck_con, &levelset,
                           Kcsc_fact, Kcsc_order, Qd, fa, La, buck_sens, shell_buck_sens, numShell_var,
                           shell_var, infile, K, itt, control.pinfo);
                */ // - EDIT
                delete K; // can delete now

                // save max constraint value
                sum = -1.0e20;
                for(j=0;j<num_buck_con;j++)
                {
                    sum = (KS_buck_con[j] > sum) ? KS_buck_con[j] : sum;
                    outFile << "\t" << KS_buck_con[j];
                }

                j = lsprob.num * itt + k; // point to place in cnstr
                cnstr[j] = sum;
                break;
            }
        }

        // stress constraint - compute K-S + solve adjoint
        for(k=0;k<lsprob.num;k++)
        {
            if(lsprob.con[k].type==6)
            {
                stress_sens = new double [num_stress_con*levelset.getTot()](); // storage for buckling sensitivities
                /*if(numShell_var>0){ shell_stress_sens = new double [num_stress_con*numShell_var](); } // shell variable senstivities

                // call stress optimization sub-function
                Stress_Con_func(numDof, numSzone, Szones, numMission, &Missions, &WingMesh, Transfer, numC8, numShell3, numShell4, ptr_dp, &fixed, KS_stress_con,
                                &levelset, Kcsc_fact, Kcsc_order, Qd, fa, La, stress_sens, shell_stress_sens,numShell_var, shell_var, infile);
                */ // -EDIT
                // save max constraint value
                sum = -1.0e20;
                for(j=0;j<num_stress_con;j++)
                {
                    sum = (KS_stress_con[j] > sum) ? KS_stress_con[j] : sum;
                    outFile << "\t" << KS_stress_con[j];
                }

                j = lsprob.num * itt + k; // point to place in cnstr
                cnstr[j] = sum;
                break;
            }
        }

        // force update of convergence file
        outFile << "\n"; // new line for next iteration
        outFile.close();
        outFile.open(cstr, ios::out | ios::app);

        // delete [] Kcsc_order; Kcsc_order = NULL; // reset ordering
        delete Kcsc_fact; // delete factorization
        delete Kcsc; // delete the stiffness matrix
        if(control.comp_mass){ delete Mcsc; }

        // -------------------------------------------------------------------- //
        //                                                                      //
        //      Step 6.6:  check for convergence                                //
        //                                                                      //
        // -------------------------------------------------------------------- //

        // Check convergence criteria
        // calculate difference between current and target constraint
        ind = 1; // check for violated constraints
        j = lsprob.num;
        for(i=0;i<j;i++)
        {
            k = j*itt + i; // pointer for cnstr
            sum = cnstr[k] - lsprob.con[i].data[0];
            if(sum > 0.0)
            {
                cout << "\nConstraint " << i+1 << " violated by " << sum;
                ind = 0;
            }
        }

        if(itt > (8+itt_in) && ind == 1)
        {
            c_max=0.0; c_min=1.0E+20;
            for(i=itt;i>itt-10;i--) // work out maximum objective change over last 10 itteration
			{
				sum = obj[i];
				c_max = (sum > c_max) ? sum : c_max; //update maximum objective value
				c_min = (sum < c_min) ? sum : c_min; //update minimum objective value
			}
			printf("\nMax obj = %lf, Min obj = %lf\n",c_max,c_min);
			sum = (c_max - c_min) / (c_max + c_min); // normalized max change
            // if max objective chnage over last 10 itteration is < gamma then stop!!
			if(sum < control.gamma) {
				printf("Converged!!! Criteria = %lf",sum);
				break;
			}
			else {
				printf("Convergence criteria = %lf",sum);
			}
        }

        gettimeofday(&clck2,NULL);
        sum = (clck1.tv_sec+(clck1.tv_usec/1000000.0)) - (clck2.tv_sec+(clck2.tv_usec/1000000.0));
        timeFile << "\t" << -sum;
        gettimeofday(&clck1,NULL);

        // -------------------------------------------------------------------- //
        //                                                                      //
        //      Step 6.7:  compute sensitivities at the boundary                //
        //                  & shell thk sensitivities if required               //
        //                                                                      //
        // -------------------------------------------------------------------- //

        numTot = levelset.getTot();
        if(lsprob.obj==1)
        {
            // first compute the weight sensitivity factors
            for(i=0;i<numMission;i++)
            {
                //if(Missions.getDyPres(i) < 1.0e-3 || control.fix_aoa)
                {
                    wgt_sens_fact[i] = 0.0; // zero if "wind off" or fixed aoa
                }
            }

            //if(control.SW_flag){ levelset.shapeSens(1, numMission, (void**)ptr_dp, (void**)ptr_dd, 0, wgt_sens_fact, wgt_case, Missions.Nfact, comp_sens,1); }
            //else{ levelset.shapeSens(1, 1, (void**)ptr_dp, (void**)ptr_dd, 0, wgt_sens_fact, wgt_case, NSW, comp_sens, 1); }
            if(control.method==1)
            {
                comp_sens = new double [numTot]();
            	if(control.SW_flag){ levelset.shapeSens(1, numMission, (void**)ptr_dp, (void**)ptr_dd, 0, wgt_sens_fact, wgt_case, Nfact, comp_sens,1); }
            	else{ levelset.shapeSens(1, numMission, (void**)ptr_dp, (void**)ptr_dd, 0, wgt_sens_fact, wgt_case, NSW, comp_sens, 1); }
            }
            // Element density derivatives
 //           else
//            {
//            	end = WingMesh.getNumEl(); // number of elements
//            	comp_sens = new double [end](); // used for raw compliance sensitivities
//
//            	if(control.SW_flag){
//					#pragma omp parallel for private (i)
//					for(i=0;i<end;i++) {
//						WingElems[i].elemSens(1, numMission, ptr_dp, ptr_dd, wgt_sens_fact, wgt_case, comp_sens, i, Nfact, 1);
//					}
//            	}
//            	else{
//					#pragma omp parallel for private (i)
//					for(i=0;i<end;i++) {
//						WingElems[i].elemSens(1, numMission, ptr_dp, ptr_dd, wgt_sens_fact, wgt_case, comp_sens, i, NSW, 1);
//					}
//            	}
//
//            	// processing of raw sensitivities
//				#pragma omp parallel for private (i,sum)
//				for(i=0;i<end;i++) {
//					sum = levelset.getVolrat(i);
//					sum *= 3.0*sum; // -3*p^2
//					comp_sens[i] *= sum;
//				}
//
//				// Filtering
//				double sum2;
//				double *sens_temp = new double [end];
//				#pragma omp parallel for private (i,j,ind,sum,sum2)
//				for(i=0;i<end;i++) {
//					sum = 0.0; sum2 = 0.0;
//					// for all elements in filter radius
//					for(j=0;j<100;j++) {
//						ind = index2(i,j,100); // indicator for element number in the filter radius
//						if(filt_ind[ind]!=-1) {
//							sum += comp_sens[filt_ind[ind]] * filt_wgt[ind] * levelset.getVolrat(filt_ind[ind]);
//							sum2 += filt_wgt[ind];
//						}
//						else {
//							break;
//						}
//					}
//					// filtered sensitivity
//					sens_temp[i] = sum / (sum2 * levelset.getVolrat(i));
//				}
//				// copy temp to main
//				cblas_dcopy(end,sens_temp,1,comp_sens,1);
//				delete [] sens_temp;
//			}

            // shell sensitivities
            // put back later - EDIT
        }

        else if(lsprob.obj==3)
        {
            j=eig_order[eig_opt];
            ptr_eig[0] = &eig_vecs[j*numDof];
            if(control.method==1)
            {
            	eig_max_sens = new double [numTot]();
            	levelset.shapeSens(1, 1, (void**)ptr_eig, (void**)ptr_eig, 0, &eig_vals[j], wgt_case, NSW, eig_max_sens, 2);
            }
            else
            {
				// TODO
            }

            // shell sensitivities
            // put back later - EDIT
        }

        // -------------------------------------------------------------------- //
        //                                                                      //
        //      Step 6.8:  Set up pointers etc. for sensitivity data            //
        //                                                                      //
        // -------------------------------------------------------------------- //

        // need:
        // total number of constraints
        // pointer array to shape sensitivities for objective and each constraint
        // pointer arrays to shell thk sensitivities for objective and each constraint
        // current constraint violations for each constraint

        numCon = 0;
        // first pointer in sens is objective
        if(lsprob.obj == 1)
        {
            shp_sens_ptr[0] = comp_sens;
            //if(numShell_var>0){ shell_sens_ptr[0] = shell_comp_sens; } - EDIT
        }
        else if(lsprob.obj == 2)
        {
            // make mass_sens
            mass_sens = new double [numTot];
            //if(numShell_var>0){shell_mass_sens_temp = new double [numShell_var];
            //    cblas_dcopy(numShell_var, shell_mass_sens, 1, shell_mass_sens_temp, 1);} // EDIT copy (as it gets modified in SLP sub-solve)
            sum = -1.0 * WingElems[0].getMat()->getDens(); // density of level set solid elems
            for(j=0;j<numTot;j++){ mass_sens[j] = sum; } // mass sens (assume all elements have same density)
            shp_sens_ptr[0] = mass_sens;
            //if(numShell_var>0){ shell_sens_ptr[0] = shell_mass_sens_temp; } // pre-computed as sum of element area x density for each design var
        }
        else if(lsprob.obj == 3)
        {
            shp_sens_ptr[0] = eig_max_sens;
            //if(numShell_var>0){ shell_sens_ptr[0] = shell_eig_max_sens; } - EDIT
        }

        // then look for constraints
        for(i=0;i<lsprob.num;i++)
        {
            // Volume of level-set region
            if(lsprob.con[i].type == 1)
            {
                con_new[numCon] = wing_vol;
                con_violate[numCon++] = lsprob.con[i].data[0] - cnstr[lsprob.num*itt + i]; // constraint violation
                // make vol_sens
                if(control.method==1)
                {
                	vol_sens = new double [numTot];
                	for(j=0;j<numTot;j++){ vol_sens[j] = -1.0; } // volume sens
                }
                else
                {
                	end = WingMesh.getNumEl();
                	vol_sens = new double [end];
                	for(j=0;j<end;j++){ vol_sens[j] = WingMesh.getElem(j)->getVol();
                		//cout << "\n" <<  vol_sens[j];
                	} // volume sens
                }
                shp_sens_ptr[numCon] = vol_sens;
               // if(numShell_var>0){ shell_sens_ptr[numCon] = shell_vol_sens; } // EDIT pre-computed as zeros (does not effect level-set volume)
            }
            // Total mass of the wing
            else if(lsprob.con[i].type == 2)
            {
                con_new[numCon] = wing_weight/9.81;
                con_violate[numCon++] = lsprob.con[i].data[0] - cnstr[lsprob.num*itt + i]; // constraint violation
                // make mass_sens
                mass_sens = new double [numTot];
                //if(numShell_var>0){shell_mass_sens_temp = new double [numShell_var];
                //    cblas_dcopy(numShell_var, shell_mass_sens, 1, shell_mass_sens_temp, 1);} // copy (as it gets modified in SLP sub-solve)
                sum = -1.0 * WingElems[0].getMat()->getDens();
                for(j=0;j<numTot;j++){ mass_sens[j] = sum; } // mass sens (assume all elements have same density)
                shp_sens_ptr[numCon] = mass_sens;
                //if(numShell_var>0){ shell_sens_ptr[numCon] = shell_mass_sens_temp; } // pre-computed as sum of element area x density for each design var
            }
            // Eigenvalue ratio
            /*else if(lsprob.con[i].type == 3)
            {
                con_violate[numCon++] = lsprob.con[i].data[0] - cnstr[lsprob.num*itt + i]; // constraint violation
                shp_sens_ptr[numCon] = eig_ratio_sens;
                if(numShell_var>0){ shell_sens_ptr[numCon] = shell_eig_ratio_sens; }
            }*/

            // Buckling constraints
            else if(lsprob.con[i].type == 5)
            {
                for(j=0;j<num_buck_con;j++)
                {
                    con_new[numCon] = KS_buck_con[j];
                    con_violate[numCon++] = -KS_buck_con[j]; // constraint violation
                    shp_sens_ptr[numCon] = &buck_sens[j*numTot];
                    //if(numShell_var>0){ shell_sens_ptr[numCon] = &shell_buck_sens[j*numShell_var]; } - EDIT
                }
            }
            // Stress constraints
            else if(lsprob.con[i].type == 6)
            {
                for(j=0;j<num_stress_con;j++)
                {
                    con_new[numCon] = KS_stress_con[j];
                    con_violate[numCon++] = -KS_stress_con[j]; // constraint violation
                    shp_sens_ptr[numCon] = &stress_sens[j*numTot];
                    //if(numShell_var>0){ shell_sens_ptr[numCon] = &shell_stress_sens[j*numShell_var]; } - EDIT
                }
            }
        }

        // if not first iteration
        err = 0; // indicate no rewind
        if(itt > itt_in && control.method==1)
        {
            // for all constraints
            k = 0; // count constraints
            for(i=0;i<numCon;i++)
            {
                if(con_active[i] != 0)
                {
                    sum = con_pred[i]/(con_new[i]-con_old[i]); // predicted / actual

                    // compute constraint relaxation factor
                    if(sum > 1.0){con_relax[i] *= 2.0;} // increase
                    if(sum < 0.5){con_relax[i] *= 0.5;} // decrease
                    if(con_relax[i] > 1.0){con_relax[i] = 1.0;} // upper limit
                    else if(con_relax[i] < 0.125){con_relax[i] = 0.125;} // lower limit

                    // rewind feature
                    if(sum<0.0){err=1;} // rewind!!

                    printf("\nRelaxation factor for %i = %f",i+1,con_relax[i]);
                    // divide shape sens by relaxation factor
                    // sum = 1.0/con_relax[i];
                    /*sum = con_relax[i];
                    k = i+1; // point to constraint shape sens
                    for(j=0;j<numTot;j++) { shp_sens_ptr[k][j] *= sum; }*/
                    con_violate[i] *= con_relax[i];
                }
            }
        }

        gettimeofday(&clck2,NULL);
        sum = (clck1.tv_sec+(clck1.tv_usec/1000000.0)) - (clck2.tv_sec+(clck2.tv_usec/1000000.0));
        timeFile << "\t" << -sum;
        gettimeofday(&clck1,NULL);

        // copy new -> old
        if(err==0 || rewind_count==0)
        {
        for(i=0;i<numCon;i++) { con_old[i] = con_new[i]; }

        // -------------------------------------------------------------------- //
        //                                                                      //
        //      Step 6.9:  Compute velocities, shell thk change, gradients      //
        //                 update the lsf & shell thk vars                      //
        //                                                                      //
        // -------------------------------------------------------------------- //

        // 1st set max and min step change in thickness variables for this iterartion
        /*#pragma omp parallel for private(i,sum,c_max,c_min)
        for(i=0;i<numShell_var;i++)
        {
            sum = 0.01*(shell_max[i] - shell_min[i]); // hard move limit
            // check move limit with real side limits
            c_max = shell_thk[i] + sum;
			c_min = shell_thk[i] - sum;
			shell_max_it[i] = (c_max < shell_max[i]) ? sum : shell_max[i] - shell_thk[i];
			shell_min_it[i] = (c_min > shell_min[i]) ? -sum : shell_min[i] - shell_thk[i];
        }*/

        if(control.method==1)
        {
			// Obtain the velocity function & step updates for shell thk vars
			//levelset.getVel(numCon+1, con_violate, shp_sens_ptr, numShell_var, shell_sens_ptr, shell_min_it, shell_max_it, shell_step, con_active, con_pred, control.pinfo);
			levelset.getVel(numCon+1, con_violate, shp_sens_ptr, 0, 0, 0, 0, 0, con_active, con_pred, control.pinfo);

			gettimeofday(&clck2,NULL);
			        sum = (clck1.tv_sec+(clck1.tv_usec/1000000.0)) - (clck2.tv_sec+(clck2.tv_usec/1000000.0));
			        timeFile << "\t" << -sum;
			        gettimeofday(&clck1,NULL);

			if(control.pinfo==2){ levelset.outBound(infile, itt, numCon+1, shp_sens_ptr); } // output the boundary surface mesh
        }
        // use MMA optimizer to update element densities
        /*
        else
        {
        	// call mma subroutine
        	for(i=0;i<lsprob.num;i++){con_violate[i]*=-1.0;}
        	k=itt+1;
        	end = WingMesh.getNumEl();
        	mmasub_(&k, &(lsprob.num), &end, &GEPS, IYFREE, levelset.getVolratio(), XMMA,
        			XMIN, XMAX, XOLD1, XOLD2, XLOW, XUPP,
        			 ALFA, BETA, A, B, C, Y, &Z, ULAM,
        			 &obj[itt], con_violate, FMAX, comp_sens, vol_sens,
        			 P, Q, P0, Q0, UU, GRADF, DSRCH, HESSF);

        	// copy data
        	cblas_dcopy(end,XOLD1,1,XOLD2,1);
        	cblas_dcopy(end,levelset.getVolratio(),1,XOLD1,1);
        	cblas_dcopy(end,XMMA,1,levelset.getVolratio(),1);
			itt++; // update iteration count
        }
        */
        // update shell thk vars
        /*#pragma omp parallel for private(i,j,k,end,sp_ptr3,sp_ptr4)
        for(i=0;i<numShell_var;i++)
        {
            if(fabs(shell_step[i])>1.0e-12)
            {
                shell_thk[i] += shell_step[i]; // update variable
                // update all elements thk associate with the variable
                end = shell_var[i]->getPos();
                for(j=0;j<end;j++)
                {
                    k = shell_var[i]->nums[j]; // element number
                    // change thickness (will multiply sub-matrcies etc.)
                    if(WingMesh.el_type[k]==18){sp_ptr3 = (Shell3*)WingMesh.Elems[k];
                        sp_ptr3->change_thk(shell_thk[i]);
                    }
                    else if(WingMesh.el_type[k]==20){sp_ptr4 = (Shell4*)WingMesh.Elems[k];
                        sp_ptr4->change_thk(shell_thk[i]);
                    }
                }
            }
        }*/
        }

        // clean shape sens & shell sens memory
        if(comp_sens){delete [] comp_sens;}
        if(mass_sens){delete [] mass_sens;}
        if(eig_max_sens){delete [] eig_max_sens;}
        if(vol_sens){delete [] vol_sens;}
        if(eig_ratio_sens){delete [] eig_ratio_sens;}
        if(buck_sens){delete [] buck_sens;}
        if(stress_sens){delete [] stress_sens;}

        if(shell_comp_sens){delete [] shell_comp_sens;}
        if(shell_mass_sens_temp){delete [] shell_mass_sens_temp;}
        if(shell_eig_max_sens){delete [] shell_eig_max_sens;}
        if(shell_eig_ratio_sens){delete [] shell_eig_ratio_sens;}
        if(shell_buck_sens){delete [] shell_buck_sens;}
        if(shell_stress_sens){delete [] shell_stress_sens;}

        if(control.method==1)
        {
			if(err==0 || rewind_count==0) // no rewind - or max rewind passed - EDIT
			{
				rewind_count = 0; // reset count
			// velocity extension
			levelset.fast_march2(1);
			if(control.pinfo==2) {
				sprintf(foname,"%s_velocity_it%i",infile,itt);
				outNodalScalarVTK(&WingMesh, numC8, levelset.getVel(), foname);
			}

			gettimeofday(&clck2,NULL);
			        sum = (clck1.tv_sec+(clck1.tv_usec/1000000.0)) - (clck2.tv_sec+(clck2.tv_usec/1000000.0));
			        timeFile << "\t" << -sum;
			        gettimeofday(&clck1,NULL);

			// lsf upwind gradient computation
			levelset.grad_upwind3();
			if(control.pinfo==2) {
				sprintf(foname,"%s_grad_upwind_it%i",infile,itt);
				outNodalScalarVTK(&WingMesh, numC8, levelset.getGradMag(), foname);
			}

			gettimeofday(&clck2,NULL);
			        sum = (clck1.tv_sec+(clck1.tv_usec/1000000.0)) - (clck2.tv_sec+(clck2.tv_usec/1000000.0));
			        timeFile << "\t" << -sum;
			        gettimeofday(&clck1,NULL);

			// update the lsf (fwd Euler scheme)
			err = levelset.update(); // update the lsf (using the computed velocities and gradients)
			if(err == 1){ reint_cnt = 0;} // re-initialize due to hit mine
			else{reint_cnt--;}  // count down to re-initialization

			itt++; // update iteration count
			}

			//otherwise, update by rewind
			else
			{
				rewind_count++; // up counter
				cout << "\nRewind for iteration " << itt;
				levelset.rewind(0.5);
			}
        }
        te = time(0);
        del_time = difftime(te, ts);
        cout << "\n Iteration took " << del_time << "sec";

        gettimeofday(&clck2,NULL);
        sum = (clck1.tv_sec+(clck1.tv_usec/1000000.0)) - (clck2.tv_sec+(clck2.tv_usec/1000000.0));
        timeFile << "\t" << -sum;
        gettimeofday(&clck1,NULL);

        timeFile.close();
        timeFile.open("timing.txt", ios::out | ios::app);
        timeFile.precision(12);

    } while (itt < control.maxItt);

    // -------------------------------------------------------------------- //
    //                                                                      //
    //      Step 7: Final post-processing                                   //
    //                                                                      //
    // -------------------------------------------------------------------- //

    timeFile.close();
    outFile.close();
    cout << "\n\nEnd of Level-Set 3D v2.3\n\n";

    //return 0;

}
