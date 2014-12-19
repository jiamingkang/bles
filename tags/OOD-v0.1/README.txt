Bath LEvel Set (BLES) version 5.30
Written by Dr. Peter D. Dunning - 7th January 2014
Bugs, questions, praise to: en2pdd@gmail.com

--------------------------------------------------------------------------------------
-- 01/12/2014 jeehanglee@gmail.com
	Not sure this code can be compiled on both 32- and 64-bit machine
	is it depending on compilers or the system itself?
	my machine is based on 32-bit, thus OS is 32-bit Snow-Leopard (MacOS 10.6.8)
	consequently gcc is 32-bit version, which brings about compile errors
	when calling functions in lapack library.
	-- Just upload and compile on Khalid's 64-bit machine.
--------------------------------------------------------------------------------------

This code is not guaranteed to be bug free or work on your system. You are free to modify this code for research purposes, with proper acknowledgment. The main references for the code are:

Dunning, P.D., Kim, H.A., 2014. "Level set based optimization using sequential linear programming," 11th World Congress on Computational Mechanics, Barcelona, Spain, 20-25 July. 

Dunning, P.D., Kim, H.A., Mullineux, G., 2011. “Investigation and improvement of sensitivity computation using the area-fraction weighted fixed grid FEM and structural optimization,” Finite Elements in Analysis and Design, 47(8): 933-941. (DOI: 10.1016/j.finel.2011.03.006)

Dunning, P., 2011. Introducing Loading Uncertainty in Level Set-Based Structural Topology Optimisation. Thesis (Doctor of Philosophy (PhD)). University of Bath.
(http://opus.bath.ac.uk/27508)

Major updates since version 5.20:
* New options for the input file:
	*fix-lsf : used to define nodes or regions where the lsf is fixed, 
		this replaces the old method of using *bound to fix the lsf
	*def_mat : used to define the material for a region,
		now up to 5 different materials can be specified and used
	*mass : used to add lumped masses to the FE model
* Eigenvalue solve. A general eigenvalue solver has been implemented
	Consistent element mass matrices used & lumped mass can be added
	5.30 can solve two natural frequency problems (see below)
	Warning! assumes that frequencies are simple and not repeated
	Eigenvalue solver is from ARPACK (required library - see below)
	One additional control - minimum area ratio for mass
* Multiple constraint handling. Improved version of SLP level-set method
	5.30 can handle 4 different objectives and 5 constraints
	Not all objective / constraint combinations are valid (see below)
	The LP sub-solve is now handled by the simplex method (from HSL - new required library)
* Converted most of the output to .vtk format to be compatible with Paraview 

The current code can solve the following optimization problems:
* Objective functions:
	1) Minmise Compliance (inc. multiple load cases)
		- can use volume & target displacement constraints
	2) Minmise Volume
		- can use compliance & target displacement constraints
	3) Maximise 1st natural frequency
		- can use volume and eigen-frequency ratio constraints
	4) Maximise Mechanical advantage (compliant mechanism design)
		- MUST use volume and an input displacement constraint
* Constraints:
	1) Maximum Volume (specified as % of design domain)
	2) Maximum Compliance
	3) Target displacement (solved as minimise least squares difference to target)
	4) Minimum Eigen frequency ratio (specified as 2nd / 1st eigen frequency)
	5) Maximum input displacement (must be included for compliant mechanism design)

For examples of how to specify each constraint see the example.txt input file

Main limitations of the code:
* 2D problems
* Elements are plane four node continuum type
* Rectangular design domain with square elements
* Isotropic material (plane stress assumption)

Required libraries:
(note that blas and lapack are often included as standard on unix type systems)
* BLAS & C BLAS, www.netlib.org/blas 
* LAPACK, www.netlib.org/lapack
* HSL MA57 (double precision version), www.hsl.rl.ac.uk/catalogue/ma57.xml (can be obtained free on an academic licence, note that BLES was developed and tested with the older MA57 code, not the newer HSL_MA57)
* HSL LA01 (revised simplex method), www.hsl.rl.ac.uk/archive/specs/la01.pdf (note may need an extra component from HSL - fd05)
* ARPACK (eigenvalue solver), www.caam.rice.edu/software/ARPACK/ (I used the serial FORTRAN version)

Tested platforms:
* Only tested on MAC OSX 10.6.8 using llvm and gfortran compilers.
* I had problems (unresolved segmentation fault) when compiling with gcc, although it may be just my setup

Building the code:
* You will need to obtain and install the required libraries mentioned above (if you do not already have them).
* There is an example makefile included with the source code, this will likely have to be customised for your system (sorry, doing clever auto build scripts is currently beyond me).

Running the code:
* After the code is built run the program by typing: ./BLESV5 input
  where input is the name of an input file (see examples)

Examples:
* Several examples are included. These are:
* The first 4 are minimisation of compliance problems, s.t. maximum volume constraint
	1) Cantilever 2:1 ratio (cant.txt)
	2) MBB beam, right half of 12:1 ratio (mbb.txt)
	3) Michell arch structure (mich.txt)
	4) L bracket (lbracket.txt)
* Four other problem examples
	5) Cantilever with displacement constraint (cantD.txt)
	6) Cantilever with compliance constraint (cant_vol.txt)
	7) Crunching compliant mechanism (crunch.txt)
	8) Frequency maximisation (freq_test2.txt)
	
