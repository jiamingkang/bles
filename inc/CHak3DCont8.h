/*
 * CHak3DCont8.h
 *
 *  Created on: 23 Feb 2015
 *      Author: knai20-admin
 */

#ifndef CHAK3DCONT8_H_
#define CHAK3DCONT8_H_

#include "inc/CHak3DElement.h"

static const Coord3D c8pos3[27] = {
    {-1,-1,0},
    {0,-1,0},
    {1,-1,0},
    {-1,0,0},
    {0,0,0},
    {1,0,0},
    {-1,1,0},
    {0,1,0},
    {1,1,0},
    {-1,-1,-1},
    {0,-1,-1},
    {1,-1,-1},
    {-1,0,-1},
    {0,0,-1},
    {1,0,-1},
    {-1,1,-1},
    {0,1,-1},
    {1,1,-1},
    {-1,-1,1},
    {0,-1,1},
    {1,-1,1},
    {-1,0,1},
    {0,0,1},
    {1,0,1},
    {-1,1,1},
    {0,1,1},
    {1,1,1},
};

static const Coord3D c8_pos[8] = { {-1,-1,-1},
    {1,-1,-1},
    {1,1,-1},
    {-1,1,-1},
    {-1,-1,1},
    {1,-1,1},
    {1,1,1},
    {-1,1,1} };

// convention for faces of a non-isoparamtic C8 element
// can access using: index(face,node,4)
static const int c8faces[24] = {
    0,1,2,3, // face 1 at z = -1
    0,3,7,4, // face 2 at x = -1
    0,1,5,4, // face 3 at y = -1
    3,2,6,7, // face 4 at y = +1
    1,2,6,5, // face 5 at x = +1
    4,5,6,7  // face 6 at z = +1
};

// array storing C8 nodal connectivity
// each node is connected to 3 others
static const int c8nodes[24] = {
    1,3,4, // nodes connected to node 0
    0,2,5, // nodes connected to node 1
    1,3,6, // nodes connected to node 2
    0,2,7, // nodes connected to node 3
    0,5,7, // nodes connected to node 4
    1,4,6, // nodes connected to node 5
    2,5,7, // nodes connected to node 6
    3,4,6  // nodes connected to node 7
};

// array storing the opposite diagonal nodes
static const int c8opp[8] = {6,7,4,5,2,3,0,1};

static const double gauss3wgt[27] = {
    200.0, 320.0, 200.0, 320.0, 512.0, 320.0, 200.0, 320.0, 200.0,
    125.0, 200.0, 125.0, 200.0, 320.0, 200.0, 125.0, 200.0, 125.0,
    125.0, 200.0, 125.0, 200.0, 320.0, 200.0, 125.0, 200.0, 125.0,
};

static const double gauss2 = 0.577350269189626;
static const double gauss3 = 0.774596669241483;

// Eight node 3D continuum element (brick)

class CHak3DCont_8: public CHak3DElement {
public:
	CHak3DCont_8();
	virtual ~CHak3DCont_8();

public:

	// create stiffness matrix
	void makeKe();
	// create mass matrix
	void makeMe();
	// get element displacements
	void getDisp(double*, double*);
	// compute C8 strain tensor at a point in the element
	void strain(Coord3D, double*, double*);
	// compute C8 stress tensor at a point in the element
	void stress(Coord3D, double*, double*);
    // compute element volume
    void makeVol();
    // add self-weight to a global load vector
    void addWeight(Coord3D, double, int*, double*);
    // compute C8 sensiticies at all gauss points
    void shpSens(int, int, double**, double**, double*, double*, double*, int, double*, double);
    // compute C8 sensitivities at all gauss points, for eigenvalue problems
    void shpSens_Eig(int numDual, int numEig, double **disp_prim, double **disp_dual,
                     double *eig, double *wgt_case, double *gSens, int pNum, double alpha);
	// output
	void KeOut() {
		int i;
		for(i=0;i<300;i++) {
			std::cout << "\n" << i << " - " << Ke[i];
		}
	}
	void MeOut() {
		int i;
		for(i=0;i<300;i++) {
			std::cout << "\n" << i << " - " << Me[i];
		}
	}
};

#endif /* CHAK3DCONT8_H_ */
