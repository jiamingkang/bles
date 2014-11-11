/*
 *  Numbering.h
 *
 *  Created by Peter Dunning
 *	Version 5.4 @ 23/04/2014
 */

#include "ls_types.h"

// Function that numbers all elements and nodes in the FG domain
void Numbering(mesh *inMesh);

// Node Co-ordinate calculation function
void Coordinates(mesh *inMesh);

// Function that orders node numbers into a 2D based on their relative positions
void NodeNums2(mesh *inMesh);

// function to number bars elements
void Bar_numbering(mesh *inMesh);

// function to free memory for a sparse matrix
void free_sp_mat(sp_mat *m);

// function to create memory for a sparse matrix
void set_sp_mat(sp_mat *m);

// function to remove dof from a sparse matrix
void rem_sp_mat(sp_mat *m, int *map, int inc);
