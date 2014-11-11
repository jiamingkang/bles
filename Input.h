/*
 *  Input.h
 *
 *  Created by Peter Dunning
 *	Version 5.4 @ 23/04/2014
 */

#define MAX_CHARS_PER_LINE 512
#define MAX_TOKENS_PER_LINE 16
#define MAX_HOLES 500

#include "ls_types.h"

int icmpfunc (const void * p1, const void * p2);

// function to read new-style input file
int read_input(char *datafile, mesh *inMesh, int *numMat, isoMat *inMat, levSet *levelset, prob *lsprob,
			   ctrl *control, int **fixDof, int *numCase, double **load, int *freeDof, sp_mat *lump_mass, bool *sw, Coord **acc);
