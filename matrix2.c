#include <stdlib.h>
#include <stdio.h>
//#include <suitesparse/umfpack.h>
#include "definitions.h"
#include "matrix2.h"
#include "csparse.c"

void SuperEquations_init(SuperEquations * eq, int nrows, int ncols, int nz)
{
	eq->A = cs_spalloc(nrows, ncols, nz, nz, 1);
	eq->b = (double *)calloc((size_t)nrows, sizeof(double));
	CHECKPOINTER(eq->b);
	//eq->x = (double *)malloc(sizeof(double) * nrows);
	//CHECKPOINTER(eq->x);
	return;
}

void SuperEquations_free(SuperEquations * eq)
{
	cs_spfree(eq->A);
	free(eq->b);
	//free(eq->x);
	return;
}

void SuperEquations_add_entry(SuperEquations * eq, int rowIdx, int colIdx, double value)
{
	//printf("%i,%i = %f\n", rowIdx, colIdx, value);
	cs_entry(eq->A, rowIdx, colIdx, value);
	return;
}

void SuperEquations_add_b_value(SuperEquations * eq, double prob, int idx)
{
	eq->b[idx] += prob;
	//printf("b[%i] = %f\n", idx, eq->b[idx]);
	return;
}

/* TODO: write functions to fill in the triplet forms of the A matrices using
 * the (migration) coalescent theory knowledge. */

/* this function solves the system of equations using UMFPACK. it assumes that
 * the matrix has been inputted in triplet form; the first thing it does is
 * convert to column form. */
void SuperEquations_solve(SuperEquations * eq)
{
	// convert to compressed column format
	eq->A = cs_triplet(eq->A);
	// use cholsol
	//if(!cs_cholsol(eq->A, eq->b, 0))
	if(!cs_qrsol(eq->A, eq->b, 0))
		PERROR("Cholesky solve error (PW)");
	return;
}

