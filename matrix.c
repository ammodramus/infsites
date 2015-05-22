#include <stdlib.h>
#include <stdio.h>
#ifndef CLUSTER
#include <suitesparse/umfpack.h>
#endif
#ifdef CLUSTER
#include <umfpack.h>
#endif
#include "definitions.h"
#include "matrix.h"

void SuperEquations_init(SuperEquations * eq, int nrows, int ncols, int nz)
{
	eq->curEntry = 0;
	eq->nrows = nrows;
	eq->ncols = ncols;
	eq->nz = nz;
	eq->Ti = (int *)malloc(sizeof(int) * nz);
	CHECKPOINTER(eq->Ti);
	eq->Tj = (int *)malloc(sizeof(int) * nz);
	CHECKPOINTER(eq->Tj);
	eq->Tx = (double *)malloc(sizeof(double) * nz);
	CHECKPOINTER(eq->Tx);
	eq->Ap = (int *)malloc(sizeof(int) * (ncols+1));
	CHECKPOINTER(eq->Ap);
	eq->Ai = (int *)malloc(sizeof(int) * nz);
	CHECKPOINTER(eq->Ai);
	eq->Ax = (double *)malloc(sizeof(int) * 2*nz);
	CHECKPOINTER(eq->Ax);
	eq->b = (double *)calloc((size_t)nrows, sizeof(double));
	CHECKPOINTER(eq->b);
	eq->x = (double *)malloc(sizeof(double) * nrows);
	CHECKPOINTER(eq->x);
	eq->symbolic = NULL;
	eq->numeric = NULL;
	return;
}

void SuperEquations_free(SuperEquations * eq)
{
	umfpack_di_free_symbolic(&(eq->symbolic));
	umfpack_di_free_numeric(&(eq->numeric));
	free(eq->Ti);
	free(eq->Tj);
	free(eq->Tx);
	free(eq->Ap);
	free(eq->Ai);
	free(eq->Ax);
	free(eq->b);
	free(eq->x);
	eq->curEntry = 0;
	return;
}

void SuperEquations_add_entry(SuperEquations * eq, int rowIdx, int colIdx, double value)
{
	if(eq->curEntry >= eq->nz)
		PERROR("attempting to add too many entries in matrix.");
	//printf("%i,%i = %f\n", rowIdx, colIdx, value);
	eq->Ti[eq->curEntry] = rowIdx;
	eq->Tj[eq->curEntry] = colIdx;
	eq->Tx[eq->curEntry] = value;
	//printf("i = %i, j = %i, z = %f\n", rowIdx, colIdx, value);
	eq->curEntry++;
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
    int i;
	int status;
	// convert triplet to column format
	status = umfpack_di_triplet_to_col(eq->nrows, eq->ncols, eq->nz, eq->Ti, eq->Tj, eq->Tx, eq->Ap, eq->Ai, eq->Ax, (int *)NULL);
	if(status != UMFPACK_OK)
	{
		fprintf(stderr, "umfpack error status = %i, triplet_to_col conversion error\n", status);
        for(i = 0; i < eq->nz; i++)
            fprintf(stderr, "A[%i, %i] = %f\n", eq->Ti[i], eq->Tj[i], eq->Tx[i]);
        fprintf(stderr, "\n");
		PERROR("UMFPACK ERROR: umfpack_di_triplet_to_col().");
	}
	// now do symbolic factorization.
	status = umfpack_di_symbolic(eq->nrows, eq->ncols, eq->Ap, eq->Ai, eq->Ax, &(eq->symbolic), (double *)NULL, (double *)NULL);
	if(status != UMFPACK_OK)
	{
		printf("umfpack error status = %i, symbolic factorization error\n", status);
		PERROR("UMFPACK ERROR.");
	}

	// numeric factorization
	status = umfpack_di_numeric(eq->Ap, eq->Ai, eq->Ax, eq->symbolic, &(eq->numeric), (double *)NULL, (double *)NULL);
	if(status != UMFPACK_OK)
	{
		printf("umfpack error status = %i, numeric factorization error\n", status);
		PERROR("UMFPACK ERROR.");
	}
	// solving
	status = umfpack_di_solve(UMFPACK_A, eq->Ap, eq->Ai, eq->Ax, eq->x, eq->b, eq->numeric, (double *)NULL, (double *)NULL);
	if(status != UMFPACK_OK)
	{
		fprintf(stderr, "umfpack error status = %i, triplet_to_col conversion error\n", status);
        for(i = 0; i < eq->nz; i++)
            fprintf(stderr, "A[%i, %i] = %f\n", eq->Ti[i], eq->Tj[i], eq->Tx[i]);
        fprintf(stderr, "\n");
		fprintf(stderr, "umfpack error status = %i, solving error\n", status);
		PERROR("UMFPACK ERROR.");
	}
	return;
}
