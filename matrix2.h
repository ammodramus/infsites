#ifndef MATRIX_H
#define MATRIX_H

#include "csparse.h"

typedef struct superequations_
{
	int curEntry;
	cs * A;
	double * b;
	double * x;
} SuperEquations; 

void SuperEquations_init(SuperEquations * eq, int nrows, int ncols, int nz);
void SuperEquations_free(SuperEquations * eq);
void SuperEquations_add_b_value(SuperEquations * eq, double prob, int idx);
void SuperEquations_add_entry(SuperEquations * eq, int rowIdx, int colIdx, double value);
void SuperEquations_solve(SuperEquations * eq);

#endif

