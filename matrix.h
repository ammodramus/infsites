#ifndef MATRIX_H
#define MATRIX_H

typedef struct superequations_
{
	int curEntry;
	int nrows;
	int ncols;
	int nz;
	int * Ti;
	int * Tj;
	double * Tx;
	int * Ap;
	int * Ai;
	double * Ax;
	double * b;
	double * x;
	void * symbolic;
	void * numeric;
} SuperEquations; 

void SuperEquations_init(SuperEquations * eq, int nrows, int ncols, int nz);
void SuperEquations_free(SuperEquations * eq);
void SuperEquations_add_b_value(SuperEquations * eq, double prob, int idx);
void SuperEquations_add_entry(SuperEquations * eq, int rowIdx, int colIdx, double value);
void SuperEquations_solve(SuperEquations * eq);

#endif
