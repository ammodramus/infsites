#ifndef BMAT_HEADER
#define BMAT_HEADER

#include <stdio.h>
#include <stdint.h>

#define DEFAULT_MAX_LINE_SIZE 10000
#define DEFAULT_MAX_NUM_LINES 1000

typedef struct bmat
{
	int ** mat;
	int nrows;
	int ncols;
} BMat;

void BMat_init(BMat * bmat, int nrows, int ncols);
void BMat_copy_matrix(BMat * bmat, BMat * bmatcopy);
void BMat_transpose_matrix(BMat * bmat, BMat * tbmat);
void BMat_free(BMat * bmat);
void BMat_print(BMat * bmat, FILE * output);
int BMat_int32tcomp_(const void * el1, const void * el2);
int BMat_pow_(int x, int p);
int BMat_get_column_binary_num_(BMat * bmat, int col);
int BMat_compare_columns_(const void * pair1, const void * pair2);
void BMat_order_columns(BMat * bmat);
int BMat_get_Lij(BMat * bmat, int i, int j);
int BMat_get_Lj(BMat * bmat, BMat * Lij, int j);
void BMat_get_haplotype_counts(BMat * bmat, int * numDuplicates, int * numUnique);
int BMat_determine_good(BMat * bmat, BMat * Lij, int * Lj);
int BMat_read_input(FILE * fin, BMat * bmat, int * numHaplotypes);

#endif
