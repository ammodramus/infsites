#ifndef BMAT_HEADER
#define BMAT_HEADER

#include <stdio.h>
#include <stdint.h>

#define DEFAULT_MAX_LINE_SIZE 10000
#define DEFAULT_MAX_NUM_LINES 1000

typedef struct bmat
{
	int32_t ** mat;
	int32_t nrows;
	int32_t ncols;
} BMat;

void BMat_init(BMat * bmat, int32_t nrows, int32_t ncols);
void BMat_copy_matrix(BMat * bmat, BMat * bmatcopy);
void BMat_transpose_matrix(BMat * bmat, BMat * tbmat);
void BMat_free(BMat * bmat);
void BMat_print(BMat * bmat, FILE * output);
int BMat_int32tcomp_(const void * el1, const void * el2);
int32_t BMat_pow_(int32_t x, int32_t p);
int32_t BMat_get_column_binary_num_(BMat * bmat, int32_t col);
int BMat_compare_columns_(const void * pair1, const void * pair2);
void BMat_order_columns(BMat * bmat);
int32_t BMat_get_Lij(BMat * bmat, int32_t i, int32_t j);
int32_t BMat_get_Lj(BMat * bmat, BMat * Lij, int32_t j);
void BMat_get_haplotype_counts(BMat * bmat, int32_t * numDuplicates, int32_t * numUnique);
int32_t BMat_determine_good(BMat * bmat, BMat * Lij, int32_t * Lj);
int32_t BMat_read_input(FILE * fin, BMat * bmat, int32_t * numHaplotypes);

#endif
