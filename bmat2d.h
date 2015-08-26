#ifndef BMAT2D_H
#define BMAT2D_H

typedef struct bmat2d_ 
{
	BMat bmat;
	int * demes;
	int nrows;
	int ncols;
} BMat2d; 

void BMat2d_init(BMat2d * b2, int nrows, int numSegSites);
void BMat2d_free(BMat2d * b2);
void BMat2d_read_input(FILE * fin, BMat2d * bmat2d);
void BMat2d_read_input_ctypes(char ** inp, int numHaplotypes, BMat2d * bmat2d);
void BMat2d_get_haplotype_counts(BMat2d * bmat2d, twoints * numDuplicates, int * numUnique);

#endif
