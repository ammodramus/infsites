#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include "definitions.h"
#include "bmat.h"

void BMat_init(BMat * bmat, int nrows, int ncols)
{
	int i, j;
	bmat->nrows = nrows;
	bmat->ncols = ncols;
	bmat->mat = (int **)malloc(sizeof(int *) * (size_t)nrows);
	CHECKPOINTER(bmat->mat);
	for(i = 0; i < nrows; i++)
	{
		bmat->mat[i] = (int *)malloc(sizeof(int) * (size_t)ncols);
		CHECKPOINTER(bmat->mat[i]);
		for(j = 0; j < ncols; j++)
			bmat->mat[i][j] = 0;
	}
	return;
}

/* assumes bmatcopy is a pointer to an uninitialized BMat */
void BMat_copy_matrix(BMat * bmat, BMat * bmatcopy)
{
	int i, j, nrows = bmat->nrows, ncols = bmat->ncols;
	BMat_init(bmatcopy, nrows, ncols);
	for(i = 0; i < nrows; i++)
		for(j = 0; j < ncols; j++)
			bmatcopy->mat[i][j] = bmat->mat[i][j];
	return;
}

/* assumes tbmat is a pointer to an uninitialized BMat */
void BMat_transpose_matrix(BMat * bmat, BMat * tbmat)
{
	int i, j, nrows = bmat->nrows, ncols = bmat->ncols;
	BMat_init(tbmat, ncols, nrows);
	for(i = 0; i < nrows; i++)
		for(j = 0; j < ncols; j++)
			tbmat->mat[j][i] = bmat->mat[i][j];
	return;
}

void BMat_free(BMat * bmat)
{
	int i;
	for(i = 0; i < bmat->nrows; i++)
		free(bmat->mat[i]);
	free(bmat->mat);
	return;
}

void BMat_print(BMat * bmat, FILE * output)
{
	int i, j;
	for(i = 0; i < bmat->nrows; i++)
	{
		for(j = 0; j < bmat->ncols; j++)
			fprintf(output, "%i", bmat->mat[i][j]);
		fprintf(output, "\n");
	}
	return;
}

int BMat_int32tcomp_(const void * el1, const void * el2)
{
	int v1 = *((int *)el1);
	int v2 = *((int *)el2);
	if(v1 > v2)
		return 1;
	if(v2 > v1)
		return -1;
	return 0;
}

int BMat_pow_(int x, int p)
{
	if (p == 0) return 1;
	if (p == 1) return x;
	int tmp = BMat_pow_(x, p/2);
	if (p % 2 == 0)
		return tmp * tmp;
	else return x * tmp * tmp;
}

int BMat_get_column_binary_num_(BMat * bmat, int col)
{
	int i, digit = 0, binNum = 0;
	for(i = bmat->nrows-1; i >= 0; i--)
	{
		if(bmat->mat[i][col])
			binNum += BMat_pow_((int)2, digit);
		digit++;
	}
	return binNum;
}

int BMat_compare_columns_(const void * pair1, const void * pair2)
{
	int v1, v2;
	v1 = *(int *)((*(void ***)pair1)[0]);
	v2 = *(int *)((*(void ***)pair2)[0]);
	if(v1 < v2)
		return 1;
	if(v2 < v1)
		return -1;
	return 0;
}

void BMat_order_columns(BMat * bmat)
{
	int i, j;
	int nrows = bmat->nrows;
	int ncols = bmat->ncols;
	int * binNums = (int *)malloc(sizeof(int) * (size_t)ncols);
	CHECKPOINTER(binNums);
	for(j = 0; j < ncols; j++)
		binNums[j] = BMat_get_column_binary_num_(bmat, j);
	BMat tbmat;
	BMat_transpose_matrix(bmat, &tbmat);
	void *** numcolpairs = (void ***)malloc(sizeof(void **) * (size_t)ncols);
	CHECKPOINTER(numcolpairs);
	for(j = 0; j < ncols; j++)
	{
		numcolpairs[j] = (void **)malloc(sizeof(void *) * (size_t)2);
		CHECKPOINTER(numcolpairs[j]);
		numcolpairs[j][0] = (void *)(&(binNums[j]));
		numcolpairs[j][1] = (void *)(tbmat.mat[j]);
	}
	// the sort
	qsort(numcolpairs, (size_t)ncols, sizeof(void **), BMat_compare_columns_);
	BMat bmatordered;
	BMat_init(&bmatordered, nrows, ncols);
	for(j = 0; j < ncols; j++)
		for(i = 0; i < nrows; i++)
			bmatordered.mat[i][j] = ((int *)(numcolpairs[j][1]))[i];
	for(j = 0; j < ncols; j++)
		for(i = 0; i < nrows; i++)
			bmat->mat[i][j] = bmatordered.mat[i][j];
	// cleaning up.
	BMat_free(&bmatordered);
	BMat_free(&tbmat);
	free(binNums);
	for(j = 0; j < ncols; j++)
		free(numcolpairs[j]);
	free(numcolpairs);
	return;
}

int BMat_get_Lij(BMat * bmat, int i, int j)
{
	int k;
	for(k = j-1; k >= 0; k--)
		if(bmat->mat[i][k])
			return k+1;		// assume 1-based indexing.
	return 0;
}

int BMat_get_Lj(BMat * bmat, BMat * Lij, int j)
{
	int i;
	int Lj = 0;
	for(i = 0; i < bmat->nrows; i++)
	{
		if(bmat->mat[i][j])
			if(Lij->mat[i][j] > Lj)
				Lj = Lij->mat[i][j];
	}
	return Lj;
}

int BMat_determine_good(BMat * bmat, BMat * Lij, int * Lj)
{
	int i,j;
	for(i = 0; i < bmat->nrows; i++)
		for(j = 0; j < bmat->ncols; j++)
			if(bmat->mat[i][j] && (Lij->mat[i][j] != Lj[j]))
				return 0;
	return 1;
}

int BMat_read_input(FILE * inp, BMat * bmat, int * numHaplotypes)
{
	int i, j, numSegSites, nrows = 0;
	char line[DEFAULT_MAX_LINE_SIZE];
	char ** lines = (char **)malloc(sizeof(char *) * DEFAULT_MAX_NUM_LINES);
	CHECKPOINTER(lines);
	lines[0] = (char *)malloc(sizeof(char) * DEFAULT_MAX_LINE_SIZE);
	CHECKPOINTER(lines[0]);
	if(fgets(line, DEFAULT_MAX_LINE_SIZE, inp) == NULL)
		PERROR("No first line in input.");
	numSegSites = (int)strlen(line)-1;
	strcpy(lines[nrows++], line);
	while(fgets(line, DEFAULT_MAX_LINE_SIZE, inp) != NULL)
	{
		lines[nrows] = (char *)malloc(sizeof(char) * DEFAULT_MAX_LINE_SIZE);
		CHECKPOINTER(lines[nrows]);
		strcpy(lines[nrows++], line);
	}
	BMat_init(bmat, nrows, numSegSites);
    int mono = 1;
	for(i = 0; i < nrows; i++)
	{
		for(j = 0; j < numSegSites; j++)
		{
			//fprintf(stdout, "%c", lines[i][j]);
			if(lines[i][j] == '\0' || lines[i][j] == '\n')
				continue;
			if(lines[i][j] != '0' && lines[i][j] != '1')
			{
				fprintf(stdout, "%c", lines[i][j]);
				PERROR("non 0/1 character in input.");
			}
            if(lines[i][j] == '1')
                mono = 0;
			bmat->mat[i][j] = lines[i][j] - '0';
		}
	}
	//BMat_print(bmat, stdout);
	for(i= 0; i < nrows; i++)
		free(lines[i]);
	free(lines);
	fclose(inp);
    *numHaplotypes = nrows;
	return mono;
}

int BMat_read_input_ctypes(char ** inp, BMat * bmat, int numHaplotypes)
{
	int i, j, numSegSites;
	char * line;

    for(i = 0; i < numHaplotypes; i++)
    {
        if(!inp[i])
        {
            fprintf(stderr, "missing line: %i\n", i);
            PERROR("missing input in BMat_read_input_ctypes()");
        }
    }

    numSegSites = (int)strlen(inp[0]);
	
	BMat_init(bmat, numHaplotypes, numSegSites);
    int mono = 1;
	for(i = 0; i < numHaplotypes; i++)
	{
        line = inp[i];
		for(j = 0; j < numSegSites; j++)
		{
			if(line[j] != '0' && line[j] != '1')
			{
				fprintf(stdout, "%c", line[j]);
				PERROR("non 0/1 character in input.");
			}
            if(line[j] == '1')
                mono = 0;
			bmat->mat[i][j] = line[j] - '0';
		}
	}
	return mono;
}

/* This is a fairly complicated way of looking at a number of rows of things
 * and determining the counts of each unique row. It makes every comparison
 * between rows i and ii except when ii has already been matched to something
 * else. 
 *
 * 0 0 0
 * 0 1 0
 * 0 0 0 
 * 0 1 1
 * 0 1 0
 *
 * The above BMat would return (2 2 1) for numDuplicates and 3 for *numUnique.
 *
 * numDuplicates is an input parameter, an array that is at least bmat->nrows long.
 * numUnique is an input parameter, the pointer to an int. */
void BMat_get_haplotype_counts(BMat * bmat, int * numDuplicates, int * numUnique)
{
	int i, ii, j, same, nrows = bmat->nrows, ncols = bmat->ncols, dupCounter = 0;
	int * checked = (int *)malloc(sizeof(int) * (size_t)bmat->nrows);
	CHECKPOINTER(checked);
	//fprintf(stdout, "nrows = %i\n", bmat->nrows);
	for(i = 0; i < nrows; i++)
	{
		checked[i] = 0;
		numDuplicates[i] = 1;
	}
	for(i = 0; i < nrows; i++)
	{
		//fprintf(stdout, "i = %i\n", i);
		if(checked[i])
			continue;
		checked[i] = 1;
		for(ii = i+1; ii < nrows; ii++)
		{
			//fprintf(stdout, "\tii = %i\n", ii);
			if(checked[ii])
			{
				//printf("\t(continuing...)\n");
				continue;
			}
			same = 1;
			for(j = 0; j < ncols; j++)
			{
				if(bmat->mat[i][j] != bmat->mat[ii][j])
				{
					same = 0;
					break;
				}
			}
			if(same)
			{
				//printf("\tsame\n");
				checked[ii] = 1;
				numDuplicates[dupCounter]++;
			}
			//else
				//printf("\tdifferent\n");
		}
		//printf("increasing dupCounter...\n");
		dupCounter++;
	}
	*numUnique = dupCounter;
	//*numUnique = dupCounter+1;
	free(checked);
	return;
}


/* test main */
/*
int main(int argc, char ** argv)
{
	int i,j;
	int nrows = 4;
	int ncols = 4;
	BMat bmat;
	BMat_init(&bmat, nrows, ncols);
	bmat.mat[1][1] = 1;
	bmat.mat[2][1] = bmat.mat[2][0] = 1;
	bmat.mat[3][1] = bmat.mat[3][2] = 1;
	for(i = 0; i < nrows; i++)
		bmat.mat[i][3] = bmat.mat[i][0];
	printf("\nat first\n");
	BMat_print(&bmat, stdout);
	int * dupcounts = (int *)malloc(sizeof(int)*(size_t)ncols);
	CHECKPOINTER(dupcounts);
	BMat_order_columns(&bmat);
	printf("\nafter sorting\n");
	BMat_print(&bmat, stdout);
	ncols = bmat.ncols;
	printf("\nafter deleting\n");
	BMat_print(&bmat, stdout);
	printf("newncols: %i\n", ncols);
	printf("dupcounts: ");
	for(j = 0; j < ncols; j++)
		printf("%i ", dupcounts[j]);
	printf("\n");
	BMat Lij;
	BMat_init(&Lij, nrows, ncols);
	for(i = 0; i < nrows; i++)
		for(j = 0; j < ncols; j++)
			Lij.mat[i][j] = BMat_get_Lij(&bmat, i, j);
	printf("\nLij:\n");
	int * Lj = (int *)malloc(sizeof(int) * (size_t)ncols);
	CHECKPOINTER(Lj);
	for(j = 0; j < ncols; j++)
	{
		Lj[j] = BMat_get_Lj(&bmat, &Lij, j);
		printf("%i ", Lj[j]);
	}
	printf("\nLij\n");
	BMat_print(&Lij, stdout);
	printf("\ngood: %i\n", BMat_determine_good(&bmat, &Lij, Lj));
	BMat_free(&Lij);
	BMat_free(&bmat);
	return 0;
}
*/
