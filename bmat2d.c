#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include "definitions.h"
#include "bmat.h"
#include "bmat2d.h"

void BMat2d_init(BMat2d * b2, int32_t nrows, int32_t numSegSites)
{
	BMat_init(&(b2->bmat), nrows, numSegSites);
	b2->demes = (int32_t *)malloc(sizeof(int32_t) * nrows);
	b2->nrows = nrows;
	b2->ncols = numSegSites;
	CHECKPOINTER(b2->demes);
	return;
}

void BMat2d_free(BMat2d * b2)
{
	BMat_free(&(b2->bmat));
	free(b2->demes);
	return;
}

void BMat2d_read_input(FILE * inp, BMat2d * bmat2d)
{
	int32_t i, j, deme, numSegSites, nrows = 0;
	char line[DEFAULT_MAX_LINE_SIZE];
	char haplotype[DEFAULT_MAX_LINE_SIZE];
	char ** lines = (char **)malloc(sizeof(char *) * DEFAULT_MAX_NUM_LINES);
	CHECKPOINTER(lines);
	int32_t * demes = (int32_t *)malloc(sizeof(int32_t) * DEFAULT_MAX_NUM_LINES);
	CHECKPOINTER(demes);
	lines[nrows] = (char *)malloc(sizeof(char) * DEFAULT_MAX_LINE_SIZE);
	CHECKPOINTER(lines[nrows]);
	if(fgets(line, DEFAULT_MAX_LINE_SIZE, inp) == NULL)
		PERROR("No first line in input.");
	sscanf(line, "%s %i\n", haplotype, &deme);
	numSegSites = (int32_t)strlen(haplotype);
	strcpy(lines[nrows], line);
	demes[nrows++] = deme;
    int32_t mono = 1;
	while(fgets(line, DEFAULT_MAX_LINE_SIZE, inp) != NULL)
	{
		lines[nrows] = (char *)malloc(sizeof(char) * DEFAULT_MAX_LINE_SIZE);
		CHECKPOINTER(lines[nrows]);
		sscanf(line, "%s %i\n", haplotype, &deme); // *space* between %s and %i
		if(deme != 0 && deme != 1)
			PERROR("Non 0/1 deme entry in two-deme input.");
		strcpy(lines[nrows], haplotype);
		demes[nrows++] = deme;
        for(i = 0; i < numSegSites; i++)
        {
            if(haplotype[i] == '1')
                mono = 0;
        }
	}
    if(!mono)
        BMat2d_init(bmat2d, nrows, numSegSites);
    else // mono (numSegSites-1 because monomorphic sites are inputted as a single column of 0's, which would otherwise be counted as a segregating site)
        BMat2d_init(bmat2d, nrows, numSegSites-1);
	for(i = 0; i < nrows; i++)
	{
		for(j = 0; j < numSegSites; j++)
		{
			if(lines[i][j] == '\0' || lines[i][j] == '\n')
				continue;
			if(lines[i][j] != '0' && lines[i][j] != '1')
			{
				fprintf(stdout, "%c", lines[i][j]);
				PERROR("non 0/1 character in input.");
			}
            if(lines[i][j] == '1')
                mono = 0;
			bmat2d->bmat.mat[i][j] = lines[i][j] - '0';
		}
		bmat2d->demes[i] = demes[i];
	}
	for(i= 0; i < nrows; i++)
		free(lines[i]);
	free(lines);
	free(demes);
	fclose(inp);
	return;
}

/* This is a fairly complicated way of looking at a number of rows of things
 * and determining the counts of each unique row. It makes every comparison
 * between rows i and ii except when ii has already been matched to something
 * else. It makes comparisons separately for each deme.
 *
 * 0 0 0
 * 0 1 0
 * 0 0 0 
 * 0 1 1
 * 0 1 0
 *
 * The above BMat would return (2 2 1) for numDuplicates and 3 for *numUnique.
 *
 * numDuplicates is an input parameter, an array of twoints that is at least bmat->nrows long.
 * numUnique is an input parameter, the pointer to an int32_t array of length 2. */
void BMat2d_get_haplotype_counts(BMat2d * bmat2d, twoints * numDuplicates, int32_t * numUnique)
{
	int32_t i, ii, j, deme, same, nrows = bmat2d->bmat.nrows, ncols = bmat2d->bmat.ncols, dupCounter;
	int32_t * demes = bmat2d->demes;
	int32_t * checked = (int32_t *)malloc(sizeof(int32_t) * (size_t)bmat2d->bmat.nrows);
	CHECKPOINTER(checked);
	//fprintf(stdout, "nrows = %i\n", bmat->nrows);
	for(deme = 0; deme < 2; deme++)
	{
		dupCounter = 0;
		for(i = 0; i < nrows; i++)
		{
			checked[i] = 0;
			numDuplicates[i][deme] = 1;
		}
		for(i = 0; i < nrows; i++)
		{
			if(demes[i] != deme)
				continue;
			//fprintf(stdout, "i = %i\n", i);
			if(checked[i])
				continue;
			checked[i] = 1;
			for(ii = i+1; ii < nrows; ii++)
			{
				if(demes[ii] != deme)
					continue;
				//fprintf(stdout, "\tii = %i\n", ii);
				if(checked[ii])
				{
					//printf("\t(continuing...)\n");
					continue;
				}
				same = 1;
				for(j = 0; j < ncols; j++)
				{
					if(bmat2d->bmat.mat[i][j] != bmat2d->bmat.mat[ii][j])
					{
						same = 0;
						break;
					}
				}
				if(same)
				{
					//printf("\tsame\n");
					checked[ii] = 1;
					numDuplicates[dupCounter][deme]++;
				}
				//else
					//printf("\tdifferent\n");
			}
			//printf("increasing dupCounter...\n");
			dupCounter++;
		}
		numUnique[deme] = dupCounter;
	}
	//*numUnique = dupCounter+1;
	free(checked);
	return;
}

/* test main */

/*
int main()
{
	int i;
	BMat2d b2;
	BMat2d_read_input("testfile2d", &b2);
	BMat_print(&(b2.bmat), stdout);
	for(i = 0; i < b2.nrows; i++)
		REPORTI(b2.demes[i]);
	return 0;
}
*/
