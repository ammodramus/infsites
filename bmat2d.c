#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include "definitions.h"
#include "bmat.h"
#include "bmat2d.h"

void BMat2d_init(BMat2d * b2, int nrows, int numSegSites)
{
	BMat_init(&(b2->bmat), nrows, numSegSites);
	b2->demes = (int *)malloc(sizeof(int) * nrows);
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
	int i, j, deme, numSegSites, nrows;
	char line[DEFAULT_MAX_LINE_SIZE];
	char haplotype[DEFAULT_MAX_LINE_SIZE];
	char ** lines = (char **)malloc(sizeof(char *) * DEFAULT_MAX_NUM_LINES);
	CHECKPOINTER(lines);
	int * demes = (int *)malloc(sizeof(int) * DEFAULT_MAX_NUM_LINES);
	CHECKPOINTER(demes);

    int mono = 1;

    lines[0] = (char *)malloc(sizeof(char) * DEFAULT_MAX_LINE_SIZE);
    CHECKPOINTER(lines[0]);

	if(fgets(line, DEFAULT_MAX_LINE_SIZE, inp) == NULL)
		PERROR("No first line in input.");

    sscanf(line, "%s %i\n", haplotype, &deme); // *space* between %s and %i
    numSegSites = (int)strlen(haplotype);
    if(deme != 0 && deme != 1)
        PERROR("Non 0/1 deme entry in two-deme input.");
    strcpy(lines[0], haplotype);
    demes[0] = deme;

    // now proceed on to processing second row
    nrows = 1;

    for(i = 0; i < numSegSites; i++)
    {
        if(haplotype[i] == '1')
            mono = 0;
    }

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
				fprintf(stderr, "bad character: %c", lines[i][j]);
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

void BMat2d_read_input_ctypes(char ** inp, int * demes, int numHaplotypes, BMat2d * bmat2d)
{
	int i, j, deme, numSegSites;
	char * line;

    int mono = 1;

    // first process first line to get number of segregating sites

    for(i = 0; i < numHaplotypes; i++)
    {
        if(inp[i] == NULL)
        {
            fprintf(stderr, "missing line: %i\n", i);
            PERROR("missing input in BMat2d_read_input_ctypes()");
        }
    }

    // first, look at first line to determine number of segregating sites
    line = inp[0];
    numSegSites = (int)strlen(line);
    REPORTI(numSegSites);
    // now determine whether monomorphic or not
    
    for(i = 0; i < numHaplotypes; i++)
    {
        for(j = 0; j < numSegSites; j++)
        {
            if(inp[i][j] == '1')
                mono = 0;
        }
    }

    if(!mono)
        BMat2d_init(bmat2d, numHaplotypes, numSegSites);
    else // mono (numSegSites-1 because monomorphic sites are inputted as a single column of 0's, which would otherwise be counted as a segregating site)
        BMat2d_init(bmat2d, numHaplotypes, numSegSites-1);

	for(i = 0; i < numHaplotypes; i++)
	{
        line = inp[i];
        deme = demes[i];
		for(j = 0; j < numSegSites; j++)
		{
			if(line[j] == '\0' || line[j] == '\n')
				continue;
			if(line[j] != '0' && line[j] != '1')
			{
				fprintf(stderr, "bad character: %c", line[j]);
				PERROR("non 0/1 character in input.");
			}
            // (from a recent C standard:)
            // In both the source and execution basic character sets, the value
            // of each character after 0 in the [...] list of decimal digits
            // shall be one greater than the value of the previous.
			bmat2d->bmat.mat[i][j] = line[j] - '0';
		}
		bmat2d->demes[i] = deme;
	}
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
 * numUnique is an input parameter, the pointer to an int array of length 2. */
void BMat2d_get_haplotype_counts(BMat2d * bmat2d, twoints * numDuplicates, int * numUnique)
{
	int i, ii, j, deme, same, nrows = bmat2d->bmat.nrows, ncols = bmat2d->bmat.ncols, dupCounter;
	int * demes = bmat2d->demes;
	int * checked = (int *)malloc(sizeof(int) * (size_t)bmat2d->bmat.nrows);
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
