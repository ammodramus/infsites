#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <stdio.h>
#include "definitions.h"
#ifndef CSPARSECOMPILE
#include "matrix.h"
#endif
#ifdef CSPARSECOMPILE
#include "matrix2.h"
#endif
#include "bmat.h"
#include "bmat2d.h"
#include "hash.h"
#include "murmur3.h"
#include "datconfig2d.h"

void DatConfig2d_init(DatConfig2d * config, int32_t configLength, DatConfig * panmictic, SuperConfig * super)
{
	config->positions = (twoints *)malloc(sizeof(twoints) * (size_t)configLength);
	CHECKPOINTER(config->positions);
	config->panmictic = panmictic;
	config->super = super;
	config->prob = 0.0; 	// uninitialized value
	config->next = NULL; 	// uninitialized value
	return;
}

void DatConfig2d_print(DatConfig2d * config, FILE * output)
{
	int32_t i;
	fprintf(output, "------printing DatConfig2d------\n");
	fprintf(output, "prob = %f\n", config->prob);
	fprintf(output, "positions:\n");
	for(i = 0; i < config->super->configLength; i++)
		fprintf(output, "%i, %i\n", config->positions[i][0], config->positions[i][1]); 
	fprintf(output, "\n");
	return;
}

void DatConfig2d_get_ref_datconfig2d(BMat2d * bmat2d, DatConfig2d * config)
{
	int32_t i, j, maxJ, deme;
	for(i = 0; i < bmat2d->bmat.ncols+1; i++)
		config->positions[i][0] = config->positions[i][1] = 0;
	for(i = 0; i < bmat2d->bmat.nrows; i++)
	{
		deme = bmat2d->demes[i];
		maxJ = 0;
		for(j = 1; j <= bmat2d->bmat.ncols; j++)
		{
			if(bmat2d->bmat.mat[i][j-1] && j > maxJ)
				maxJ = j;
		}
		//printf("maxJ = %i, deme = %i\n", maxJ, deme);
		config->positions[maxJ][deme] += 1;
	}
	//for(i = 0; i < bmat2d->bmat.ncols+1; i++)
		//printf("config->positions[%i][0] = %i, config->positions[%i][1] = %i\n", i, config->positions[i][0], i, config->positions[i][1]);
	return;
}

void DatConfig2d_free(DatConfig2d * config)
{
	free(config->positions);
	return;
}

/* positions is an array for the panmictic configuration, and configLength is the length of the configuration */
//void SuperConfig_init(SuperConfig * super, int32_t * positions, int32_t * numChildren, int32_t configLength)

/*
typedef struct superconfig_
{
	DatConfig panmictic;
	DatConfig2d ** configs2d;
	//int32_t curNumConfigs2d; -- no longer necessary, calculating the index each time.
	int32_t numConfigs2d;
	int32_t configLength;
	int32_t * numChildren;
	int32_t * positionMultipliers;
	// array of int32_t's to multiply against, used to get the matrix /
	// DatConfig2d array index from a particular DatConfig2d's positions array.
	SuperEquations eq;
	struct superconfig_ * next;
} SuperConfig;
*/
void SuperConfig_init_static(SuperConfig * super, int32_t * numChildren, int32_t configLength, int32_t numThetas)
{
	super->configLength = configLength;
	DatConfig_init(&(super->panmictic), configLength, numChildren, numThetas);
	super->positionMultipliers = (int32_t *)calloc(configLength, sizeof(int32_t));
	CHECKPOINTER(super->positionMultipliers);
	return;
}

void SuperConfig_init(SuperConfig * super, DatConfig * panmictic, DataSet2d * ds)
{
	int32_t i, j, configLength = panmictic->length;
	int32_t * numChildren = panmictic->numChildren;
	int32_t * positions = panmictic->positions;
	for(i = 0; i < configLength; i++)
	{
		super->panmictic.positions[i] = positions[i];
		super->panmictic.active[i] = panmictic->active[i];
		for(j = 0; j < panmictic->numChildren[i]; j++)
			super->panmictic.satisfied[i][j] = panmictic->satisfied[i][j];
	}
	super->numChildren = numChildren;
	super->numConfigs2d = SuperConfig_get_num_configs2d(positions, configLength);
	super->configs2d = (DatConfig2d **)malloc(sizeof(DatConfig2d *) * (size_t)(super->numConfigs2d));
	CHECKPOINTER(super->configs2d);
	for(i = 0; i < super->numConfigs2d; i++)
	{
		super->configs2d[i] = (DatConfig2d *)malloc(sizeof(DatConfig2d));
		CHECKPOINTER(super->configs2d[i]);
		DatConfig2d_init(super->configs2d[i], super->configLength, &(super->panmictic), super);
	}
	super->next = NULL;
	SuperConfig_get_position_multipliers(positions, configLength, super->positionMultipliers);
	// matrix initialization
	SuperConfig_fill_out_matrix(super, ds->migRates, ds->theta);
	//SuperConfig_fill_in_datconfigs2d(super);
	SuperConfig_fill_in_datconfigs2d_nonrecursive(super);
	return;
}

/* migRates is an input parameter, a malloc'ed array of length 2 */
void SuperConfig_fill_out_matrix(SuperConfig * super, double * migRates, double theta)
{
	// nz initialized to numConfigs2d for the diagonal entries, which aren't counted later
	int32_t j, idx, n, n0, n1, nz = super->numConfigs2d;
	double totalRate;
	twoints * positions = (twoints *)calloc(super->configLength, sizeof(twoints));
	CHECKPOINTER(positions);
	twoints * positionsBucket = (twoints *)malloc(sizeof(twoints) * super->configLength);
	CHECKPOINTER(positionsBucket);
	for(idx = 0; idx < super->numConfigs2d; idx++)
	{
		SuperConfig_index_to_positions(idx, super, positions);
		nz += SuperConfig_calculate_num_row_entries_(positions, super->configLength);
	}
	// allocate matrix with nz...
	SuperEquations_init(&(super->eq), (int)(super->numConfigs2d), (int)(super->numConfigs2d), (int)nz);
	// fill out matrix
	for(idx = 0; idx < super->numConfigs2d; idx++)
	{
		SuperConfig_index_to_positions(idx, super, positions);
		n0 = 0;
		n1 = 0;
		for(j = 0; j < super->configLength; j++)
		{
			n0 += positions[j][0];
			n1 += positions[j][1];
		}
		n = n0 + n1;
		// n.b. all rates are multiplied by 2.
		totalRate = n0*(n0-1.0) + n1*(n1-1.0) + n*theta + n0*migRates[0] + n1*migRates[1];
		// fill out the matrix row...
		SuperConfig_fill_matrix_row(super, idx, positions, positionsBucket, migRates, totalRate);
	}
	free(positions);
	free(positionsBucket);
	return;
}

void SuperConfig_fill_matrix_row(SuperConfig * super, int32_t rowIdx, twoints * positions, twoints * positionsBucket, double * migRates, double totalRate)
{
	int32_t i, j, deme, colIdx;
	double migProb;
	for(deme = 0; deme < 2; deme++)
	{
		for(i = 0; i < super->configLength; i++)
		{
			if(positions[i][deme] == 0)
				continue;
			for(j = 0; j < super->configLength; j++)
			{
				positionsBucket[j][0] = positions[j][0];
				positionsBucket[j][1] = positions[j][1];
			}
			positionsBucket[i][deme]--;
			positionsBucket[i][!deme]++;
			colIdx = SuperConfig_get_index(positionsBucket, super->positionMultipliers, super->configLength);
			//printf("*** %i %i (deme = %i, migRates[deme] = %f, totalRate = %f) ***\n", positions[i][0], positions[i][1], deme, migRates[deme], totalRate);
			migProb = (double)(positions[i][deme]) * migRates[deme] / totalRate;
			SuperEquations_add_entry(&(super->eq), rowIdx, colIdx, -1.0 * migProb);
		}
	}
	SuperEquations_add_entry(&(super->eq), rowIdx, rowIdx, 1.0);
	return;
}

int32_t SuperConfig_calculate_num_row_entries_(twoints * positions, int32_t configLength)
{
	int32_t i, deme, numEntries = 0;
	for(i = 0; i < configLength; i++)
	{
		for(deme = 0; deme < 2; deme++)
		{
			if(positions[i][deme] > 0)
				numEntries++;
		}
	}
	return numEntries;
}



void SuperConfig_free(SuperConfig * config)
{
	int32_t i;
	for(i = 0; i < config->numConfigs2d; i++)
	{
		DatConfig2d_free(config->configs2d[i]);
		free(config->configs2d[i]);
	}
	free(config->configs2d);
	DatConfig_free(&(config->panmictic));
	free(config->positionMultipliers);
	/* TODO: free the matrices */
	return;
}

void SuperConfig_free_static(SuperConfig * config)
{
	DatConfig_free(&(config->panmictic));
	free(config->positionMultipliers);
	return;
}

/* this function frees only the arrays that change in size from SuperConfig to
 * SuperConfig within a dataset:
 *
 * numConfigs2d, configs2d
 * */
void SuperConfig_free_dynamic(SuperConfig * config)
{
	int32_t i;
	for(i = 0; i < config->numConfigs2d; i++)
	{
		DatConfig2d_free(config->configs2d[i]);
		free(config->configs2d[i]);
	}
	free(config->configs2d);
	config->next = NULL;
	SuperEquations_free(&(config->eq));
	return;
}
	
/* positions is the position-counts for the panmictic configuration of the SuperConfig
 * configLength is the length of that positions array. 
 * (See the fourth bullet of the left column of p. 615 of Wu's paper.) */
int32_t SuperConfig_get_num_configs2d(int32_t * positions, int32_t configLength)
{
	int32_t i, numConfigs = 1;
	for(i = 0; i < configLength; i++)
	{
		if(positions[i] > 0)
			numConfigs *= positions[i]+1;
	}
	return numConfigs;
}

int32_t SuperConfig_get_index(twoints * positions2d, int32_t * positionMultipliers, int32_t configLength)
{
	int32_t i, idx = 0;
	for(i = 0; i < configLength; i++)
	{
		if(positions2d[i][0] > 0)
			idx += positions2d[i][0] * positionMultipliers[i];
	}
	return idx;
}

void SuperConfig_index_to_positions(int32_t idx, SuperConfig * super, twoints * positions)
{
	int32_t i, remainingIdx = idx;
	for(i = super->configLength-1; i >= 0; i--)
	{
		if(super->positionMultipliers[i] > 0)
		{
			// n.b. integer division
			positions[i][0] = remainingIdx / super->positionMultipliers[i];
			positions[i][1] = super->panmictic.positions[i] - positions[i][0];
			remainingIdx %= super->positionMultipliers[i];
		}
		else
			positions[i][0] = positions[i][1] = 0;		// debug test.  this fixed some uninitialized-memory issues.
	}
	return;
}

//int32_t * SuperConfig_get_position_multipliers(int32_t * panmicticPositions, int32_t configLength)
// positionMultipliers is an input argument
void SuperConfig_get_position_multipliers(int32_t * panmicticPositions, int32_t configLength, int32_t * positionMultipliers)
{
	int32_t i, currentMultiplier = 1;
	CHECKPOINTER(positionMultipliers);
	for(i = 0; i < configLength; i++)
	{
		positionMultipliers[i] = 0;		// a big bug! wasn't previously resetting to 0.
		if(panmicticPositions[i] > 0)
		{
			positionMultipliers[i] = currentMultiplier;
			currentMultiplier *= panmicticPositions[i] + 1;
		}
	}
	return;
}

void SuperConfig_add_datconfig2d(SuperConfig * super, twoints * positions)
{
	int32_t i, idx;
	// test results seem to indicate that this is not necessary, but the proof requires more attention
	idx = SuperConfig_get_index(positions, super->positionMultipliers, super->configLength);
	DatConfig2d * config = super->configs2d[idx];
	for(i = 0; i < super->configLength; i++)
	{
		config->positions[i][0] = positions[i][0];
		config->positions[i][1] = positions[i][1];
	}
	/* everything else was taken care of when it was initialized. */
	return;
}

void SuperConfig_add_datconfig2d_nonrecursive(SuperConfig * super, twoints * positions, int32_t idx)
{
	int32_t i;
	DatConfig2d * config = super->configs2d[idx];
	for(i = 0; i < super->configLength; i++)
	{
		config->positions[i][0] = positions[i][0];
		config->positions[i][1] = positions[i][1];
	}
	/* everything else was taken care of when it was initialized. */
	return;
}

void SuperConfig_fill_in_datconfigs2d_recursive_(SuperConfig * super, twoints * positions, const int32_t * reference, int32_t len)
{
	int32_t i, j;
	// add it
	SuperConfig_add_datconfig2d(super, positions); 
	// get the next
	for(i = 0; i < len; i++)
	{
		if(positions[i][0] < reference[i])
		{
			positions[i][0]++;
			positions[i][1]--;
			for(j = 0; j < i; j++)
			{
				positions[j][0] = 0;
				positions[j][1] = reference[j];
			}
			SuperConfig_fill_in_datconfigs2d_recursive_(super, positions, reference, len);
		}
	}
	return;
}

/* TODO: make this use SuperConfig_index_to_positions for speed and non-recursion simplicity */
void SuperConfig_fill_in_datconfigs2d(SuperConfig * super)
{
	// reference is super->panmictic.positions
	int32_t i;
	const int32_t * reference = super->panmictic.positions;
	twoints * bucket = (twoints *)malloc(sizeof(twoints) * super->configLength);
	CHECKPOINTER(bucket);
	for(i = 0; i < super->configLength; i++)
	{
		bucket[i][0] = 0;
		bucket[i][1] = reference[i];
	}
	SuperConfig_fill_in_datconfigs2d_recursive_(super, bucket, reference, super->configLength);
	free(bucket);
	return;
}

void SuperConfig_fill_in_datconfigs2d_nonrecursive(SuperConfig * super)
{
	// reference is super->panmictic.positions
	int32_t i, j, numConfigs2d = super->numConfigs2d;
	twoints * bucket = (twoints *)malloc(sizeof(twoints) * super->configLength);
	CHECKPOINTER(bucket);
	for(i = 0; i < super->configLength; i++)		// debug initialization
		bucket[i][0] = bucket[i][1] = -1;
	for(i = 0; i < numConfigs2d; i++)
	{
		SuperConfig_index_to_positions(i, super, bucket);
		for(j = 0; j < super->configLength; j++)
			if(bucket[j][0] == -1 || bucket[j][1] == -1)
				printf("hi\n");
		SuperConfig_add_datconfig2d_nonrecursive(super, bucket, i);
	}
	free(bucket);
	return;
}

void SuperConfig_print(SuperConfig * config, FILE * output)
{
	int32_t i;
	fprintf(output, "--printing SuperConfig--\n");
	fprintf(output, "numConfigs2d = %i\n", config->numConfigs2d);
	fprintf(output, "configLength = %i\n", config->configLength);
	for(i = 0; i < config->configLength; i++)
		fprintf(output, "numChildren[%i] = %i\n", i, config->numChildren[i]);
	for(i = 0; i < config->configLength; i++)
		fprintf(output, "positionMultipliers[%i] = %i\n", i, config->positionMultipliers[i]);
	for(i = 0; i < config->numConfigs2d; i++)
		DatConfig2d_print(config->configs2d[i], output);
	return;
}

/* (end of SuperConfig functions) */

/////////////////////////
/* HashSuper functions */
/////////////////////////

void HashSuper_init(HashSuper * table, int32_t elementLength)
{
	int32_t i;
	table->elements = (SuperConfig **)malloc(sizeof(SuperConfig *) * DEFAULT_HASH_TABLE_SIZE);
	CHECKPOINTER(table->elements);
	for(i = 0; i < DEFAULT_HASH_TABLE_SIZE; i++)
		table->elements[i] = NULL;
	table->size = DEFAULT_HASH_TABLE_SIZE;
	table->numElements = 0;
	table->elementLength = elementLength;
	return;
}

void HashSuper_free(HashSuper * table)
{
	free(table->elements);
	return;
}

void HashSuper_reset(HashSuper * table)
{
	int32_t i;
	for(i = 0; i < table->size; i++)
		table->elements[i] = NULL;
	return;
}

void HashSuper_insert_config(SuperConfig * config, HashSuper * table)
{
	uint32_t hashIdx = HashSuper_get_hash_idx(config->panmictic.positions, table->elementLength);
	//printf("[HASH insertion into %p] hashIdx = %i, config = %p\n", table, hashIdx, config);
	SuperConfig * insertionPoint;
	config->next = NULL;
	if(table->elements[hashIdx] == NULL)
	{
		table->elements[hashIdx] = config;
		return;
	}
	//insertionPoint = config;  // a former bug
	insertionPoint = table->elements[hashIdx];
	//printf("[HASH using open addressing in %p]\n", table);
	while(insertionPoint->next != NULL)
		insertionPoint = insertionPoint->next;
	insertionPoint->next = config;
	return;
}


/* this function is only used when SuperConfigs are being *inserted* into the
 * table, not when it's being searched later on. the function computes the hash
 * and checks all SuperConfigs at that address. if a matching SuperConfig isn't
 * found at that address (either because the table entry is NULL or because none
 * of the open-addressed SuperConfigs match, then a pointer to the insertion
 * point is returned. If the SuperConfig *is* found, then the insertion point is
 * returned NULL. */
SuperConfig ** HashSuper_get_insertion_point(int32_t * positions, HashSuper * table)
{
	uint32_t hashIdx = HashSuper_get_hash_idx(positions, table->elementLength);
	SuperConfig * foundConfig = table->elements[hashIdx];
	if(foundConfig == NULL)
		return &(table->elements[hashIdx]);
	SuperConfig ** insertionPoint = &(table->elements[hashIdx]);
	while((*insertionPoint)->next != NULL)
	{
		if(HashSuper_configs_equal(positions, (*insertionPoint)->panmictic.positions, table->elementLength))
			return NULL;
		*insertionPoint = (*insertionPoint)->next;
	}
	if(HashSuper_configs_equal(positions, (*insertionPoint)->panmictic.positions, table->elementLength))
		return NULL;
	return insertionPoint;
}

int32_t HashSuper_check_in_table(int32_t * positions, HashSuper * table)
{
	uint32_t hashIdx = HashSuper_get_hash_idx(positions, table->elementLength);
	SuperConfig * foundConfig = table->elements[hashIdx];
	if(foundConfig == NULL)
		return 0;
	while(foundConfig->next != NULL)
	{
		if(HashSuper_configs_equal(positions, foundConfig->panmictic.positions, table->elementLength))
			return 1;
		foundConfig = foundConfig->next;
	}
	if(HashSuper_configs_equal(positions, foundConfig->panmictic.positions, table->elementLength))
		return 1;
	return 0;
}

SuperConfig * HashSuper_find_superconfig(int32_t * positions, HashSuper * table)
{
	uint32_t hashIdx = HashSuper_get_hash_idx(positions, table->elementLength);
	//fprintf(stdout, "hashIdx = %i\n", hashIdx);
	SuperConfig * foundConfig = table->elements[hashIdx];
	//printf("[HASH search of %p] hashIdx = %i, foundConfig = %p\n", table, hashIdx, foundConfig);
	// no elements in this bucket
	if(foundConfig == NULL)
		return NULL;
	// one element in this bucket
	if(foundConfig->next == NULL)
	{
		if(!HashSuper_configs_equal(positions, foundConfig->panmictic.positions, table->elementLength))
				printf("not equal!\n");
		return foundConfig;
	}
	// multiple elements in this bucket
	// (must check to see which is the correct one)
	do {
		if(HashSuper_configs_equal(foundConfig->panmictic.positions, positions, table->elementLength))
			return foundConfig;
		foundConfig = foundConfig->next;
	} while(foundConfig != NULL);
	PERROR("broken open addressing in HashSuper_find_config().");
}

int32_t HashSuper_configs_equal(int32_t * positions1, int32_t * positions2, int32_t configLength)
{
	int32_t i;
	for(i = 0; i < configLength; i++)
		if(positions1[i] != positions2[i])
			return 0;
	return 1;
}

uint32_t HashSuper_get_hash_idx_old(int32_t * positions, int32_t length)
{
	uint32_t hash = jenkins_one_at_a_time_hash((char *)positions, sizeof(int32_t) * length);
	hash %= DEFAULT_HASH_TABLE_SIZE;
	return hash;
}

uint32_t HashSuper_get_hash_idx(int32_t * positions, int32_t length)
{
	uint32_t hash128[4];
	uint32_t seed = 11;
	MurmurHash3_x64_128((const void *)positions, (int)(length * sizeof(int32_t)), seed, hash128);
	hash128[0] %= DEFAULT_HASH_TABLE_SIZE;
	return hash128[0];
}

/* (End HashSuper functions) */

///////////////////////////////
/* SuperCollection functions */
///////////////////////////////

void SuperCollection_init(SuperCollection * collection, int32_t configLength, int32_t * numChildren, int32_t numThetas)
{
	int32_t i;
	collection->numChildren = numChildren;
	collection->curNumSuperConfigs = 0;
	collection->maxNumSuperConfigs = DEFAULT_SUPER_COLLECTION_NUMSUPERCONFIGS;
	collection->configLength = configLength;
	collection->superConfigs = (SuperConfig **)malloc(sizeof(SuperConfig *) * DEFAULT_SUPER_COLLECTION_NUMSUPERCONFIGS);
	CHECKPOINTER(collection->superConfigs);
	for(i = 0; i < DEFAULT_SUPER_COLLECTION_NUMSUPERCONFIGS; i++)
	{
		collection->superConfigs[i] = (SuperConfig *)malloc(sizeof(SuperConfig));
		CHECKPOINTER(collection->superConfigs[i]);
		SuperConfig_init_static(collection->superConfigs[i], numChildren, configLength, numThetas);
	}
	HashSuper_init(&(collection->hashSuper), configLength);
	return;
}

void SuperCollection_reset(SuperCollection * collection)
{
	int32_t i;
	for(i = 0; i < collection->curNumSuperConfigs; i++)
		SuperConfig_free_dynamic(collection->superConfigs[i]);
	collection->curNumSuperConfigs = 0;
	return;
}

void SuperCollection_free(SuperCollection * collection)
{
	int32_t i;
	for(i = 0; i < collection->maxNumSuperConfigs; i++)
	{
		SuperConfig_free_static(collection->superConfigs[i]);
		free(collection->superConfigs[i]);
	}
	free(collection->superConfigs);
	HashSuper_free(&(collection->hashSuper));
	return;
}

void SuperCollection_add_SuperConfig_space(SuperCollection * collection, int32_t numNewSuperConfigs, int32_t numThetas)
{
	int32_t i, newNumSuperConfigs = collection->maxNumSuperConfigs + numNewSuperConfigs;
	collection->superConfigs = (SuperConfig **)realloc((void *)(collection->superConfigs), sizeof(SuperConfig *) * (size_t)newNumSuperConfigs);
	CHECKPOINTER(collection->superConfigs);
	for(i = collection->maxNumSuperConfigs; i < newNumSuperConfigs; i++)
	{
		collection->superConfigs[i] = (SuperConfig *)malloc(sizeof(SuperConfig));
		CHECKPOINTER(collection->superConfigs);
		// TODO: getting configLength from the first superConfigs is kind of sloppy.
		SuperConfig_init_static(collection->superConfigs[i], collection->numChildren, collection->superConfigs[0]->configLength, numThetas);
	}
	collection->maxNumSuperConfigs = newNumSuperConfigs;
	return;
}

SuperConfig * SuperCollection_get_SuperConfig(SuperCollection * collection, int32_t numThetas)
{
	if(collection->curNumSuperConfigs >= collection->maxNumSuperConfigs)
		SuperCollection_add_SuperConfig_space(collection, DEFAULT_SUPER_COLLECTION_INCREASE, numThetas);
	SuperConfig * config = collection->superConfigs[collection->curNumSuperConfigs];
	(collection->curNumSuperConfigs)++;
	return config;
}

/* this function is used to add a SuperConfig to a SuperCollection when there 
 * is a transfer of SuperACs between stages. In the previous stage, the SuperACs
 * of that stage are modified to make derivative SuperACs, which are added to 
 * this collection. They're cloned, so this one here needs to be cloned into a
 * SuperConfig withdrawn from collection. */
//void SuperCollection_add_SuperConfig(SuperCollection * collection, SuperConfig * cloneConfig)
void SuperCollection_add_SuperConfig(SuperCollection * collection, DatConfig * panmictic, struct dataset2d_ * ds)
{
	// this may not work...
	int32_t numThetas = panmictic->numProbs;
	//SuperConfig ** insertionPoint;
	// updated: now have to check whether it's in the collection already:
	if(!HashSuper_check_in_table(panmictic->positions, &(collection->hashSuper)))		// TODO: check/rewrite function to check and get insertion point (or NULL) at same time (see below).
	{
		SuperConfig * config = SuperCollection_get_SuperConfig(collection, numThetas);
		SuperConfig_init(config, panmictic, ds);
		HashSuper_insert_config(config, &(collection->hashSuper));
	}
	/*
	insertionPoint = HashSuper_get_insertion_point(panmictic->positions, &(collection->hashSuper));
	if(insertionPoint != NULL)
	{
		SuperConfig * config = SuperCollection_get_SuperConfig(collection);
		SuperConfig_init(config, panmictic, ds);
		*insertionPoint = config;
	}
	*/
	return;
}

void SuperCollection_print(SuperCollection * collection, FILE * output)
{
	int32_t i;
	fprintf(output, "----Printing SuperCollection----\n");
	fprintf(output, "curNumSuperConfigs = %i\n", collection->curNumSuperConfigs);
	fprintf(output, "maxNumSuperConfigs = %i\n", collection->maxNumSuperConfigs);
	fprintf(output, "configLength = %i\n", collection->configLength);
	for(i = 0; i < collection->configLength; i++)
		fprintf(output, "numChildren[%i] = %i\n", i, collection->numChildren[i]);
	fprintf(output, "(printing SuperConfigs)\n");
	for(i = 0; i < collection->curNumSuperConfigs; i++)
		SuperConfig_print(collection->superConfigs[i], output);
	return;
}

/* (end of SuperCollection functions) */

/* a test */
/*
int32_t main()
{
	SuperConfig super;
	SuperCollection collection;
	int32_t * positions = (int32_t *)malloc(sizeof(int32_t) * 3);
	CHECKPOINTER(positions);
	positions[0] = 2; positions[1] = 3; positions[2] = 2;
	int32_t * numChildren = (int32_t *)malloc(sizeof(int32_t) * 3);
	CHECKPOINTER(numChildren);
	numChildren[0] = 4; numChildren[1] = 1; numChildren[2] = 2;		// meaningless values
	// no longer up to date:
	SuperConfig_init(&super, positions, numChildren, 3);
	SuperCollection_init(&collection, 3, numChildren);
	SuperCollection_add_SuperConfig(&collection, &super);
	SuperCollection_print(&collection, stdout);
	SuperConfig_free(&super);
	SuperCollection_free(&collection);
	free(positions);
	free(numChildren);
	return 0;
}
*/
