#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include "definitions.h"
#include "bmat.h"
#include "datconfig.h"
#include "nodelist.h"
#include "dataset.h"
#include "hash.h"

void DatConfig_init(DatConfig * df, int32_t length, int32_t * numChildren, int32_t numProbs)
{
	int32_t i;
	df->length = length;
	df->numChildren = numChildren;
	df->positions = (int32_t *)calloc((size_t)length, sizeof(int32_t));
	CHECKPOINTER(df->positions);
	df->satisfied = (int32_t **)malloc(sizeof(int32_t *) * (size_t)(length));
	CHECKPOINTER(df->satisfied);
	for(i = 0; i < length; i++)
	{
		df->satisfied[i] = (int32_t *)calloc((size_t)(numChildren[i]), sizeof(int32_t));
		CHECKPOINTER(df->satisfied[i]);
	}
	df->active = (int32_t *)calloc((size_t)length, sizeof(int32_t));
	CHECKPOINTER(df->active);
	df->numProbs = numProbs;
	df->probs = (double *)malloc(sizeof(double) * numProbs);
	return;
}

void DatConfig_free(DatConfig * df)
{
	int32_t i;
	for(i = 0; i < df->length; i++)
		free(df->satisfied[i]);
	free(df->satisfied);
	free(df->positions);
	free(df->active);
	free(df->probs);
	return;
}


void DatConfig_get_ref_config(BMat * bmat, DatConfig * config)
{
	int32_t i, j, maxJ;
	//for(i = 0; i < bmat->ncols; i++)
	for(i = 0; i < bmat->ncols+1; i++)		// if single-deme stops working, could be here.
		config->positions[i] = config->active[i] = 0;
	for(i = 0; i < bmat->nrows; i++)
	{
		maxJ = 0;
		for(j = 1; j <= bmat->ncols; j++)
		{
			if(bmat->mat[i][j-1] && j > maxJ)
				maxJ = j;
		}
		config->positions[maxJ] += 1;
	}
	for(i = 0; i < config->length; i++)
	{
		for(j = 0; j < config->numChildren[i]; j++)
			config->satisfied[i][j] = 0;
	}
	return;
}

void DatConfig_set_root_config(DatConfig * df)
{
	int32_t i, j;
	for(i = 0; i < df->length; i++)
		df->positions[i] = df->active[i] = 0;
	df->positions[0] = df->active[0] = 1;
	for(i = 0; i < df->length; i++)
	{
		for(j = 0; j < df->numChildren[i]; j++)
			df->satisfied[i][j] = 0;
	}
	for(i = 0; i < df->numProbs; i++)
		df->probs[i] = 1.0;
	return;
}

void DatConfig_print(DatConfig * df, FILE * output, int32_t tabCount)
{
	int32_t i,j;
	for(i = 0; i < tabCount; i++)
		fprintf(output, "\t");
	fprintf(output, "\n----Printing DatConfig---- (length: %i, prob %f)\n", df->length, df->prob);
	for(i = 0; i < tabCount; i++)
		fprintf(output, "\t");
	fprintf(output, "numChildren: ");
	for(i = 0; i < df->length; i++)
		fprintf(output, "%i ", df->numChildren[i]);
	fprintf(output, "\n");
	for(i = 0; i < tabCount; i++)
		fprintf(output, "\t");
	fprintf(output,"positions (satisfied): ");
	for(i = 0; i < df->length; i++)
	{
		fprintf(output, "%i ( ", df->positions[i]);
		for(j = 0; j < df->numChildren[i]; j++)
			fprintf(output, "%i ", df->satisfied[i][j]);
		fprintf(output, ") ");
	}
	fprintf(output, "\n");
	for(i = 0; i < tabCount; i++)
		fprintf(output, "\t");
	fprintf(output, "active: ");
	for(i = 0; i < df->length; i++)
		fprintf(output, "%i ", df->active[i]);
	fprintf(output, "\n\n");
	return;
}

/* copies the dimensions of df into dfClone */
void DatConfig_copy_dimensions(DatConfig * df, DatConfig * dfClone)
{
	DatConfig_init(dfClone, df->length, df->numChildren, df->numProbs);
	return;
}

/* this function copies the contents of df in to dfClone. it assumes that 
 * dfClone's vectors (positions, satisfied, numChildren, active) are the same
 * dimensions as df's (using DatConfig_copy_dimensions()) */
void DatConfig_copy_config(DatConfig * df, DatConfig * dfClone)
{
	int32_t i,j;
	dfClone->length = df->length;
	dfClone->numChildren = df->numChildren;
	dfClone->prob = df->prob;
	for(i = 0; i < df->length; i++)
	{
		dfClone->positions[i] = df->positions[i];
		dfClone->active[i] = df->active[i];
		for(j = 0; j < df->numChildren[i]; j++)
			dfClone->satisfied[i][j] = df->satisfied[i][j];
	}
	return;
}

/////////////////////////
/* HashTable functions */
/////////////////////////

void HashTable_init(HashTable * table, int32_t elementLength)
{
	int32_t i;
	table->elements = (DatConfig **)malloc(sizeof(DatConfig *) * DEFAULT_HASH_TABLE_SIZE);
	CHECKPOINTER(table->elements);
	for(i = 0; i < DEFAULT_HASH_TABLE_SIZE; i++)
		table->elements[i] = NULL;
	table->numElements = 0;
	table->elementLength = elementLength;
	table->size = DEFAULT_HASH_TABLE_SIZE;
	return;
}

void HashTable_free(HashTable * table)
{
	free(table->elements);
	return;
}

void HashTable_reset(HashTable * table)
{
	int32_t i;
	//for(i = 0; i < table->numElements; i++)    // a former bug (14 Jan 2014)
	for(i = 0; i < table->size; i++)
		table->elements[i] = NULL;
	return;
}

void HashTable_insert_config(DatConfig * config, HashTable * table)
{
	int32_t hashIdx = HashTable_get_hash_idx(config->positions, table->elementLength);
	DatConfig * insertionPoint;
	config->next = NULL;
	if(table->elements[hashIdx] == NULL)
	{
		table->elements[hashIdx] = config;
		return;
	}
	// insertionPoint = config; // a former bug
	insertionPoint = table->elements[hashIdx];
	while(insertionPoint->next != NULL)
		insertionPoint = insertionPoint->next;
	insertionPoint->next = config;
	return;
}

DatConfig * HashTable_find_config(DatConfig * config, HashTable * table)
{
	uint32_t hashIdx = HashTable_get_hash_idx(config->positions, table->elementLength);
	//fprintf(stdout, "hashIdx = %i\n", hashIdx);
	DatConfig * foundConfig;
	foundConfig = table->elements[hashIdx];
	// no elements in this bucket
	if(foundConfig == NULL)
		return NULL;
	// one element in this bucket
	//if(foundConfig->next == NULL)
		//return foundConfig;    // a bug: we don't know whether foundConfig is the same as config.
	// multiple elements in this bucket
	// (must check to see which is the correct one)
	do {
		//printf("open addressing...\n");
		if(HashTable_configs_equal(foundConfig, config))
			return foundConfig;
		foundConfig = foundConfig->next;
	} while(foundConfig != NULL);
	//printf("would have been buggy here.\n");
	return NULL;
	PERROR("broken open addressing in HashTable_find_config().");
}

int32_t HashTable_configs_equal(DatConfig * config1, DatConfig * config2)
{
	int32_t i;
	if(config1->length != config2->length)
		PERROR("comparing two configs of different lengths");
	for(i = 0; i < config1->length; i++)
		if(config1->positions[i] != config2->positions[i])
			return 0;
	return 1;
}

uint32_t HashTable_get_hash_idx(int32_t * positions, int32_t length)
{
	uint32_t hash = jenkins_one_at_a_time_hash((char *)positions, sizeof(int32_t) * length);
	hash %= DEFAULT_HASH_TABLE_SIZE;
	return hash;
}

/* (end of HashTable functions) */


////////////////////////////////
/* ConfigCollection functions */
////////////////////////////////

void ConfigCollection_init(ConfigCollection * cc, int32_t configLength, int32_t * numChildren, int32_t numThetas)
{
	int32_t i;
	cc->configs = (DatConfig **)malloc(sizeof(DatConfig *) * (size_t)(DEFAULT_CONFIGCOLLECTION_NUMCONFIGS));
	CHECKPOINTER(cc->configs);
	for(i = 0; i < DEFAULT_CONFIGCOLLECTION_NUMCONFIGS; i++)
	{
		cc->configs[i] = (DatConfig *)malloc(sizeof(DatConfig));
		CHECKPOINTER(cc->configs[i]);
		DatConfig_init(cc->configs[i], configLength, numChildren, numThetas);
	}
	cc->curNumConfigs = 0;
	cc->maxNumConfigs = DEFAULT_CONFIGCOLLECTION_NUMCONFIGS;
	cc->configLength = configLength;
	cc->numChildren = numChildren;
	cc->numThetas = numThetas;
	HashTable_init(&(cc->hashTable), configLength);
	return;
}

void ConfigCollection_free(ConfigCollection * cc)
{
	int32_t i;
	for(i = 0; i < cc->maxNumConfigs; i++)
	{
		DatConfig_free(cc->configs[i]);
		free(cc->configs[i]);
	}
	free(cc->configs);
	HashTable_free(&(cc->hashTable));
	return;
}

void ConfigCollection_add_config_space(ConfigCollection * cc, int32_t numNewConfigs)
{
	int32_t i, newNumConfigs = cc->maxNumConfigs + numNewConfigs;
	cc->configs = (DatConfig **)realloc((void *)cc->configs, sizeof(DatConfig *) * (size_t)newNumConfigs);
	CHECKPOINTER(cc->configs);
	for(i = cc->maxNumConfigs; i < newNumConfigs; i++)
	{
		cc->configs[i] = (DatConfig *)malloc(sizeof(DatConfig));
		CHECKPOINTER(cc->configs[i]);
		DatConfig_init(cc->configs[i], cc->configLength, cc->numChildren, cc->numThetas);
	}
	cc->maxNumConfigs = newNumConfigs;
	return;
}

DatConfig * ConfigCollection_get_empty_config(ConfigCollection * cc)
{
	DatConfig * nextConfig;
	if(cc->curNumConfigs+1 >= cc->maxNumConfigs)
		ConfigCollection_add_config_space(cc, DEFAULT_CONFIGCOLLECTION_INCREASESIZE);
	nextConfig = cc->configs[cc->curNumConfigs];
	nextConfig->next = NULL;		// reset the HashTable open addressing.
	cc->curNumConfigs++;
	return nextConfig;
}

void ConfigCollection_add_config(ConfigCollection * cc, DatConfig * config)
{
	DatConfig * hashConfig;
	DatConfig * cloneConfig;
	// look for config
	hashConfig = HashTable_find_config(config, &(cc->hashTable));
	if(hashConfig == NULL)
	{
		// first, clone the config into a config we are storing...
		cloneConfig = ConfigCollection_get_empty_config(cc);
		DatConfig_copy_config(config, cloneConfig);
		// then insert the config we're storing into the hash table.
		HashTable_insert_config(cloneConfig, &(cc->hashTable));
	}
	else
	{
		//fprintf(stdout, "(adding prob...) \"old prob\" = %f, config->prob = %f\n", hashConfig->prob, config->prob);
		//DatConfig_print(hashConfig, stdout, 1);
		hashConfig->prob += config->prob;
	}
	return;
}

/* For now, just reset the curNumConfigs index after each resetting of a
 * collection, between stages. */
void ConfigCollection_reset(ConfigCollection * cc)
{
	cc->curNumConfigs = 0;
	return;
}

void ConfigCollection_print(ConfigCollection * cc, FILE * output)
{
	int32_t i;
	fprintf(output,"\n====Printing ConfigCollection====\n\n");
	fprintf(output, "curNumConfigs = %i\n", cc->curNumConfigs);
	for(i = 0; i < cc->curNumConfigs; i++)
	{
		fprintf(output, "config[%i]:\n", i);
		DatConfig_print(cc->configs[i], stdout, 1);
	}
	fprintf(output,"=================================\n\n");
	return;
}

double ConfigCollection_get_final_prob(ConfigCollection * cc, int32_t thetaIdx)
{
	if(cc->curNumConfigs != 1)
	{
		ConfigCollection_print(cc, stdout);
		PERROR("Attempting to print zero or multiple probabilities in ConfigCollection_print_prob.");
	}
	return cc->configs[0]->probs[thetaIdx];
}

/* (end of ConfigCollection functions) */

/*
int	main()
{
	int32_t i;
	ConfigCollection cc;
	int32_t * numChildren = (int32_t *)malloc(sizeof(int32_t) * 6);
	numChildren[0] = 1;
	numChildren[1] = 1;
	numChildren[2] = 4;
	numChildren[3] = 5;
	numChildren[4] = 3;
	numChildren[5] = 4;
	ConfigCollection_init(&cc, 6, numChildren);
	DatConfig config1, config2, config3;
	DatConfig_init(&config1, 6, numChildren);
	DatConfig_init(&config2, 6, numChildren);
	DatConfig_init(&config3, 6, numChildren);
	for(i = 0; i < 6; i++)
	{
		printf("%i ", config1.numChildren[i]);
		config1.positions[i] = config3.positions[i] = i;
		config2.positions[i] = i+1;
	}
	printf("\n");
	ConfigCollection_add_config(&cc, &config1);
	ConfigCollection_add_config(&cc, &config2);
	ConfigCollection_add_config(&cc, &config3);
	printf("\n");
	return;
}
*/
