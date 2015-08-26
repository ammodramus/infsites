#ifndef DATCONFIG_HEADER
#define DATCONFIG_HEADER

#include <stdio.h>
#include <stdint.h>
#include "bmat.h"
#include "nodelist.h"

typedef struct datconfig_ 
{
	int length;
	int * positions;
	int ** satisfied;
	int * numChildren;		// to be shared by all DatConfigs in the same DataSet
	int * active;
	double prob;
	/* <experimental> */
	int numProbs;
	double * probs;
	/* </experimental> */
	struct datconfig_ * next;
} DatConfig;


typedef struct hashtable_
{
	DatConfig ** elements;
	int elementLength;
	int numElements;
	int size;
} HashTable;

typedef struct configcollection_
{
	DatConfig ** configs;
	int curNumConfigs;
	int maxNumConfigs;
	int configLength;
	int * numChildren;
	HashTable hashTable;
	int numThetas;
} ConfigCollection;



void DatConfig_init(DatConfig * df, int length, int * numChildren, int numProbs);
void DatConfig_free(DatConfig * df);
void DatConfig_get_ref_config(BMat * bmat, DatConfig * config);
void DatConfig_print(DatConfig * df, FILE * output, int tabCount);
void DatConfig_copy_dimensions(DatConfig * df, DatConfig * dfClone);
void DatConfig_set_root_config(DatConfig * df);
void DatConfig_copy_config(DatConfig * df, DatConfig * dfClone);
void DatConfig_set_initial_node_indices(DatConfig * config, int * initialNodeIndices);

void HashTable_init(HashTable * table, int elementLength);
void HashTable_free(HashTable * table);
void HashTable_reset(HashTable * table);
void HashTable_insert_config(DatConfig * config, HashTable * table);
DatConfig * HashTable_find_config(DatConfig * config, HashTable * table);
int HashTable_configs_equal(DatConfig * config1, DatConfig * config2);
uint32_t HashTable_get_hash_idx(int * positions, int length);

#define DEFAULT_CONFIGCOLLECTION_NUMCONFIGS 1000
#define DEFAULT_CONFIGCOLLECTION_INCREASESIZE 1000

void ConfigCollection_init(ConfigCollection * cc, int configLength, int * numChildren, int numThetas);
void ConfigCollection_free(ConfigCollection * cc);
void ConfigCollection_add_config_space(ConfigCollection * cc, int numNewConfigs);
DatConfig * ConfigCollection_get_empty_config(ConfigCollection * cc);
void ConfigCollection_add_config(ConfigCollection * cc, DatConfig * config);
void ConfigCollection_reset(ConfigCollection * cc);
void ConfigCollection_print(ConfigCollection * cc, FILE * output);
double ConfigCollection_get_final_prob(ConfigCollection * cc, int thetaIdx);

#endif
