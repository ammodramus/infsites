#ifndef DATCONFIG_HEADER
#define DATCONFIG_HEADER

#include <stdio.h>
#include <stdint.h>
#include "bmat.h"
#include "nodelist.h"

typedef struct datconfig_ 
{
	int32_t length;
	int32_t * positions;
	int32_t ** satisfied;
	int32_t * numChildren;		// to be shared by all DatConfigs in the same DataSet
	int32_t * active;
	double prob;
	struct datconfig_ * next;
} DatConfig;


typedef struct hashtable_
{
	DatConfig ** elements;
	int32_t elementLength;
	int32_t numElements;
	int32_t size;
} HashTable;

typedef struct configcollection_
{
	DatConfig ** configs;
	int32_t curNumConfigs;
	int32_t maxNumConfigs;
	int32_t configLength;
	int32_t * numChildren;
	HashTable hashTable;
} ConfigCollection;



void DatConfig_init(DatConfig * df, int32_t length, int32_t * numChildren);
void DatConfig_free(DatConfig * df);
void DatConfig_get_ref_config(BMat * bmat, DatConfig * config);
void DatConfig_print(DatConfig * df, FILE * output, int32_t tabCount);
void DatConfig_copy_dimensions(DatConfig * df, DatConfig * dfClone);
void DatConfig_set_root_config(DatConfig * df);
void DatConfig_copy_config(DatConfig * df, DatConfig * dfClone);

void HashTable_init(HashTable * table, int32_t elementLength);
void HashTable_free(HashTable * table);
void HashTable_reset(HashTable * table);
void HashTable_insert_config(DatConfig * config, HashTable * table);
DatConfig * HashTable_find_config(DatConfig * config, HashTable * table);
int32_t HashTable_configs_equal(DatConfig * config1, DatConfig * config2);
uint32_t HashTable_get_hash_idx(int32_t * positions, int32_t length);

#define DEFAULT_CONFIGCOLLECTION_NUMCONFIGS 1000
#define DEFAULT_CONFIGCOLLECTION_INCREASESIZE 1000

void ConfigCollection_init(ConfigCollection * cc, int32_t configLength, int32_t * numChildren);
void ConfigCollection_free(ConfigCollection * cc);
void ConfigCollection_add_config_space(ConfigCollection * cc, int32_t numNewConfigs);
DatConfig * ConfigCollection_get_empty_config(ConfigCollection * cc);
void ConfigCollection_add_config(ConfigCollection * cc, DatConfig * config);
void ConfigCollection_reset(ConfigCollection * cc);
void ConfigCollection_print(ConfigCollection * cc, FILE * output);
double ConfigCollection_get_final_prob(ConfigCollection * cc);

#endif
