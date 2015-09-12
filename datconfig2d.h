#ifndef DATCONFIG2D_HEADER
#define DATCONFIG2D_HEADER

#include "datconfig.h"
#include "definitions.h"
#include "hash.h"
#include "matrix.h"

#define DEFAULT_SUPER_COLLECTION_NUMSUPERCONFIGS 1000
#define DEFAULT_SUPER_COLLECTION_INCREASE 1000

typedef struct datconfig2d_ 
{
	twoints * positions;			// the counts of haplotypes in each of the two demes
	DatConfig * panmictic;	// the corresponding panmictic configuration
	double prob;					// probability of this data configuration
	struct superconfig_ * super; // pointer to the SuperAC to which this 2d config belongs.
	struct datconfig2d_ * next;		// for storing in the (2D) hash table.
	int numThetas;
	int numMigRates;
	double ** probs;
} DatConfig2d;

typedef struct superconfig_
{
	DatConfig panmictic;
	DatConfig2d ** configs2d;
	int numConfigs2d;
	int configLength;
	int * numChildren;
	int * positionMultipliers;
	int numThetas;
	int numMigRates;
	SuperEquations ** eqs;
	struct superconfig_ * next;
} SuperConfig;

typedef struct hashsuper_
{
	SuperConfig ** elements;
	int elementLength;
	int numElements;
	int size;
} HashSuper;

/* SuperCollection
 * 		A SuperCollection is a collection of SuperConfigs. There will be two
 * 		SuperCollections, one for the present stage and one for the next stage.
 * 		Each SuperCollection has a collection of SuperConfigs, representing the
 * 		SuperACs from Wu 2011. It is not generally known ahead of time how many
 * 		SuperConfigs there will be in a SuperCollection, so each is initialized
 * 		with an initial size and resized as necessary. 
 * 		
 * 		*/
typedef struct supercollection_
{
	SuperConfig ** superConfigs;
	// the collection of SuperConfigs, initialized to a size of
	// DEFAULT_SUPER_COLLECTION_NUMSUPERCONFIGS, incremented by
	// DEFAULT_SUPER_COLLECTION_INCREASE when necessary.
	int curNumSuperConfigs;
	// current number of SuperConfigs in the collection.
	// these are added by SuperCollection_add_SuperConfig, and this counter is
	// incremented when one is added. it is reset to zero when the collection is
	// reset (TODO).
	int maxNumSuperConfigs;
	// current number of slots allocated for SuperConfigs
	// initialized to DEFAULT_SUPER_COLLECTION_NUMSUPERCONFIGS, incremented by
	// DEFAULT_SUPER_COLLECTION_INCREASE when necessary.
	int configLength;			// panmictic config length? may not be necessary
	int * numChildren;			// panmictic numChildren for NodeList ... may not be necessary
	HashSuper hashSuper;		// for storing Super configurations in a given stage
	int numThetas;
	int numMigRates;
} SuperCollection;

struct dataset2d_;

#include "dataset2d.h"

void HashSuper_init(HashSuper * table, int elementLength);
void HashSuper_free(HashSuper * table);
void HashSuper_reset(HashSuper * table);
void HashSuper_insert_config(SuperConfig * config, HashSuper * table);
SuperConfig * HashSuper_find_superconfig(int * positions, HashSuper * table);
int HashSuper_configs_equal(int * positions1, int * positions2, int configLength);
uint32_t HashSuper_get_hash_idx_old(int * positions, int length);
uint32_t HashSuper_get_hash_idx(int * positions, int length);

void DatConfig2d_init(DatConfig2d * config, int configLength, DatConfig * panmictic, SuperConfig * super, int numThetas, int numMigRates);
void DatConfig2d_print(DatConfig2d * config, FILE * output);
void DatConfig2d_get_ref_datconfig2d(BMat2d * bmat2d, DatConfig2d * config);
void DatConfig2d_free(DatConfig2d * config);

void SuperConfig_init(SuperConfig * super, DatConfig * panmictic, struct dataset2d_ * ds);
void SuperConfig_init_static(SuperConfig * super, int * numChildren, int configLength, int numThetas, int numMigRates);
void SuperConfig_fill_out_matrices(SuperConfig * super, double * thetas, double * migRates);
void SuperConfig_fill_matrix_row(SuperConfig * super, SuperEquations * eq, int rowIdx, twoints * positions, twoints * positionsBucket, double migRate, double totalRate);
void SuperConfig_index_to_positions(int idx, SuperConfig * super, twoints * positions);
int SuperConfig_calculate_num_row_entries_(twoints * positions, int configLength);
void SuperConfig_free(SuperConfig * config);
void SuperConfig_free_dynamic(SuperConfig * config);
void SuperConfig_get_position_multipliers(int * panmicticPositions, int configLength, int * positionMultipliers);
int SuperConfig_get_num_configs2d(int * positions, int configLength);
int SuperConfig_get_index(twoints * positions2d, int * positionMultipliers, int configLength);
void SuperConfig_add_datconfig2d(SuperConfig * super, twoints * positions);
void SuperConfig_add_datconfig2d_nonrecursive(SuperConfig * super, twoints * positions, int idx);
void SuperConfig_fill_in_datconfigs2d_recursive_(SuperConfig * super, twoints * positions, const int * reference, int len);
void SuperConfig_fill_in_datconfigs2d(SuperConfig * super);
void SuperConfig_fill_in_datconfigs2d_nonrecursive(SuperConfig * super);
void SuperConfig_print(SuperConfig * config, FILE * output);

void SuperCollection_init(SuperCollection * collection, int configLength, int * numChildren, int numThetas, int numMigRates);
void SuperCollection_reset(SuperCollection * collection);
void SuperCollection_add_SuperConfig_space(SuperCollection * collection, int numNewSuperConfigs, int numThetas);
SuperConfig * SuperCollection_get_SuperConfig(SuperCollection * collection, int numThetas);
//void SuperCollection_add_SuperConfig(SuperCollection * collection, SuperConfig * cloneConfig);
void SuperCollection_add_SuperConfig(SuperCollection * collection, DatConfig * panmictic, struct dataset2d_ * ds);
void SuperCollection_free(SuperCollection * collection);
void SuperCollection_print(SuperCollection * collection, FILE * output);
#endif
