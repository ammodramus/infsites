#ifndef DATCONFIG2D_HEADER
#define DATCONFIG2D_HEADER

#include "datconfig.h"
#include "definitions.h"
#include "hash.h"
//#include "matrix2.h"  // CSPARSE EDIT
#ifndef CSPARSECOMPILE
#include "matrix.h"
#endif
#ifdef CSPARSECOMPILE
#include "matrix2.h"
#endif

#define DEFAULT_SUPER_COLLECTION_NUMSUPERCONFIGS 1000
#define DEFAULT_SUPER_COLLECTION_INCREASE 1000

typedef struct datconfig2d_ 
{
	twoints * positions;			// the counts of haplotypes in each of the two demes
	DatConfig * panmictic;	// the corresponding panmictic configuration
	double prob;					// probability of this data configuration
	struct superconfig_ * super; // pointer to the SuperAC to which this 2d config belongs.
	struct datconfig2d_ * next;		// for storing in the (2D) hash table.
	int32_t numThetas;
	int32_t numMigRates;
	double ** probs;
} DatConfig2d;

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
	int32_t numThetas;
	int32_t numMigRates;
	// for multiple thetas and migration rates
	SuperEquations ** eqs;
	struct superconfig_ * next;
} SuperConfig;

typedef struct hashsuper_
{
	SuperConfig ** elements;
	int32_t elementLength;
	int32_t numElements;
	int32_t size;
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
	int32_t curNumSuperConfigs;
	// current number of SuperConfigs in the collection.
	// these are added by SuperCollection_add_SuperConfig, and this counter is
	// incremented when one is added. it is reset to zero when the collection is
	// reset (TODO).
	int32_t maxNumSuperConfigs;
	// current number of slots allocated for SuperConfigs
	// initialized to DEFAULT_SUPER_COLLECTION_NUMSUPERCONFIGS, incremented by
	// DEFAULT_SUPER_COLLECTION_INCREASE when necessary.
	int32_t configLength;			// panmictic config length? may not be necessary
	int32_t * numChildren;			// panmictic numChildren for NodeList ... may not be necessary
	HashSuper hashSuper;		// for storing Super configurations in a given stage
	int32_t numThetas;
	int32_t numMigRates;
} SuperCollection;

struct dataset2d_;

#include "dataset2d.h"

void HashSuper_init(HashSuper * table, int32_t elementLength);
void HashSuper_free(HashSuper * table);
void HashSuper_reset(HashSuper * table);
void HashSuper_insert_config(SuperConfig * config, HashSuper * table);
SuperConfig * HashSuper_find_superconfig(int32_t * positions, HashSuper * table);
int32_t HashSuper_configs_equal(int32_t * positions1, int32_t * positions2, int32_t configLength);
uint32_t HashSuper_get_hash_idx_old(int32_t * positions, int32_t length);
uint32_t HashSuper_get_hash_idx(int32_t * positions, int32_t length);

void DatConfig2d_init(DatConfig2d * config, int32_t configLength, DatConfig * panmictic, SuperConfig * super, int32_t numThetas, int32_t numMigRates);
void DatConfig2d_print(DatConfig2d * config, FILE * output);
void DatConfig2d_get_ref_datconfig2d(BMat2d * bmat2d, DatConfig2d * config);
void DatConfig2d_free(DatConfig2d * config);

void SuperConfig_init(SuperConfig * super, DatConfig * panmictic, struct dataset2d_ * ds);
void SuperConfig_init_static(SuperConfig * super, int32_t * numChildren, int32_t configLength, int32_t numThetas, int32_t numMigRates);
void SuperConfig_fill_out_matrices(SuperConfig * super, double * thetas, double * migRates);
void SuperConfig_fill_matrix_row(SuperConfig * super, SuperEquations * eq, int32_t rowIdx, twoints * positions, twoints * positionsBucket, double migRate, double totalRate);
void SuperConfig_index_to_positions(int32_t idx, SuperConfig * super, twoints * positions);
int32_t SuperConfig_calculate_num_row_entries_(twoints * positions, int32_t configLength);
void SuperConfig_free(SuperConfig * config);
void SuperConfig_free_dynamic(SuperConfig * config);
void SuperConfig_get_position_multipliers(int32_t * panmicticPositions, int32_t configLength, int32_t * positionMultipliers);
int32_t SuperConfig_get_num_configs2d(int32_t * positions, int32_t configLength);
int32_t SuperConfig_get_index(twoints * positions2d, int32_t * positionMultipliers, int32_t configLength);
void SuperConfig_add_datconfig2d(SuperConfig * super, twoints * positions);
void SuperConfig_add_datconfig2d_nonrecursive(SuperConfig * super, twoints * positions, int32_t idx);
void SuperConfig_fill_in_datconfigs2d_recursive_(SuperConfig * super, twoints * positions, const int32_t * reference, int32_t len);
void SuperConfig_fill_in_datconfigs2d(SuperConfig * super);
void SuperConfig_fill_in_datconfigs2d_nonrecursive(SuperConfig * super);
void SuperConfig_print(SuperConfig * config, FILE * output);

void SuperCollection_init(SuperCollection * collection, int32_t configLength, int32_t * numChildren, int32_t numThetas, int32_t numMigRates);
void SuperCollection_reset(SuperCollection * collection);
void SuperCollection_add_SuperConfig_space(SuperCollection * collection, int32_t numNewSuperConfigs);
SuperConfig * SuperCollection_get_SuperConfig(SuperCollection * collection);
//void SuperCollection_add_SuperConfig(SuperCollection * collection, SuperConfig * cloneConfig);
void SuperCollection_add_SuperConfig(SuperCollection * collection, DatConfig * panmictic, struct dataset2d_ * ds);
void SuperCollection_free(SuperCollection * collection);
void SuperCollection_print(SuperCollection * collection, FILE * output);
#endif
