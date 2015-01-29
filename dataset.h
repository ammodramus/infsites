#ifndef DATASET_HEADER
#define DATASET_HEADER

#include <stdint.h>
#include "bmat.h"
#include "node.h"
#include "nodelist.h"
#include "datconfig.h"

typedef struct dataset_
{
	int32_t numSegSites;
	int32_t numSamples;
	BMat * bmat;
	BMat * Lij;
	int32_t * Lj;
	NodeList nodeList;
	DatConfig refConfig;
	ConfigCollection * collection[2];
	double migRate;		// later
	double theta;
} DataSet;

void DataSet_init(DataSet * ds, BMat * inputbmat, double theta);
void DataSet_init_print(DataSet * ds, BMat * inputbmat, double theta);
void DataSet_free(DataSet * ds);
void DataSet_transfer_config_collections(ConfigCollection * donor, ConfigCollection * recipient, DataSet * ds);
void DataSet_donate_deriv_configs(DatConfig * curConfig, ConfigCollection * recipient, Node ** idxToNode, DataSet * ds);
int32_t DataSet_binomial_coeff_(int32_t n, int32_t k);
int32_t DataSet_get_prob_multiplier(BMat * bmat);

#endif
