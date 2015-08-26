#ifndef DATASET_HEADER
#define DATASET_HEADER

#include <stdint.h>
#include "bmat.h"
#include "node.h"
#include "nodelist.h"
#include "datconfig.h"

typedef struct dataset_
{
	int numSegSites;
	int numSamples;
	BMat * bmat;
	BMat * Lij;
	int * Lj;
	NodeList nodeList;
	DatConfig refConfig;
	ConfigCollection * collection[2];
	double migRate;		// later
	double theta;
	int numThetas;
	double * thetas;
    int * initialNodes;
    double probMultiplier;
    int printAll;
    int ordered;
} DataSet;

void DataSet_init(DataSet * ds, BMat * inputbmat, int numThetas, double * thetas, int printAll, int ordered);
void DataSet_free(DataSet * ds);
void DataSet_transfer_config_collections(ConfigCollection * donor, ConfigCollection * recipient, DataSet * ds);
void DataSet_donate_deriv_configs(DatConfig * curConfig, ConfigCollection * recipient, Node ** idxToNode, DataSet * ds);
int DataSet_binomial_coeff_(int n, int k);
int DataSet_get_prob_multiplier(BMat * bmat);
int DataSet_get_prob_multiplier2(DatConfig * refconfig);

#endif
