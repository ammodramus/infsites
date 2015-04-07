#ifndef DATASET2D_H
#define DATASET2D_H

#include <stdint.h>
#include "bmat.h"
#include "bmat2d.h"
#include "nodelist.h"
#include "datconfig2d.h"

//struct datconfig_;
//struct datconfig2d_;
//struct supercollection_;

typedef struct dataset2d_
{
	int32_t numSegSites;
	int32_t numSamples;
	int32_t recipientCollection;
	BMat2d * bmat2d;
	BMat * Lij;
	int32_t * Lj;
	NodeList nodeList;
	struct datconfig_ refConfig;
	struct datconfig2d_ refConfig2d;
	struct supercollection_ * collection[2];
	double migRatesPair[2];
	double theta;
	int32_t numThetas;
	double * thetas;
	int32_t numMigRates;
	double * migRates;
    int32_t * initialNodes;
} DataSet2d;


void DataSet2d_init(DataSet2d * ds, BMat2d * inputbmat, int32_t numThetas, double * thetas, int32_t numMigRates, double * migRates);
void DataSet2d_free(DataSet2d * ds);
void DataSet2d_iterate_stages(SuperCollection * donor, SuperCollection * recipient, DataSet2d * ds);
void DataSet2d_solve_equations(SuperCollection * collection);
void DataSet2d_donate_deriv_superconfigs(SuperConfig * super, SuperCollection * recipient, DataSet2d * ds);
void DataSet2d_transfer_config_collections(SuperCollection * donor, SuperCollection * recipient, DataSet2d * ds);
void DataSet2d_link_supercollections(SuperCollection * donor, SuperCollection * recipient, DataSet2d * ds);
void DataSet2d_link_datconfig2ds(SuperConfig * donorConfig, SuperCollection * recipient, DataSet2d * ds);
void DataSet2d_link_probabilities(DatConfig2d * config, SuperCollection * recipient, DataSet2d * ds);
int32_t DataSet2d_binomial_coeff_(int32_t n, int32_t k);
int32_t DataSet2d_get_prob_multiplier(DataSet2d * ds);
void DataSet2d_donate_deriv_superconfigs(SuperConfig * super, SuperCollection * recipient, DataSet2d * ds);

#endif
