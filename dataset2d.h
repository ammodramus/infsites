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
	int numSegSites;
	int numSamples;
	int recipientCollection;
	BMat2d * bmat2d;
	BMat * Lij;
	int * Lj;
	NodeList nodeList;
	struct datconfig_ refConfig;
	struct datconfig2d_ refConfig2d;
	struct supercollection_ * collection[2];
	double migRatesPair[2];
	double theta;
	int numThetas;
	double * thetas;
	int numMigRates;
	double * migRates;
    int * initialNodes;
    int ordered;
} DataSet2d;


void DataSet2d_init(DataSet2d * ds, BMat2d * inputbmat, int numThetas, double * thetas, int numMigRates, double * migRates, int printAll, int ordered, int genetree);
void DataSet2d_solve_ctypes(DataSet2d * ds, BMat2d * inputbmat, int numThetas, double * thetas, int numMigRates, double * migRates, int ordered, int genetree, double * samplingProbs);
void DataSet2d_free(DataSet2d * ds);
void DataSet2d_iterate_stages(SuperCollection * donor, SuperCollection * recipient, DataSet2d * ds, int printAll);
void DataSet2d_solve_equations(SuperCollection * collection);
void DataSet2d_donate_deriv_superconfigs(SuperConfig * super, SuperCollection * recipient, DataSet2d * ds);
void DataSet2d_transfer_config_collections(SuperCollection * donor, SuperCollection * recipient, DataSet2d * ds);
void DataSet2d_link_supercollections(SuperCollection * donor, SuperCollection * recipient, DataSet2d * ds);
void DataSet2d_link_datconfig2ds(SuperConfig * donorConfig, SuperCollection * recipient, DataSet2d * ds);
void DataSet2d_link_probabilities(DatConfig2d * config, SuperCollection * recipient, DataSet2d * ds);
int DataSet2d_binomial_coeff_(int n, int k);
int DataSet2d_get_prob_multiplier(DataSet2d * ds);
int DataSet2d_get_prob_multiplier2(DatConfig2d * config);
void DataSet2d_donate_deriv_superconfigs(SuperConfig * super, SuperCollection * recipient, DataSet2d * ds);

#endif
