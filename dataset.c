#include <stdlib.h>
#include <stdint.h>
#include "definitions.h"
#include "bmat.h"
#include "node.h"
#include "nodelist.h"
#include "datconfig.h"
#include "dataset.h"
#include "hash.h"

//void DataSet_init(DataSet * ds, BMat * inputbmat, double theta)
void DataSet_init(DataSet * ds, BMat * inputbmat, int32_t numThetas, double * thetas)
{
	int32_t i,j,k, zeroFirst, numStages, probMultiplier;
	double finalProb;
	ds->bmat = inputbmat;
	ds->numSegSites = inputbmat->ncols;
	ds->numSamples = inputbmat->nrows;
	numStages = ds->numSegSites + ds->numSamples;
	//ds->theta = theta;
	ds->numThetas = numThetas;
	ds->thetas = thetas;
	ds->collection[0] = (ConfigCollection *)malloc(sizeof(ConfigCollection));
	CHECKPOINTER(ds->collection[0]);
	ds->collection[1] = (ConfigCollection *)malloc(sizeof(ConfigCollection));
	CHECKPOINTER(ds->collection[1]);
	/* there are one more nodes in the phylogeny than
	 * segregating sites. */
	//printf("\nBefore ordering columns\n");
	//BMat_print(ds->bmat, stdout);
	BMat_order_columns(ds->bmat);
	//printf("\nAfter ordering columns\n");
	//BMat_print(ds->bmat, stdout);
	probMultiplier = DataSet_get_prob_multiplier(ds->bmat);
	ds->Lij = (BMat *)malloc(sizeof(BMat));
	CHECKPOINTER(ds->Lij);
	BMat_init(ds->Lij, ds->bmat->nrows, ds->bmat->ncols);
	for(i = 0; i < ds->bmat->nrows; i++)
		for(j = 0; j < ds->bmat->ncols; j++)
			ds->Lij->mat[i][j] = BMat_get_Lij(ds->bmat, i, j);
	ds->Lj = (int32_t *)malloc(sizeof(int32_t) * (size_t)ds->bmat->ncols);
	CHECKPOINTER(ds->Lj);
	for(j = 0; j < ds->bmat->ncols; j++)
		ds->Lj[j] = BMat_get_Lj(ds->bmat, ds->Lij, j);
	if(!BMat_determine_good(ds->bmat, ds->Lij, ds->Lj))
		PERROR("Data does not conform to infinite-sites mutation model.");
	NodeList_create_phylogeny(&(ds->nodeList), ds->bmat, ds->Lj);
	//NodeList_print(&(ds->nodeList), stdout);
	NodeList_get_num_children(&(ds->nodeList));
	NodeList_get_idxToNode(&(ds->nodeList));
	DatConfig_init(&(ds->refConfig), ds->bmat->ncols+1, ds->nodeList.numChildren, ds->numThetas);
	DatConfig_get_ref_config(ds->bmat, &(ds->refConfig));
	DatConfig rootConfig;
	DatConfig_init(&rootConfig, ds->bmat->ncols+1, ds->nodeList.numChildren, ds->numThetas);
	DatConfig_set_root_config(&rootConfig);
	/* initialize the ConfigCollections, and then add the root configuration to
	 * the first Collection, and off we go! */
	ConfigCollection_init(ds->collection[0], ds->bmat->ncols+1, ds->nodeList.numChildren, ds->numThetas);
	ConfigCollection_init(ds->collection[1], ds->bmat->ncols+1, ds->nodeList.numChildren, ds->numThetas);
	ConfigCollection_add_config(ds->collection[0], &rootConfig);
	zeroFirst = 0;
	for(i = 0; i < numStages-1; i++)
	{
		DataSet_transfer_config_collections(ds->collection[zeroFirst], ds->collection[!zeroFirst], ds);
		//printf("stage = %i\n", i);
		//fprintf(stdout,"ds->collection[%i]->curNumConfigs = %i\n", zeroFirst, ds->collection[zeroFirst]->curNumConfigs);
		//fprintf(stdout,"ds->collection[%i]->curNumConfigs = %i\n", !zeroFirst, ds->collection[!zeroFirst]->curNumConfigs);
		zeroFirst = !zeroFirst;
	}
	for(k = 0; k < ds->numThetas; k++)
	{
		finalProb = ConfigCollection_get_final_prob(ds->collection[zeroFirst], k) * (double)probMultiplier;
		fprintf(stdout, "%e\n", finalProb);
	}

	// freeing memory
	ConfigCollection_free(ds->collection[0]);
	ConfigCollection_free(ds->collection[1]);
	free(ds->collection[0]);
	free(ds->collection[1]);
	DatConfig_free(&(ds->refConfig));
	DatConfig_free(&rootConfig);
	//NodeList_free(&(ds->nodeList));
	return;
}

void DataSet_free(DataSet * ds)
{
	BMat_free(ds->Lij);
	NodeList_free(&(ds->nodeList));
	free(ds->Lij);
	free(ds->Lj);
	return;
}

void DataSet_transfer_config_collections(ConfigCollection * donor, ConfigCollection * recipient, DataSet * ds)
{
	int32_t i;
	ConfigCollection_reset(recipient);
	for(i = 0; i < donor->curNumConfigs; i++)
		DataSet_donate_deriv_configs(donor->configs[i], recipient, ds->nodeList.idxToNode, ds);
	HashTable_reset(&(donor->hashTable));
	return;
}

void DataSet_donate_deriv_configs(DatConfig * curConfig, ConfigCollection * recipient, Node ** idxToNode, DataSet * ds)
{
	int32_t i, j, k, n, ntot, ntotCoal, nref, nunsat, numChildren, mutIdx, lingering;
	double transitionProb;
	// TODO: write separate functions for copying dimensions and copying
	// contents of a DatConfig
	DatConfig deriv;
	DatConfig_copy_dimensions(curConfig, &deriv);
	ntot = 0;
	for(j = 0; j < curConfig->length; j++)
		ntot += curConfig->positions[j];
	ntotCoal = ntot + 1;
	for(i = 0; i < curConfig->length; i++)
	{
		numChildren = curConfig->numChildren[i];
		if(curConfig->active[i])
		{
			//////////////
			// COALESCENCE
			nunsat = 0;
			n = curConfig->positions[i];
			for(j = 0; j < numChildren; j++)
				if(curConfig->satisfied[i][j] == 0)
					nunsat++;
			nref = ds->refConfig.positions[i];
			//fprintf(stdout, "nunsat = %i, nref = %i, n = %i\n", nunsat, nref, n);
			// coalescence conditions
			if(n < nref + nunsat)
			{
				DatConfig_copy_config(curConfig, &deriv);
				deriv.positions[i] += 1;
				for(k = 0; k < ds->numThetas; k++)
				{
					REPORTF(ds->thetas[k]);
					transitionProb =  (double)((n+1.0)*n)/(double)((ntotCoal)*ds->thetas[k] + ntotCoal*(ntotCoal-1.0));
					REPORTF(transitionProb);
					deriv.probs[k] = curConfig->probs[k] * transitionProb;
				}
				printf("\n");
				//fprintf(stdout, "Coalescing : curConfig->prob = %f, transitionProb = %f\n", curConfig->prob, transitionProb);
				if(nunsat == 0 && nref == deriv.positions[i]) 	// check for completion
					deriv.active[i] = 0;
				ConfigCollection_add_config(recipient, &deriv);		// this function adds the probability.
			}
			///////////
			// MUTATION
			if(nunsat > 0)
			{
				// lingering is a boolean int to see whether there is still work to be done
				// at the focal node [index i] after this mutation event
				lingering = (nunsat > 1 || nref > 0);
				// mutation conditions
				if( (lingering && n > 1) || !lingering )
				{
					for(j = 0; j < numChildren; j++)
					{
						if(curConfig->satisfied[i][j] == 0)
						{
							DatConfig_copy_config(curConfig, &deriv);
							mutIdx = idxToNode[i]->children[j]->mut;
							deriv.positions[i]--;
							deriv.positions[mutIdx]++;
							deriv.active[mutIdx] = 1;
							deriv.satisfied[i][j] = 1;
							if(nunsat == 1 && deriv.positions[i] == nref)
								deriv.active[i] = 0;
							//fprintf(stdout, "curConfig->prob = %f, transitionProb = %f\n", curConfig->prob, (double)(ds->theta)/(double)(ds->theta + deriv.positions[i]+1.0-1.0));
							for(k = 0; k < ds->numThetas; k++)
							{
								transitionProb = (double)(ds->thetas[k]/2.0)/(double)(ntot * ds->thetas[k]/2.0 + ntot*(ntot-1)/2);
								//fprintf(stdout, "Mutating : curConfig->prob = %f, transitionProb = %f\n", curConfig->prob, transitionProb);
								//deriv.prob = curConfig->prob * transitionProb;
								deriv.probs[k] = curConfig->probs[k] * transitionProb;
							}
							ConfigCollection_add_config(recipient, &deriv);
						}
					}
				}
			}
		}
	}
	DatConfig_free(&deriv);
	return;
}

int32_t DataSet_binomial_coeff_(int32_t n, int32_t k)
{
	int32_t r = 1, d = n - k;
	if(d > k)
	{
		k = d;
		d = n - k;
	}
	while(n > k)
	{
		if(r >= INT32T_MAX / n)
			return 0;
		r *= n--;
		while(d > 1 && !(r % d))
			r /= d--;
	}
	return r;
}

int32_t DataSet_get_prob_multiplier(BMat * bmat)
{
	int32_t i, mult = 1, numRemaining = 0;
	int32_t * numDuplicates = (int32_t *)malloc(sizeof(int32_t) * bmat->nrows);
	CHECKPOINTER(numDuplicates);
	int32_t numUnique;
	BMat_get_haplotype_counts(bmat, numDuplicates, &numUnique);
	//fprintf(stdout, "numDuplicates: ");
	for(i = 0; i < numUnique; i++)
	{
		numRemaining += numDuplicates[i];
		//fprintf(stdout, "%i ", numDuplicates[i]);
	}
	//fprintf(stdout, "\n");
	for(i = 0; i < numUnique; i++)
	{
		mult *= DataSet_binomial_coeff_(numRemaining, numDuplicates[i]);
		numRemaining -= numDuplicates[i];
	}
	free(numDuplicates);
	return mult;
}
