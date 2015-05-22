#include <stdlib.h>
#include <stdint.h>
#include "definitions.h"
#include "bmat.h"
#include "node.h"
#include "nodelist.h"
#include "datconfig.h"
#include "dataset.h"
#include "hash.h"

void DataSet_init(DataSet * ds, BMat * inputbmat, int32_t numThetas, double * thetas, int32_t printAll, int32_t ordered)
{
	int32_t i,j,k, zeroFirst, numStages;
	ds->bmat = inputbmat;
	ds->numSegSites = inputbmat->ncols;
	ds->numSamples = inputbmat->nrows;
	numStages = ds->numSegSites + ds->numSamples;
	ds->numThetas = numThetas;
	ds->thetas = thetas;

    ds->printAll = printAll;

	ds->collection[0] = (ConfigCollection *)malloc(sizeof(ConfigCollection));
	CHECKPOINTER(ds->collection[0]);
	ds->collection[1] = (ConfigCollection *)malloc(sizeof(ConfigCollection));
	CHECKPOINTER(ds->collection[1]);
	/* there are one more nodes in the phylogeny than
	 * segregating sites. */
	BMat_order_columns(ds->bmat);
    ds->ordered = ordered;

    // check for and then create a perfect phylogeny according to the
    // algorithms of Gusfield (1991) "Efficient Algorithms for Inferring
    // Evolutionary Trees."
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
	NodeList_get_num_children(&(ds->nodeList));
	NodeList_get_idxToNode(&(ds->nodeList));
	DatConfig_init(&(ds->refConfig), ds->bmat->ncols+1, ds->nodeList.numChildren);
	DatConfig_get_ref_config(ds->bmat, &(ds->refConfig));


	DatConfig rootConfig;
	DatConfig_init(&rootConfig, ds->bmat->ncols+1, ds->nodeList.numChildren);
	DatConfig_set_root_config(&rootConfig);

    ds->initialNodes = (int32_t *)malloc(sizeof(int32_t) * (size_t)(ds->refConfig.length));
    CHECKPOINTER(ds->initialNodes);
    DatConfig_set_initial_node_indices(&ds->refConfig, ds->initialNodes);

	/* initialize the ConfigCollections, and then add the root configuration to
	 * the first Collection, and off we go! */
	ConfigCollection_init(ds->collection[0], ds->bmat->ncols+1, ds->nodeList.numChildren);
	ConfigCollection_init(ds->collection[1], ds->bmat->ncols+1, ds->nodeList.numChildren);
	ConfigCollection_add_config(ds->collection[0], &rootConfig);
	zeroFirst = 0;
	for(i = 0; i < numStages-1; i++)
	{
		DataSet_transfer_config_collections(ds->collection[zeroFirst], ds->collection[!zeroFirst], ds);
		zeroFirst = !zeroFirst;
	}


    if(!ds->printAll)
    {
        if(!ds->ordered)
            ds->probMultiplier = DataSet_get_prob_multiplier(ds->bmat);
        else
            ds->probMultiplier = 1.0;
        double finalProb;
        for(k = 0; k < ds->numThetas; k++)
        {
            finalProb = ConfigCollection_get_final_prob(ds->collection[zeroFirst], k) * (double)ds->probMultiplier;
            fprintf(stdout, "%.16e\n", finalProb);
        }
    }

	probMultiplier = DataSet_get_prob_multiplier(ds->bmat);
	finalProb = ConfigCollection_get_final_prob(ds->collection[zeroFirst]) * (double)probMultiplier;

	fprintf(stdout, "prob multiplier: %i\n", probMultiplier);
	fprintf(stdout, "Final prob: %e\n", finalProb);

	ConfigCollection_free(ds->collection[0]);
	ConfigCollection_free(ds->collection[1]);
	free(ds->collection[0]);
	free(ds->collection[1]);
	DatConfig_free(&(ds->refConfig));
	DatConfig_free(&rootConfig);
    free(ds->initialNodes);
	return;
}

void DataSet_free(DataSet * ds)
{
	BMat_free(ds->Lij);
	NodeList_free(&(ds->nodeList));
	free(ds->Lij);
	free(ds->Lj);
    free(ds->thetas);
	return;
}

void DataSet_print_good_probabilities(ConfigCollection * collection, DataSet * ds)
{
    int32_t i, k, good;
    DatConfig * config;
    for(i = 0; i < collection->curNumConfigs; i++)
    {
        config = collection->configs[i];
        good = 1;
        for(k = 0; k < config->length; k++)
        {
            if( (config->positions[k] > 0 && ds->initialNodes[k] == 0) || (config->positions[k] == 0 && ds->initialNodes[k] == 1) )
            {
                good = 0;
                break;
            }
        }
        if(good)
        {
            // must calculate multinomial multiplier for each good DatConfig
            double probMultiplier;
            if(!ds->ordered)
                probMultiplier = (double)DataSet_get_prob_multiplier2(config);
            else
                probMultiplier = 1.0;
            for(k = 0; k < config->length-1; k++)
                printf("%i ", config->positions[k]);
            printf("%i; ", config->positions[config->length-1]);
            for(k = 0; k < ds->numThetas-1; k++)
                printf("%.16e ", config->probs[k] * probMultiplier);
            printf("%.16e\n", config->probs[ds->numThetas-1] * probMultiplier);
        }
    }
    return;
}

void DataSet_transfer_config_collections(ConfigCollection * donor, ConfigCollection * recipient, DataSet * ds)
{
	int32_t i;
	ConfigCollection_reset(recipient);
	for(i = 0; i < donor->curNumConfigs; i++)
		DataSet_donate_deriv_configs(donor->configs[i], recipient, ds->nodeList.idxToNode, ds);
    if(ds->printAll)
        DataSet_print_good_probabilities(recipient, ds);
	HashTable_reset(&(donor->hashTable));
	return;
}

void DataSet_donate_deriv_configs(DatConfig * curConfig, ConfigCollection * recipient, Node ** idxToNode, DataSet * ds)
{
	int32_t i, j, n, ntot, ntotCoal, nref, nunsat, numChildren, mutIdx, lingering;
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
				transitionProb =  (double)((n+1.0)*n)/(double)((ntotCoal)*ds->theta + ntotCoal*(ntotCoal-1.0));
				deriv.prob = curConfig->prob * transitionProb;
				//fprintf(stdout, "Coalescing : curConfig->prob = %f, transitionProb = %f\n", curConfig->prob, transitionProb);
				if(nunsat == 0 && nref == deriv.positions[i]) 	// check for completion
					deriv.active[i] = 0;
				ConfigCollection_add_config(recipient, &deriv);		// this function adds the probability.
                
                /*
                // check positions, if only the initial positions are filled, print the probability
                int32_t good = 1;
                for(k = 0; k < deriv.length; k++)
                {
                    if( (deriv.positions[k] > 0 && ds->initialNodes[k] == 0) || (deriv.positions[k] == 0 && ds->initialNodes[k] == 1) )
                    {
                        good = 0;
                        break;
                    }
                }
                if(good)
                {
                    for(k = 0; k < deriv.length; k++)
                        printf("%i ", deriv.positions[k]);
                    for(k = 0; k < ds->numThetas-1; k++)
                        printf("%f ", deriv.probs[k] * ds->probMultiplier);
                    printf("%f\n", deriv.probs[ds->numThetas-1] * ds->probMultiplier);
                }
                */


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
							transitionProb = (double)(ds->theta/2.0)/(double)(ntot * ds->theta/2.0 + ntot*(ntot-1)/2);
							//fprintf(stdout, "Mutating : curConfig->prob = %f, transitionProb = %f\n", curConfig->prob, transitionProb);
							deriv.prob = curConfig->prob * transitionProb;
							ConfigCollection_add_config(recipient, &deriv);
                            /*
                            int32_t good = 1;
                            for(k = 0; k < deriv.length; k++)
                            {
                                if( (deriv.positions[k] > 0 && ds->initialNodes[k] == 0) || (deriv.positions[k] == 0 && ds->initialNodes[k] == 1) )
                                {
                                    good = 0;
                                    break;
                                }
                            }
                            if(good)
                            {
                                for(k = 0; k < deriv.length; k++)
                                    printf("%i ", deriv.positions[k]);
                                for(k = 0; k < ds->numThetas-1; k++)
                                    printf("%f ", deriv.probs[k] * ds->probMultiplier);
                                printf("%f\n", deriv.probs[ds->numThetas-1] * ds->probMultiplier);
                            }
                            */
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

// these two functions just return the multinomial coefficient of the integers
// in refconfig->positions.

int32_t DataSet_get_prob_multiplier(BMat * bmat)
{
	int32_t i, mult = 1, numRemaining = 0;
	int32_t * numDuplicates = (int32_t *)malloc(sizeof(int32_t) * bmat->nrows);
	CHECKPOINTER(numDuplicates);
	int32_t numUnique;
	BMat_get_haplotype_counts(bmat, numDuplicates, &numUnique);
	for(i = 0; i < numUnique; i++)
		numRemaining += numDuplicates[i];
	for(i = 0; i < numUnique; i++)
	{
		mult *= DataSet_binomial_coeff_(numRemaining, numDuplicates[i]);
		numRemaining -= numDuplicates[i];
	}
	free(numDuplicates);
	return mult;
}

int32_t DataSet_get_prob_multiplier2(DatConfig * refconfig)
{
    int32_t i, mult = 1, totalHaps = 0, numRemaining;
    for(i = 0; i < refconfig->length; i++)
        totalHaps += refconfig->positions[i];
    numRemaining = totalHaps;
    for(i = 0; i < refconfig->length; i++)
    {
        if(refconfig->positions[i] == 0)
            continue;
        mult *= DataSet_binomial_coeff_(numRemaining, refconfig->positions[i]);
        numRemaining -= refconfig->positions[i];
    }
    return mult;
}
