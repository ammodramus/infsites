#include <stdlib.h>
#include <stdint.h>
#include "definitions.h"
#include "bmat.h"
#include "bmat2d.h"
#include "node.h"
#include "nodelist.h"
#include "datconfig.h"
#include "datconfig2d.h"
#include "dataset.h"
#include "dataset2d.h"
#include "hash.h"

void DataSet2d_init(DataSet2d * ds, BMat2d * inputbmat, int numThetas, double * thetas, int numMigRates, double * migRates, int printAll, int ordered, int genetree)
{
	int i,j, numStages, probMultiplier, configLength, finalIdx;
	double finalProb;
	ds->bmat2d = inputbmat;
	ds->numSegSites = inputbmat->ncols;
	ds->numSamples = inputbmat->nrows;
	ds->recipientCollection = 1;
	numStages = ds->numSegSites + ds->numSamples;

    ds->ordered = ordered;

	ds->thetas = thetas;
	ds->numThetas = numThetas;
	ds->migRates = migRates;
	ds->numMigRates = numMigRates;

	ds->collection[0] = (SuperCollection *)malloc(sizeof(SuperCollection));
	CHECKPOINTER(ds->collection[0]);
	ds->collection[1] = (SuperCollection *)malloc(sizeof(SuperCollection));
	CHECKPOINTER(ds->collection[1]);
	BMat_order_columns(&(ds->bmat2d->bmat));

    if(ds->ordered)
        probMultiplier = 1;
    else // unordered
        probMultiplier = DataSet2d_get_prob_multiplier(ds);

	ds->Lij = (BMat *)malloc(sizeof(BMat));
	CHECKPOINTER(ds->Lij);
	BMat_init(ds->Lij, ds->bmat2d->nrows, ds->bmat2d->ncols);
	for(i = 0; i < ds->bmat2d->nrows; i++)
		for(j = 0; j < ds->bmat2d->ncols; j++)
			ds->Lij->mat[i][j] = BMat_get_Lij(&(ds->bmat2d->bmat), i, j);
	ds->Lj = (int *)malloc(sizeof(int) * (size_t)ds->bmat2d->ncols);
	CHECKPOINTER(ds->Lj);
	for(j = 0; j < ds->bmat2d->ncols; j++)
		ds->Lj[j] = BMat_get_Lj(&(ds->bmat2d->bmat), ds->Lij, j);
	if(!BMat_determine_good(&(ds->bmat2d->bmat), ds->Lij, ds->Lj))
		PERROR("Data does not conform to infinite-sites mutation model.");
	NodeList_create_phylogeny(&(ds->nodeList), &(ds->bmat2d->bmat), ds->Lj);
	configLength = ds->nodeList.numNodes;
	NodeList_get_num_children(&(ds->nodeList));
	NodeList_get_idxToNode(&(ds->nodeList));
	DatConfig_init(&(ds->refConfig), configLength, ds->nodeList.numChildren, ds->numThetas);
	DatConfig_get_ref_config(&(ds->bmat2d->bmat), &(ds->refConfig));

    // get initialNodes
    ds->initialNodes = (int *)malloc(sizeof(int) * (size_t)(ds->refConfig.length));
    CHECKPOINTER(ds->initialNodes);
    DatConfig_set_initial_node_indices(&ds->refConfig, ds->initialNodes);

	DatConfig2d_init(&(ds->refConfig2d), configLength, &(ds->refConfig), (SuperConfig *)NULL, ds->numThetas, ds->numMigRates); 
	DatConfig2d_get_ref_datconfig2d(ds->bmat2d, &(ds->refConfig2d));

	DatConfig rootPanmictic;
	DatConfig_init(&rootPanmictic, configLength, ds->nodeList.numChildren, ds->numThetas);
	DatConfig_set_root_config(&rootPanmictic);


	SuperCollection_init(ds->collection[0], configLength, ds->refConfig.numChildren, ds->numThetas, ds->numMigRates);
	SuperCollection_init(ds->collection[1], configLength, ds->refConfig.numChildren, ds->numThetas, ds->numMigRates);
	SuperCollection_add_SuperConfig(ds->collection[0], &rootPanmictic, ds);

	for(i = 0; i < ds->numThetas; i++)
	{
		for(j = 0; j < ds->numMigRates; j++)
		{
			ds->collection[0]->superConfigs[0]->configs2d[0]->probs[i][j] = 1.0;
			ds->collection[0]->superConfigs[0]->configs2d[1]->probs[i][j] = 1.0;
		}
	}

	for(j = 0; j < numStages-1; j++)
		DataSet2d_iterate_stages(ds->collection[!(ds->recipientCollection)], ds->collection[ds->recipientCollection], ds, printAll);

    // multiply theta by 2.0 again if mimicking genetree
    double thetaMultiplier = genetree ? 2.0 : 1.0;
    if(!printAll)
    {
        finalIdx = SuperConfig_get_index(ds->refConfig2d.positions, ds->collection[!(ds->recipientCollection)]->superConfigs[0]->positionMultipliers, configLength);
        printf("theta\tM\tprob\n");
        for(i = 0; i < ds->numThetas; i++)
        {
            for(j = 0; j < ds->numMigRates; j++)
            {
                finalProb = ds->collection[!(ds->recipientCollection)]->superConfigs[0]->configs2d[finalIdx]->probs[i][j] * (double)probMultiplier;
                printf("%f\t%f\t%.16e\n", thetaMultiplier * ds->thetas[i], ds->migRates[j], finalProb);
            }
        }
    }
	SuperCollection_reset(ds->collection[!(ds->recipientCollection)]);
	SuperCollection_reset(ds->collection[(ds->recipientCollection)]);
	DatConfig_free(&rootPanmictic);
	DataSet2d_free(ds);
	return;
}

void DataSet2d_solve_ctypes(DataSet2d * ds, BMat2d * inputbmat, int numThetas, double * thetas, int numMigRates, double * migRates, int ordered, int genetree, double * samplingProbs)
{
	int i,j, numStages, probMultiplier, configLength, finalIdx;
	double finalProb;
	ds->bmat2d = inputbmat;
	ds->numSegSites = inputbmat->ncols;
	ds->numSamples = inputbmat->nrows;
	ds->recipientCollection = 1;
	numStages = ds->numSegSites + ds->numSamples;

    ds->ordered = ordered;

	ds->thetas = thetas;
	ds->numThetas = numThetas;
	ds->migRates = migRates;
	ds->numMigRates = numMigRates;

	ds->collection[0] = (SuperCollection *)malloc(sizeof(SuperCollection));
	CHECKPOINTER(ds->collection[0]);
	ds->collection[1] = (SuperCollection *)malloc(sizeof(SuperCollection));
	CHECKPOINTER(ds->collection[1]);
	BMat_order_columns(&(ds->bmat2d->bmat));

    if(ds->ordered)
        probMultiplier = 1;
    else
        probMultiplier = DataSet2d_get_prob_multiplier(ds);

	ds->Lij = (BMat *)malloc(sizeof(BMat));
	CHECKPOINTER(ds->Lij);
	BMat_init(ds->Lij, ds->bmat2d->nrows, ds->bmat2d->ncols);
	for(i = 0; i < ds->bmat2d->nrows; i++)
		for(j = 0; j < ds->bmat2d->ncols; j++)
			ds->Lij->mat[i][j] = BMat_get_Lij(&(ds->bmat2d->bmat), i, j);
	ds->Lj = (int *)malloc(sizeof(int) * (size_t)ds->bmat2d->ncols);
	CHECKPOINTER(ds->Lj);
	for(j = 0; j < ds->bmat2d->ncols; j++)
		ds->Lj[j] = BMat_get_Lj(&(ds->bmat2d->bmat), ds->Lij, j);
	if(!BMat_determine_good(&(ds->bmat2d->bmat), ds->Lij, ds->Lj))
		PERROR("Data does not conform to infinite-sites mutation model.");
	NodeList_create_phylogeny(&(ds->nodeList), &(ds->bmat2d->bmat), ds->Lj);
	configLength = ds->nodeList.numNodes;
	NodeList_get_num_children(&(ds->nodeList));
	NodeList_get_idxToNode(&(ds->nodeList));
	DatConfig_init(&(ds->refConfig), configLength, ds->nodeList.numChildren, ds->numThetas);
	DatConfig_get_ref_config(&(ds->bmat2d->bmat), &(ds->refConfig));

    // get initialNodes
    ds->initialNodes = (int *)malloc(sizeof(int) * (size_t)(ds->refConfig.length));
    CHECKPOINTER(ds->initialNodes);
    DatConfig_set_initial_node_indices(&ds->refConfig, ds->initialNodes);

	DatConfig2d_init(&(ds->refConfig2d), configLength, &(ds->refConfig), (SuperConfig *)NULL, ds->numThetas, ds->numMigRates); 
	DatConfig2d_get_ref_datconfig2d(ds->bmat2d, &(ds->refConfig2d));

	DatConfig rootPanmictic;
	DatConfig_init(&rootPanmictic, configLength, ds->nodeList.numChildren, ds->numThetas);
	DatConfig_set_root_config(&rootPanmictic);


	SuperCollection_init(ds->collection[0], configLength, ds->refConfig.numChildren, ds->numThetas, ds->numMigRates);
	SuperCollection_init(ds->collection[1], configLength, ds->refConfig.numChildren, ds->numThetas, ds->numMigRates);
	SuperCollection_add_SuperConfig(ds->collection[0], &rootPanmictic, ds);

	for(i = 0; i < ds->numThetas; i++)
	{
		for(j = 0; j < ds->numMigRates; j++)
		{
			ds->collection[0]->superConfigs[0]->configs2d[0]->probs[i][j] = 1.0;
			ds->collection[0]->superConfigs[0]->configs2d[1]->probs[i][j] = 1.0;
		}
	}

	for(j = 0; j < numStages-1; j++)
		DataSet2d_iterate_stages(ds->collection[!(ds->recipientCollection)], ds->collection[ds->recipientCollection], ds, 0); // last argument (printAll) is zero in this ctypes function

    finalIdx = SuperConfig_get_index(ds->refConfig2d.positions, ds->collection[!(ds->recipientCollection)]->superConfigs[0]->positionMultipliers, configLength);
    int spIdx = 0;
    for(i = 0; i < ds->numThetas; i++)
    {
        for(j = 0; j < ds->numMigRates; j++)
        {
            finalProb = ds->collection[!(ds->recipientCollection)]->superConfigs[0]->configs2d[finalIdx]->probs[i][j] * (double)probMultiplier;
            samplingProbs[spIdx] = finalProb;
            spIdx++;
        }
    }

	SuperCollection_reset(ds->collection[!(ds->recipientCollection)]);
	SuperCollection_reset(ds->collection[(ds->recipientCollection)]);
	DatConfig_free(&rootPanmictic);
	DataSet2d_free(ds);
	return;
}

void DataSet2d_solve_ctypes_all(DataSet2d * ds, BMat2d * inputbmat, int numThetas, double * thetas, int numMigRates, double * migRates, int ordered, int genetree, int ** recIdxs, double ** samplingProbs)
{
	int i,j, numStages, configLength;
	ds->bmat2d = inputbmat;
	ds->numSegSites = inputbmat->ncols;
	ds->numSamples = inputbmat->nrows;
	ds->recipientCollection = 1;
	numStages = ds->numSegSites + ds->numSamples;

    ds->ordered = ordered;

	ds->thetas = thetas;
	ds->numThetas = numThetas;
	ds->migRates = migRates;
	ds->numMigRates = numMigRates;

	ds->collection[0] = (SuperCollection *)malloc(sizeof(SuperCollection));
	CHECKPOINTER(ds->collection[0]);
	ds->collection[1] = (SuperCollection *)malloc(sizeof(SuperCollection));
	CHECKPOINTER(ds->collection[1]);
	BMat_order_columns(&(ds->bmat2d->bmat));
	ds->Lij = (BMat *)malloc(sizeof(BMat));
	CHECKPOINTER(ds->Lij);
	BMat_init(ds->Lij, ds->bmat2d->nrows, ds->bmat2d->ncols);
	for(i = 0; i < ds->bmat2d->nrows; i++)
		for(j = 0; j < ds->bmat2d->ncols; j++)
			ds->Lij->mat[i][j] = BMat_get_Lij(&(ds->bmat2d->bmat), i, j);
	ds->Lj = (int *)malloc(sizeof(int) * (size_t)ds->bmat2d->ncols);
	CHECKPOINTER(ds->Lj);
	for(j = 0; j < ds->bmat2d->ncols; j++)
		ds->Lj[j] = BMat_get_Lj(&(ds->bmat2d->bmat), ds->Lij, j);
	if(!BMat_determine_good(&(ds->bmat2d->bmat), ds->Lij, ds->Lj))
		PERROR("Data does not conform to infinite-sites mutation model.");
	NodeList_create_phylogeny(&(ds->nodeList), &(ds->bmat2d->bmat), ds->Lj);
	configLength = ds->nodeList.numNodes;
	NodeList_get_num_children(&(ds->nodeList));
	NodeList_get_idxToNode(&(ds->nodeList));
	DatConfig_init(&(ds->refConfig), configLength, ds->nodeList.numChildren, ds->numThetas);
	DatConfig_get_ref_config(&(ds->bmat2d->bmat), &(ds->refConfig));

    // get initialNodes
    ds->initialNodes = (int *)malloc(sizeof(int) * (size_t)(ds->refConfig.length));
    CHECKPOINTER(ds->initialNodes);
    DatConfig_set_initial_node_indices(&ds->refConfig, ds->initialNodes);

	DatConfig2d_init(&(ds->refConfig2d), configLength, &(ds->refConfig), (SuperConfig *)NULL, ds->numThetas, ds->numMigRates); 
	DatConfig2d_get_ref_datconfig2d(ds->bmat2d, &(ds->refConfig2d));

	DatConfig rootPanmictic;
	DatConfig_init(&rootPanmictic, configLength, ds->nodeList.numChildren, ds->numThetas);
	DatConfig_set_root_config(&rootPanmictic);


	SuperCollection_init(ds->collection[0], configLength, ds->refConfig.numChildren, ds->numThetas, ds->numMigRates);
	SuperCollection_init(ds->collection[1], configLength, ds->refConfig.numChildren, ds->numThetas, ds->numMigRates);
	SuperCollection_add_SuperConfig(ds->collection[0], &rootPanmictic, ds);

	for(i = 0; i < ds->numThetas; i++)
	{
		for(j = 0; j < ds->numMigRates; j++)
		{
			ds->collection[0]->superConfigs[0]->configs2d[0]->probs[i][j] = 1.0;
			ds->collection[0]->superConfigs[0]->configs2d[1]->probs[i][j] = 1.0;
		}
	}

    int curRecIdx = 0;
	for(j = 0; j < numStages-1; j++)
		DataSet2d_iterate_stages_ctypes(ds->collection[!(ds->recipientCollection)], ds->collection[ds->recipientCollection], ds, recIdxs, samplingProbs, &curRecIdx);

	SuperCollection_reset(ds->collection[!(ds->recipientCollection)]);
	SuperCollection_reset(ds->collection[(ds->recipientCollection)]);
	DatConfig_free(&rootPanmictic);
	DataSet2d_free(ds);
	return;
}

void DataSet2d_free(DataSet2d * ds)
{
	BMat_free(ds->Lij);
	DatConfig_free(&(ds->refConfig));
	SuperCollection_free(ds->collection[0]);
	SuperCollection_free(ds->collection[1]);
	NodeList_free(&(ds->nodeList));
	DatConfig2d_free(&(ds->refConfig2d));
    if(ds->initialNodes)
        free(ds->initialNodes);
	free(ds->collection[0]);
	free(ds->collection[1]);
	free(ds->Lij);
	free(ds->Lj);
	return;
}

void DataSet2d_print_good_probabilities(SuperCollection * recipient, DataSet2d * ds)
{
    int i, j, k, l, good;
    DatConfig * curConfig;
    SuperConfig * curSuperConfig;
    DatConfig2d * curConfig2d;
    int numPositions;
    double probMultiplier;
    for(i = 0; i < recipient->curNumSuperConfigs; i++)
    {
        curConfig = &(recipient->superConfigs[i]->panmictic);
        curSuperConfig = recipient->superConfigs[i];
        numPositions = curConfig->length;
        good = 1;
        for(j = 0; j < curConfig->length; j++)
        {
            if((curConfig->positions[j] > 0 && ds->initialNodes[j] == 0) || (curConfig->positions[j] == 0 && ds->initialNodes[j] == 1))
            {
                good = 0;
                break;
            }
        }
        if(good)
        {
            for(j = 0; j < curSuperConfig->numConfigs2d; j++)
            {
                curConfig2d = curSuperConfig->configs2d[j];
                if(!ds->ordered)
                    probMultiplier = (double)DataSet2d_get_prob_multiplier2(curConfig2d);
                else
                    probMultiplier = 1.0;

                for(k = 0; k < numPositions-1; k++)
                    printf("%i ", curConfig2d->positions[k][0]);
                printf("%i|", curConfig2d->positions[numPositions-1][0]);
                for(k = 0; k < numPositions-1; k++)
                    printf("%i ", curConfig2d->positions[k][1]);
                printf("%i;", curConfig2d->positions[numPositions-1][1]);
                for(k = 0; k < ds->numThetas; k++)
                {
                    for(l = 0; l < ds->numMigRates; l++)
                    {
                        if(k == ds->numThetas-1 && l == ds->numMigRates-1)
                            printf("%.16e\n", curConfig2d->probs[k][l] * probMultiplier);
                        else
                            printf("%.16e ", curConfig2d->probs[k][l] * probMultiplier);
                    }
                }
            }
        }
    }
    return;
}

void DataSet2d_record_good_probabilities_ctypes(SuperCollection * recipient, DataSet2d * ds, int ** recIdxs, double ** samplingProbs, int * p_curRecIdx)
{
    int i, j, k, l, good, spIdx;
    DatConfig * curConfig;
    SuperConfig * curSuperConfig;
    DatConfig2d * curConfig2d;
    int numPositions;
    double probMultiplier;
    for(i = 0; i < recipient->curNumSuperConfigs; i++)
    {
        curConfig = &(recipient->superConfigs[i]->panmictic);
        curSuperConfig = recipient->superConfigs[i];
        numPositions = curConfig->length;
        good = 1;
        for(j = 0; j < curConfig->length; j++)
        {
            if((curConfig->positions[j] > 0 && ds->initialNodes[j] == 0) || (curConfig->positions[j] == 0 && ds->initialNodes[j] == 1))
            {
                good = 0;
                break;
            }
        }
        if(good)
        {
            for(j = 0; j < curSuperConfig->numConfigs2d; j++)
            {
                curConfig2d = curSuperConfig->configs2d[j];
                if(!ds->ordered)
                    probMultiplier = (double)DataSet2d_get_prob_multiplier2(curConfig2d);
                else
                    probMultiplier = 1.0;

                for(k = 0; k < numPositions; k++)
                {
                    recIdxs[*p_curRecIdx][k] = curConfig2d->positions[k][0];
                    recIdxs[*p_curRecIdx][k+numPositions] = curConfig2d->positions[k][1];
                }
                spIdx = 0;
                for(k = 0; k < ds->numThetas; k++)
                {
                    for(l = 0; l < ds->numMigRates; l++)
                    {
                        samplingProbs[*p_curRecIdx][spIdx] = curConfig2d->probs[k][l] * probMultiplier;
                        spIdx++;
                    }
                }
                (*p_curRecIdx)++;
            }
        }
    }
    return;
}

void DataSet2d_iterate_stages(SuperCollection * donor, SuperCollection * recipient, DataSet2d * ds, int printAll)
{
	DataSet2d_transfer_config_collections(donor, recipient, ds);
	DataSet2d_link_supercollections(donor, recipient, ds);
	DataSet2d_solve_equations(recipient);
    if(printAll)
        DataSet2d_print_good_probabilities(recipient, ds);
	ds->recipientCollection = !(ds->recipientCollection);
	return;
}

void DataSet2d_iterate_stages_ctypes(SuperCollection * donor, SuperCollection * recipient, DataSet2d * ds, int ** recIdxs, double ** samplingProbs, int * p_curRecIdx)
{
	DataSet2d_transfer_config_collections(donor, recipient, ds);
	DataSet2d_link_supercollections(donor, recipient, ds);
	DataSet2d_solve_equations(recipient);
    DataSet2d_record_good_probabilities_ctypes(recipient, ds, recIdxs, samplingProbs, p_curRecIdx);
	ds->recipientCollection = !(ds->recipientCollection);
	return;
}

void DataSet2d_solve_equations(SuperCollection * collection)
{
	int i, j, k, l;
	for(k = 0; k < collection->numThetas; k++)
	{
		for(l = 0; l < collection->numMigRates; l++)
		{
			for(i = 0; i < collection->curNumSuperConfigs; i++)
			{
				SuperEquations_solve(&(collection->superConfigs[i]->eqs[k][l]));
				for(j = 0; j < collection->superConfigs[i]->numConfigs2d; j++)
				{
#ifndef CSPARSECOMPILE
					//printf("collection->superConfigs[%i]->eq.x[%i] = %f\n", i, j, collection->superConfigs[i]->eq.x[j]);
					collection->superConfigs[i]->configs2d[j]->probs[k][l] = collection->superConfigs[i]->eqs[k][l].x[j];
#endif
#ifdef CSPARSECOMPILE
					//printf("collection->superConfigs[%i]->eq.b[%i] = %f\n", i, j, collection->superConfigs[i]->eq.b[j]);
					collection->superConfigs[i]->configs2d[j]->probs[k][l] = collection->superConfigs[i]->eqs[k][l].b[j];   //CSPARSE EDIT
#endif
				}
			}
		}
	}
	return;
}

void DataSet2d_donate_deriv_superconfigs(SuperConfig * super, SuperCollection * recipient, DataSet2d * ds)
{
	int i, j, n, ntot, nref, nunsat, numChildren, mutIdx, lingering;
	DatConfig * panmictic = &(super->panmictic);
	DatConfig deriv;
	DatConfig_copy_dimensions(panmictic, &deriv);
	ntot = 0;
	for(j = 0; j < panmictic->length; j++)
		ntot += panmictic->positions[j];
	for(i = 0; i < super->configLength; i++)
	{
		numChildren = super->numChildren[i];
		if(panmictic->active[i])
		{
			//////////////
			// COALESCENCE
			nunsat = 0;
			n = panmictic->positions[i];
			for(j = 0; j < numChildren; j++)
				if(panmictic->satisfied[i][j] == 0)
					nunsat++;
			nref = ds->refConfig.positions[i];
			// coalescence conditions
			if(n < nref + nunsat)
			{
				DatConfig_copy_config(panmictic, &deriv);
				deriv.positions[i] += 1;
				//transitionProb =  (double)((n+1.0)*n)/(double)((ntotCoal)*ds->theta + ntotCoal*(ntotCoal-1.0));
				//deriv.prob = curConfig->prob * transitionProb;
				//fprintf(stdout, "Coalescing : curConfig->prob = %f, transitionProb = %f\n", curConfig->prob, transitionProb);
				if(nunsat == 0 && nref == deriv.positions[i]) 	// check for completion
					deriv.active[i] = 0;
				//ConfigCollection_add_config(recipient, &deriv);		// this function adds the probability.
				SuperCollection_add_SuperConfig(recipient, &deriv, ds);
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
						if(panmictic->satisfied[i][j] == 0)
						{
							DatConfig_copy_config(panmictic, &deriv);
							mutIdx = ds->nodeList.idxToNode[i]->children[j]->mut;
							deriv.positions[i]--;
							deriv.positions[mutIdx]++;
							deriv.active[mutIdx] = 1;
							deriv.satisfied[i][j] = 1;
							if(nunsat == 1 && deriv.positions[i] == nref)
								deriv.active[i] = 0;
							SuperCollection_add_SuperConfig(recipient, &deriv, ds);
						}
					}
				}
			}
		}
	}
	DatConfig_free(&deriv);
	return;
}

void DataSet2d_transfer_config_collections(SuperCollection * donor, SuperCollection * recipient, DataSet2d * ds)
{
	int i;
	SuperCollection_reset(recipient);
	for(i = 0; i < donor->curNumSuperConfigs; i++)
		DataSet2d_donate_deriv_superconfigs(donor->superConfigs[i], recipient, ds);
	HashSuper_reset(&(donor->hashSuper));
	return;
}

// *linking probabilities* doesn't require the duplication of DatConfig2d's... just the 2d-positions.
void DataSet2d_link_probabilities(DatConfig2d * config, SuperCollection * recipient, DataSet2d * ds)
{
	int i, j, k, thetaIdx, migRateIdx, n, n0, n1, ntot, ntot0, ntot1, ntotCoal, ntotCoal0, ntotCoal1, nref, nunsat, curNumChildren, mutIdx, lingering, twoDemeIdx;
	int configLength = config->panmictic->length;
	int * numChildren = config->panmictic->numChildren;
	double totalRate, transitionProb;
	SuperConfig * superConfig;
	twoints * positions2d = (twoints *)malloc(sizeof(twoints) * configLength);
	CHECKPOINTER(positions2d);
	int * positions = (int *)malloc(sizeof(int) * configLength);
	CHECKPOINTER(positions);
	HashSuper * hashSuper = &(ds->collection[ds->recipientCollection]->hashSuper);
	ntot = ntot0 = ntot1 = 0;
	for(j = 0; j < configLength; j++)
	{
		ntot += config->panmictic->positions[j];
		ntot0 += config->positions[j][0];
		ntot1 += config->positions[j][1];
	}
	ntotCoal = ntot + 1;
	ntotCoal0 = ntot0 + 1;
	ntotCoal1 = ntot1 + 1;
	for(i = 0; i < configLength; i++)
	{
		curNumChildren = numChildren[i];
		if(config->panmictic->active[i])
		{
			nunsat = 0;
			n0 = config->positions[i][0];
			n1 = config->positions[i][1];
			n = n0 + n1;
			for(j = 0; j < curNumChildren; j++)
				if(config->panmictic->satisfied[i][j] == 0)
					nunsat++;
			nref = ds->refConfig.positions[i];
			//fprintf(stdout, "nunsat = %i, nref = %i, n = %i\n", nunsat, nref, n);
			// coalescence conditions
			if(n0 > 0)
			{
				////////////////
				// COALESCENCE 0
				if(n < nref + nunsat)
				{
					for(k = 0; k < configLength; k++)
					{
						positions[k] = config->panmictic->positions[k];
						positions2d[k][0] = config->positions[k][0];
						positions2d[k][1] = config->positions[k][1];
					}
					positions2d[i][0] += 1;
					positions[i] += 1;
					superConfig = HashSuper_find_superconfig(positions, hashSuper);
					if(superConfig == NULL)
					{
						printf("debug output (coalescence 0):\n");
						for(k = 0; k < configLength; k++)
							printf("positions[%i] = %i\n", k, positions[k]);
						for(k = 0; k < ds->collection[ds->recipientCollection]->curNumSuperConfigs; k++)
							SuperConfig_print(ds->collection[ds->recipientCollection]->superConfigs[k], stdout);
						PERROR("superConfig not found.");
					}
					twoDemeIdx = SuperConfig_get_index(positions2d, superConfig->positionMultipliers, superConfig->configLength);
					for(thetaIdx = 0; thetaIdx < superConfig->numThetas; thetaIdx++)
					{
						for(migRateIdx = 0; migRateIdx < superConfig->numMigRates; migRateIdx++)
						{
							totalRate = (ntotCoal)*ds->thetas[thetaIdx] + ntotCoal0*(ntotCoal0-1.0) + ntot1*(ntot1-1.0) + ntotCoal0*ds->migRates[migRateIdx] + ntot1*ds->migRates[migRateIdx];
							transitionProb = (n0+1.0)*n0 / totalRate;
							superConfig->eqs[thetaIdx][migRateIdx].b[twoDemeIdx] += transitionProb * config->probs[thetaIdx][migRateIdx];
						}
					}
				}
				/////////////
				// MUTATION 0
				if(nunsat > 0)
				{
					// lingering is a boolean int to see whether there is still work to be done
					// at the focal node [index i] after this mutation event
					lingering = (nunsat > 1 || nref > 0);
					// mutation conditions
					//if( (lingering && n0 > 1) || !lingering )
					if( (lingering && n > 1) || !lingering )
					{
						for(j = 0; j < curNumChildren; j++)
						{
							if(config->panmictic->satisfied[i][j] == 0)
							{
								for(k = 0; k < configLength; k++)
								{
									positions[k] = config->panmictic->positions[k];
									positions2d[k][0] = config->positions[k][0];
									positions2d[k][1] = config->positions[k][1];
								}
								mutIdx = ds->nodeList.idxToNode[i]->children[j]->mut;
								positions[i]--;
								positions2d[i][0]--;
								positions[mutIdx]++;
								positions2d[mutIdx][0]++;
								superConfig = HashSuper_find_superconfig(positions, hashSuper);
								if(superConfig == NULL)
								{
									printf("debug output: (mutation 0)\n");
									for(k = 0; k < configLength; k++)
										printf("positions[%i] = %i\n", k, positions[k]);
									for(k = 0; k < ds->collection[ds->recipientCollection]->curNumSuperConfigs; k++)
										SuperConfig_print(ds->collection[ds->recipientCollection]->superConfigs[k], stdout);
									PERROR("superConfig not found.");
								}
								twoDemeIdx = SuperConfig_get_index(positions2d, superConfig->positionMultipliers, superConfig->configLength);
								//totalRate = (ntotCoal)*ds->theta + ntot0*(ntot0-1.0) + ntot1*(ntot1-1.0) + ntot0*ds->migRatesPair[0] + ntot1*ds->migRatesPair[1];
								for(thetaIdx = 0; thetaIdx < superConfig->numThetas; thetaIdx++)
								{
									for(migRateIdx = 0; migRateIdx < superConfig->numMigRates; migRateIdx++)
									{
										totalRate = (ntot)*ds->thetas[thetaIdx] + ntot0*(ntot0-1.0) + ntot1*(ntot1-1.0) + ntot0*ds->migRates[migRateIdx] + ntot1*ds->migRates[migRateIdx];
										transitionProb = ds->thetas[thetaIdx] / totalRate;
										superConfig->eqs[thetaIdx][migRateIdx].b[twoDemeIdx] += transitionProb * config->probs[thetaIdx][migRateIdx];
									}
								}
							}
						}
					}
				}
			}
			if(n1 > 0)
			{
				////////////////
				// COALESCENCE 1
				if(n < nref + nunsat)
				{
					for(k = 0; k < configLength; k++)
					{
						positions[k] = config->panmictic->positions[k];
						positions2d[k][0] = config->positions[k][0];
						positions2d[k][1] = config->positions[k][1];
					}
					positions2d[i][1] += 1;
					positions[i] += 1;
					superConfig = HashSuper_find_superconfig(positions, hashSuper);
					//for(j = 0; j < configLength; j++)
						//printf("positions[%i] = %i (%i, %i)\n", j, positions[j], positions2d[j][0], positions2d[j][1]);
					if(superConfig == NULL)
					{
						printf("debug output: (coalescence 1)\n");
						for(k = 0; k < configLength; k++)
							printf("positions[%i] = %i\n", k, positions[k]);
						for(j = 0; j < ds->collection[ds->recipientCollection]->curNumSuperConfigs; j++)
							SuperConfig_print(ds->collection[ds->recipientCollection]->superConfigs[j], stdout);
						PERROR("superConfig not found.");
					}
					for(thetaIdx = 0; thetaIdx < superConfig->numThetas; thetaIdx++)
					{
						for(migRateIdx = 0; migRateIdx < superConfig->numMigRates; migRateIdx++)
						{
							totalRate = (ntotCoal)*ds->thetas[thetaIdx] + ntot0*(ntot0-1.0) + ntotCoal1*(ntotCoal1-1.0) + ntot0*ds->migRates[migRateIdx] + ntotCoal1*ds->migRates[migRateIdx];
							twoDemeIdx = SuperConfig_get_index(positions2d, superConfig->positionMultipliers, superConfig->configLength);
							transitionProb = (n1+1.0)*n1 / totalRate;
							superConfig->eqs[thetaIdx][migRateIdx].b[twoDemeIdx] += transitionProb * config->probs[thetaIdx][migRateIdx];
						}
					}
				}
				/////////////
				// MUTATION 1
				if(nunsat > 0)
				{
					// lingering is a boolean int to see whether there is still work to be done
					// at the focal node [index i] after this mutation event
					lingering = (nunsat > 1 || nref > 0);
					// mutation conditions
					//if( (lingering && n1 > 1) || !lingering )
					if( (lingering && n > 1) || !lingering )
					{
						for(j = 0; j < curNumChildren; j++)
						{
							if(config->panmictic->satisfied[i][j] == 0)
							{
								for(k = 0; k < configLength; k++)
								{
									positions[k] = config->panmictic->positions[k];
									positions2d[k][0] = config->positions[k][0];
									positions2d[k][1] = config->positions[k][1];
								}
								mutIdx = ds->nodeList.idxToNode[i]->children[j]->mut;
								positions[i]--;
								positions2d[i][1]--;
								positions[mutIdx]++;
								positions2d[mutIdx][1]++;
								superConfig = HashSuper_find_superconfig(positions, hashSuper);
								if(superConfig == NULL)
								{
									printf("debug output: (mutation 1)\n");
									for(k = 0; k < configLength; k++)
										printf("positions[%i] = %i\n", k, positions[k]);
									for(j = 0; j < ds->collection[ds->recipientCollection]->curNumSuperConfigs; j++)
										SuperConfig_print(ds->collection[ds->recipientCollection]->superConfigs[j], stdout);
									PERROR("superConfig not found.");
								}
								twoDemeIdx = SuperConfig_get_index(positions2d, superConfig->positionMultipliers, superConfig->configLength);
								//totalRate = (ntotCoal)*ds->theta + ntotCoal0*(ntotCoal0-1.0) + ntot0*ds->migRatesPair[0] + ntot1*ds->migRatesPair[1];
								//totalRate = (ntot)*ds->theta + ntot0*(ntot0-1.0) + ntot0*ds->migRatesPair[0] + ntot1*ds->migRatesPair[1];
								for(thetaIdx = 0; thetaIdx < superConfig->numThetas; thetaIdx++)
								{
									for(migRateIdx = 0; migRateIdx < superConfig->numMigRates; migRateIdx++)
									{
										totalRate = (ntot)*ds->thetas[thetaIdx] + ntot0*(ntot0-1.0) + ntot1*(ntot1-1.0) + ntot0*ds->migRates[migRateIdx] + ntot1*ds->migRates[migRateIdx];
										transitionProb = ds->thetas[thetaIdx] / totalRate;
										superConfig->eqs[thetaIdx][migRateIdx].b[twoDemeIdx] += transitionProb * config->probs[thetaIdx][migRateIdx];
									}
								}
							}
						}
					}
				}
			}
		}
	}
	free(positions);
	free(positions2d);
	return;
}

/*
int DataSet_get_prob_multiplier(BMat * bmat)
{
	int i, mult = 1, numRemaining = 0;
	int * numDuplicates = (int *)malloc(sizeof(int) * bmat->nrows);
	CHECKPOINTER(numDuplicates);
	int numUnique;
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
*/
int DataSet2d_binomial_coeff_(int n, int k)
{
	int r = 1, d = n - k;
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

int DataSet2d_get_prob_multiplier(DataSet2d * ds)
{
	int i, deme, mult = 1, numRemaining;
	twoints * numDuplicates = (twoints *)malloc(sizeof(twoints) * ds->bmat2d->bmat.nrows);
	CHECKPOINTER(numDuplicates);
	int numUnique[2];
	BMat2d_get_haplotype_counts(ds->bmat2d, numDuplicates, numUnique);
	for(deme = 0; deme < 2; deme++)
	{
		numRemaining = 0;
		for(i = 0; i < numUnique[deme]; i++)
			numRemaining += numDuplicates[i][deme];
		//fprintf(stdout, "\n");
		for(i = 0; i < numUnique[deme]; i++)
		{
			mult *= DataSet2d_binomial_coeff_(numRemaining, numDuplicates[i][deme]);
			numRemaining -= numDuplicates[i][deme];
		}
	}
	free(numDuplicates);
	return mult;
}

int DataSet2d_get_prob_multiplier2(DatConfig2d * config)
{
    int i, j, deme, mult = 1, totalHaps = 0, numRemaining;

    for(deme = 0; deme < 2; deme++)
    {
        totalHaps = 0;
        for(j = 0; j < config->panmictic->length; j++)
            totalHaps += config->positions[j][deme];
        numRemaining = totalHaps;
        for(i = 0; i < config->panmictic->length; i++)
        {
            if(config->positions[i][deme] == 0)
                continue;
            mult *= DataSet_binomial_coeff_(numRemaining, config->positions[i][deme]);
            numRemaining -= config->positions[i][deme];
        }
    }
    return mult;
}

void DataSet2d_link_datconfig2ds(SuperConfig * donorConfig, SuperCollection * recipient, DataSet2d * ds)
{
	int i, numConfigs2d = donorConfig->numConfigs2d;
	for(i = 0; i < numConfigs2d; i++)
		DataSet2d_link_probabilities(donorConfig->configs2d[i], recipient, ds);
	return;
}

void DataSet2d_link_supercollections(SuperCollection * donor, SuperCollection * recipient, DataSet2d * ds)
{
	int i, numSuperConfigs = donor->curNumSuperConfigs;
	for(i = 0; i < numSuperConfigs; i++)
		DataSet2d_link_datconfig2ds(donor->superConfigs[i], recipient, ds);
	return;
}

// a test main
/*
int main()
{
	BMat2d b2;
	DataSet2d ds;
	double theta = 1.0;
	double * migRatesPair = (double *)malloc(sizeof(double) * 2);
	migRatesPair[0] = 0.5;
	migRatesPair[1] = 0.5;
	theta /= 2.0; 		// to make probabilities maximally compatible with genetree, which defines theta as 4*N_{tot}*mu = 4*N*D*mu, which here is twice 4*N*mu.
	BMat2d_read_input("testfile2d", &b2);
	DataSet2d_init(&ds, &b2, theta, migRatesPair);
	BMat2d_free(&b2);
	free(migRatesPair);
	return 0;
}
*/
