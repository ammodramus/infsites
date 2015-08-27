#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <getopt.h>
#include <errno.h>
#include "definitions.h"
#include "bmat.h"
#include "bmat2d.h"
#include "dataset.h"
#include "dataset2d.h"

typedef struct solveroptions_
{
	char filenameIn[200];
	char filenameOut[200];
	int numDemes;
	double * thetas;
	int numThetas;
	double * migRates;
	int numMigRates;
	FILE * fin;
	FILE * fout;
    int all;
    int ordered;
    int genetree;
} SolverOptions;

static struct option long_options[] =
{
	{"help", no_argument, 0, 'h'},
	{"input", required_argument, 0, 'i'},
	{"numdemes", required_argument, 0, 'D'},
	{"theta-file", required_argument, 0, 't'},
	{"migration-file", required_argument, 0, 'M'},
	{"output", required_argument, 0, 'o'},
	{"ordered", no_argument, 0, 'O'},
	{"unordered", no_argument, 0, 'u'},
	{"no-genetree", no_argument, 0, 1},
	{0, 0, 0, 0}
};

void Solver_print_usage()
{
	fprintf(stderr, "\nUsage: solver [OPTIONS] -i datafile\n");
	return;
}

void SolverOptions_parse_options(int argc, char ** argv, SolverOptions * opt)
{
	int c, i, optionIndex, success, numThetas = -1, numMigRates = -1, thetaIdx, migrationIdx, thetaFileProvided = 0, migFileProvided = 0, stdinSet = 0;
	char ch, * line;
	FILE * thetain, * migrationin;
	double * thetas = NULL, * migRates = NULL;
	opt->numDemes = -1;
    opt->ordered = 1;
    opt->all = 0;
    opt->genetree = 1;
	strcpy(opt->filenameIn, "");
	strcpy(opt->filenameOut, "");
	//opt->theta = -1.0;
	while(1)
	{
		c = getopt_long(argc, argv, "hi:D:t:M:o:saO", long_options, &optionIndex);
		if(c == -1)
			break;
		switch(c)
		{
			case 'h':
				Solver_print_usage();
				break;
			case 'i':
				strncpy(opt->filenameIn, optarg, sizeof(opt->filenameIn));
				break;
			case 'D':
				success = (int)sscanf(optarg, "%i", &(opt->numDemes));
				if(success == 0)
					PERROR("Invalid input for number of demes (--numdemes, -D)");
				break;
			case 'o':
				strncpy(opt->filenameOut, optarg, sizeof(opt->filenameOut));
				break;
			case 't':
				thetaFileProvided = 1;
				// deal with theta-file
				thetain = fopen(optarg, "r");
				if(!thetain)
				{
					fprintf(stderr, "Could not open %s for reading theta parameters.\n", optarg);
					exit(1);
				}
				// count number of lines/theta parameters
				numThetas = 0;
				while(!feof(thetain))
				{
					ch = fgetc(thetain);
					if(ch == '\n')
						numThetas++;
				}
				fseek(thetain, SEEK_SET, 0);
				thetas = (double *)calloc((size_t)numThetas, sizeof(double));
				CHECKPOINTER(thetas);
				// back to beginning of file to set theta parameters
				line = (char *)malloc(sizeof(char) * DEFAULT_MAX_LINE_SIZE);
				CHECKPOINTER(line);
				thetaIdx = 0;
				while(fgets(line, DEFAULT_MAX_LINE_SIZE, thetain) != NULL)
				{
					success = sscanf(line, "%lf\n", &(thetas[thetaIdx++]));
					if(!success)
						PERROR("Invalid theta in theta-file.");
				}
				free(line);
				break;
			case 'M':
				migFileProvided = 1;
				// deal with theta-file
				migrationin = fopen(optarg, "r");
				if(!migrationin)
				{
					fprintf(stderr, "Could not open %s for reading migration parameters.\n", optarg);
					exit(1);
				}
				// count number of lines/migration parameters
				numMigRates = 0;
				while(!feof(migrationin))
				{
					ch = fgetc(migrationin);
					if(ch == '\n')
						numMigRates++;
				}
				fseek(migrationin, SEEK_SET, 0);
				migRates = (double *)calloc((size_t)numMigRates, sizeof(double));
				CHECKPOINTER(migRates);
				// back to beginning of file to set migration parameters
				line = (char *)malloc(sizeof(char) * DEFAULT_MAX_LINE_SIZE);
				CHECKPOINTER(line);
				migrationIdx = 0;
				while(fgets(line, DEFAULT_MAX_LINE_SIZE, migrationin) != NULL)
				{
					success = sscanf(line, "%lf\n", &(migRates[migrationIdx++]));
					if(!success)
						PERROR("Invalid migration in migration-file.");
				}
				free(line);
				break;
            case 's':
                // stdin input
                stdinSet = 1;
                break;
            case 'a':
                opt->all = 1;
                break;
            case 'O':
                opt->ordered = 1;
                break;
            case 'u':
                opt->ordered = 0;
                break;
            case 1:
                opt->genetree = 0;
                break;
			default:
				printf("Bad option: %c\n", c);
				exit(1);
				break;
		}
	}
	if(!thetaFileProvided && !stdinSet)
	{
		// now read thetas
		numThetas = argc - optind;
		thetas = (double *)calloc((size_t)numThetas, sizeof(double));
		CHECKPOINTER(thetas);
		thetaIdx = 0;
		for(i = optind; i < argc; i++)
		{
			//printf("argv[%i]: %s\n", i, argv[i]);
			success = (int)sscanf(argv[i], "%lf", &(thetas[thetaIdx]));
			if(!success)
				PERROR("Could not read theta in parsing options.");
			thetaIdx++;
		}
	}
	if(!migFileProvided && opt->numDemes == 2 && !stdinSet)
		PERROR("No migration file provided.\n");
    if((thetaFileProvided || migFileProvided) && stdinSet)
        PERROR("Cannot provide migration or theta file if reading input from STDIN (-s).");
    if(stdinSet && opt->numDemes == 1)
    {
        // thetas first, then locus
        char * line;
        line = (char *)malloc(sizeof(char) * DEFAULT_MAX_LINE_SIZE);
        CHECKPOINTER(line);

        int MAX_NUM_THETAS = 1000;
        thetas = (double *)malloc(sizeof(double) * MAX_NUM_THETAS);
        CHECKPOINTER(thetas);

        thetaIdx = 0;
        while(fgets(line, DEFAULT_MAX_LINE_SIZE, stdin) != NULL)
        {
            if(line[0] == '-')
                break;
            success = sscanf(line, "%lf\n", &(thetas[thetaIdx++]));
            if(!success)
                PERROR("Invalid theta in theta-file.");
        }

        numThetas = thetaIdx;

        thetas = (double *)realloc(thetas, sizeof(double)*numThetas);

        // locus
        opt->fin = stdin;
        free(line);
    }
    if(stdinSet && opt->numDemes == 2)
    {
        // thetas first, then migration rates, then locus
        char * line;
        line = (char *)malloc(sizeof(char) * DEFAULT_MAX_LINE_SIZE);
        CHECKPOINTER(line);

        int MAX_NUM_THETAS = 1000;
        thetas = (double *)malloc(sizeof(double) * MAX_NUM_THETAS);
        CHECKPOINTER(thetas);

        thetaIdx = 0;
        while(fgets(line, DEFAULT_MAX_LINE_SIZE, stdin) != NULL)
        {
            if(line[0] == '-')
                break;
            success = sscanf(line, "%lf\n", &(thetas[thetaIdx++]));
            if(!success)
                PERROR("Invalid theta in theta-file.");
        }

        numThetas = thetaIdx;

        thetas = (double *)realloc(thetas, sizeof(double)*numThetas);

        // migration rates
        int MAX_NUM_MIGRATES = 1000;

        migRates = (double *)malloc(sizeof(double) * MAX_NUM_MIGRATES);
        CHECKPOINTER(migRates);

        int migRateIdx = 0;
        while(fgets(line, DEFAULT_MAX_LINE_SIZE, stdin) != NULL)
        {
            if(line[0] == '-')
                break;
            success = sscanf(line, "%lf\n", &(migRates[migRateIdx++]));
            if(!success)
                PERROR("Invalid migration rate in input.");
        }

        numMigRates = migRateIdx;

        migRates = (double *)realloc(migRates, sizeof(double)*numMigRates);

        // locus
        opt->fin = stdin;
        free(line);
    }

	opt->thetas = thetas;
	opt->numThetas = numThetas;

	opt->migRates = migRates;
	opt->numMigRates = numMigRates;

	if(strlen(opt->filenameIn) == 0 && !stdinSet)
		PERROR("No input file specified");
    if(!stdinSet)
        opt->fin = fopen(opt->filenameIn, "r");
	if(opt->fin == NULL)
		PERROR("Cannot open input file");
	if(strlen(opt->filenameOut) == 0)
		opt->fout = stdout;
	else
	{
		opt->fout = fopen(opt->filenameOut, "w");
		if(opt->fout == NULL)
			PERROR("Cannot open output file");
	}
	if(opt->numDemes == -1)
		PERROR("Number of demes not specified (--numdemes, -D)");
	return;
}

int check_mono_D1(FILE * inp)
{
	char line[DEFAULT_MAX_LINE_SIZE];
	if(fgets(line, DEFAULT_MAX_LINE_SIZE, inp) == NULL)
		return -1;
	rewind(inp);
	unsigned char c;
	int numLines = 0;
	while(1)
	{
		c = (unsigned char)fgetc(inp);
		if(c == '1')
			return 0;
		if(c == '\n')
			numLines++;
		if(feof(inp))
			break;
	}
	//return numLines+1;    // can't have extra lines at the end.
	return numLines;    // normal text files end with a \n
}

void solve_D1(FILE * fin, int numThetas, double * thetas, int printAll, int ordered)
{
	int i, j, k;
	double theta, prob;
	BMat bmat;
	//int mono = check_mono_D1(fin);
    
    int numHaplotypes = -1;
    int mono = BMat_read_input(fin, &bmat, &numHaplotypes);

	if(!mono)
	{
		DataSet ds;
		DataSet_init(&ds, &bmat, numThetas, thetas, printAll, ordered);
		DataSet_free(&ds);
        free(ds.thetas);
		BMat_free(&bmat);
	}
	else // mono
	{
        BMat_free(&bmat);
        if(printAll)
        {
            for(i = numHaplotypes; i > 0; i--)
            {
                printf("%i; ", i);
                for(j = 0; j < numThetas-1; j++)
                {
                    theta = thetas[j];
                    prob = 1.0;
                    for(k = 1; k < i; k++)
                        prob *= (double)k / ((double)k + theta);
                        printf("%.16e ", prob);
                }
                prob = 1.0;
                theta = thetas[numThetas-1];
                for(k = 1; k < i; k++)
                    prob *= (double)k / ((double)k + theta);
                printf("%.16e\n", prob);
            }
        }
        else // !printAll
        {
            for(k = 0; k < numThetas; k++)
            {
                theta = thetas[k];
                prob = 1.0;
                for(i = numHaplotypes-1; i > 0; i--)
                    prob *= (double)i / ((double)i + theta);
                printf("%e\n", prob);
            }
        }
    }
	return;
}

void solve_D1_ctypes(char ** inp, const int numHaplotypes, const int numThetas, double * thetas, const int ordered, double * samplingProbs)
{
	int i, k;
	double theta, prob;
	BMat bmat;
	//int mono = check_mono_D1(fin);
    
    int mono = BMat_read_input_ctypes(inp, &bmat, numHaplotypes);

	if(!mono)
	{
		DataSet ds;
		DataSet_init_ctypes(&ds, &bmat, numThetas, thetas, ordered, samplingProbs); // 0 for printAll argument
		DataSet_free(&ds);
		BMat_free(&bmat);
	}
	else // mono
	{
        for(k = 0; k < numThetas; k++)
        {
            theta = thetas[k];
            prob = 1.0;
            for(i = numHaplotypes-1; i > 0; i--)
                prob *= (double)i / ((double)i + theta);
            samplingProbs[k] = prob;
        }
        BMat_free(&bmat);
    }
	return;
}

void solve_D2(FILE * fin, int numThetas, double * thetas, int numMigRates, double * migRates, int printAll, int ordered, int genetree)
{
	BMat2d b2;
	DataSet2d ds;
	int k;

    ////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////
    //// to make probabilities maximally compatible with *genetree*, which     /////
    //// defines theta as 4*N_{tot}*mu = 4*N*D*mu, which here is twice 4*N*mu. /////
    ////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////
    if(genetree)
    {
        for(k = 0; k < numThetas; k++)
            thetas[k] /= 2.0;
    }
	BMat2d_read_input(fin, &b2);
	DataSet2d_init(&ds, &b2, numThetas, thetas, numMigRates, migRates, printAll, ordered, genetree);
	BMat2d_free(&b2);
	free(thetas);
	free(migRates);
	return;
}

void test_ctypes(double * test)
{
    int i;
    for(i = 0; i < 3; i++)
        test[i] = 1.0;
    return;
}

void test_ctypes_2(char ** inp, int n)
{
    int i;

    for(i = 0; i < n; i++)
    {
        char * line = inp[i];
        printf("%s\n", line);
    }

    return;
}

// fills array of sampling probs for original reconfig only (will make another version that fills reconfig probs for all reconfigs)
void solve_D2_ctypes(char ** input, int * demes, int numHaplotypes, int numThetas, double * thetas, int numMigRates, double * migRates, int ordered, int genetree, double * samplingProbs)
{
	BMat2d b2;
	DataSet2d ds;
	int k;

    ////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////
    //// to make probabilities maximally compatible with *genetree*, which     /////
    //// defines theta as 4*N_{tot}*mu = 4*N*D*mu, which here is twice 4*N*mu. /////
    ////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////
    if(genetree)
    {
        for(k = 0; k < numThetas; k++)
            thetas[k] /= 2.0;
    }

	BMat2d_read_input_ctypes(input, demes, numHaplotypes, &b2);
	DataSet2d_solve_ctypes(&ds, &b2, numThetas, thetas, numMigRates, migRates, ordered, genetree, samplingProbs);
	BMat2d_free(&b2);
	return;
}

void solve_D2_ctypes_all(char ** input, int * demes, int numHaplotypes, int numThetas, double * thetas, int numMigRates, double * migRates, int ordered, int genetree, int ** recIdxs, double ** samplingProbs)
{
	BMat2d b2;
	DataSet2d ds;
	int k;

    ////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////
    //// to make probabilities maximally compatible with *genetree*, which     /////
    //// defines theta as 4*N_{tot}*mu = 4*N*D*mu, which here is twice 4*N*mu. /////
    ////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////
    if(genetree)
    {
        for(k = 0; k < numThetas; k++)
            thetas[k] /= 2.0;
    }

	BMat2d_read_input_ctypes(input, demes, numHaplotypes, &b2);
	DataSet2d_solve_ctypes_all(&ds, &b2, numThetas, thetas, numMigRates, migRates, ordered, genetree, recIdxs, samplingProbs);
	BMat2d_free(&b2);
	return;
}

void test_solve_D2_ctypes()
{
    int i;
    char ** seqs = (char **)malloc(sizeof(char *) * 8);
    seqs[0] = (char *)malloc(sizeof(char) * 10);
    seqs[1] = (char *)malloc(sizeof(char) * 10);
    seqs[2] = (char *)malloc(sizeof(char) * 10);
    seqs[3] = (char *)malloc(sizeof(char) * 10);
    seqs[4] = (char *)malloc(sizeof(char) * 10);
    seqs[5] = (char *)malloc(sizeof(char) * 10);
    seqs[6] = (char *)malloc(sizeof(char) * 10);
    seqs[7] = (char *)malloc(sizeof(char) * 10);

    strcpy(seqs[0], "000 0\n");
    strcpy(seqs[1], "000 0\n");
    strcpy(seqs[2], "100 0\n");
    strcpy(seqs[3], "100 0\n");
    strcpy(seqs[4], "010 1\n");
    strcpy(seqs[5], "010 1\n");
    strcpy(seqs[6], "011 1\n");
    strcpy(seqs[7], "011 1\n");
    double * sp = (double *)malloc(sizeof(double) * 8);
    double thetas[2] = {1.0, 2.0};
    double migRates[2] = {3.0, 4.0};

    int demes[] = {0,0,0,0,1,1,1,1};

    solve_D2_ctypes(seqs, (int *)demes, 8, 2, thetas, 2, migRates, 1, 1, sp);

    for(i = 0; i < 4; i++)
        printf("%e\n", sp[i]);
    return;
}

void test_solve_D2_ctypes_all()
{
    int i, j;
    char ** seqs = (char **)malloc(sizeof(char *) * 4);
    seqs[0] = (char *)malloc(sizeof(char) * 10);
    seqs[1] = (char *)malloc(sizeof(char) * 10);
    seqs[2] = (char *)malloc(sizeof(char) * 10);
    seqs[3] = (char *)malloc(sizeof(char) * 10);

    strcpy(seqs[0], "100\n");
    strcpy(seqs[1], "101\n");
    strcpy(seqs[2], "100\n");
    strcpy(seqs[3], "010\n");
    double thetas[2] = {1.0, 2.0};
    double migRates[2] = {3.0, 4.0};

    int demes[] = {0,0,1,1};

    int numReconfigs = 20;

    int recIdxs[20][8] = {0};
    double sampProbs[20][8] = {0.0};

    double * sp = (double *)malloc(sizeof(double) * 8);
    solve_D2_ctypes(seqs, (int *)demes, 4, 2, thetas, 2, migRates, 1, 1, sp);

    for(i = 0; i < 4; i++)
        printf("%f ", sp[i]);
    printf("\n");

    solve_D2_ctypes_all(seqs, (int *)demes, 4, 2, (double *)thetas, 2, (double *)migRates, 1, 1, (int **)recIdxs, (double **)sampProbs);

    for(i = 0; i < numReconfigs; i++)
    {
        for(j = 0; j < 8; j++)
            printf("%i ", recIdxs[i][j]);
        printf(": ");
        for(j = 0; j < 4; j++)
            printf("%f ", sampProbs[i][j]);
        printf("\n");
    }
    return;
}

void SolverOptions_run_program(SolverOptions * opt)
{
	if(opt->numDemes == 2)
		solve_D2(opt->fin, opt->numThetas, opt->thetas, opt->numMigRates, opt->migRates, opt->all, opt->ordered, opt->genetree);
	else if(opt->numDemes == 1)
		solve_D1(opt->fin, opt->numThetas, opt->thetas, opt->all, opt->ordered);
	return;
}

int main(int argc, char ** argv)
{
    test_solve_D2_ctypes_all();
    /*
	SolverOptions opt;
	SolverOptions_parse_options(argc, argv, &opt);
	SolverOptions_run_program(&opt);
    */
	return 0;
}
