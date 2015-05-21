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
	int32_t numDemes;
	double * thetas;
	int32_t numThetas;
	double * migRates;
	int32_t numMigRates;
	FILE * fin;
	FILE * fout;
    int32_t all;
    int32_t ordered;
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
	{0, 0, 0, 0}
};

void Solver_print_usage()
{
	fprintf(stderr, "\nUsage: solver [OPTIONS] -i datafile\n");
	return;
}

void SolverOptions_parse_options(int32_t argc, char ** argv, SolverOptions * opt)
{
	int32_t c, i, optionIndex, success, numThetas = -1, numMigRates = -1, thetaIdx, migrationIdx, thetaFileProvided = 0, migFileProvided = 0, stdinSet = 0;
	char ch, * line;
	FILE * thetain, * migrationin;
	double * thetas = NULL, * migRates = NULL;
	opt->numDemes = -1;
    opt->ordered = 0;
    opt->all = 0;
	strcpy(opt->filenameIn, "");
	strcpy(opt->filenameOut, "");
	//opt->theta = -1.0;
	while(1)
	{
		c = getopt_long(argc, argv, "hi:D:t:M:o:sa", long_options, &optionIndex);
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

        int32_t MAX_NUM_THETAS = 1000;
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

        int32_t MAX_NUM_THETAS = 1000;
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
        int32_t MAX_NUM_MIGRATES = 1000;

        migRates = (double *)malloc(sizeof(double) * MAX_NUM_MIGRATES);
        CHECKPOINTER(migRates);

        int32_t migRateIdx = 0;
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

int32_t check_mono_D1(FILE * inp)
{
	char line[DEFAULT_MAX_LINE_SIZE];
	if(fgets(line, DEFAULT_MAX_LINE_SIZE, inp) == NULL)
		return -1;
	rewind(inp);
	unsigned char c;
	int32_t numLines = 0;
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

void solve_D1(FILE * fin, int32_t numThetas, double * thetas, int32_t ordered)
{
	int32_t i, j, k;
	double theta, prob;
	BMat bmat;
	//int32_t mono = check_mono_D1(fin);
    
    int32_t numHaplotypes = -1;
    int32_t mono = BMat_read_input(fin, &bmat, &numHaplotypes);

	if(!mono)
	{
		DataSet ds;
		DataSet_init(&ds, &bmat, numThetas, thetas, ordered);
		DataSet_free(&ds);
		BMat_free(&bmat);
	}
	else // mono
	{
        BMat_free(&bmat);
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
	
	return;
}

void solve_D2(FILE * fin, int32_t numThetas, double * thetas, int32_t numMigRates, double * migRates, int32_t printAll, int32_t ordered)
{
	BMat2d b2;
	DataSet2d ds;
	int32_t k;

    ////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////
    //// to make probabilities maximally compatible with *genetree*, which     /////
    //// defines theta as 4*N_{tot}*mu = 4*N*D*mu, which here is twice 4*N*mu. /////
    ////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////
	for(k = 0; k < numThetas; k++)
		thetas[k] /= 2.0;
	BMat2d_read_input(fin, &b2);
	DataSet2d_init(&ds, &b2, numThetas, thetas, numMigRates, migRates, printAll, ordered);
	BMat2d_free(&b2);
	free(thetas);
	free(migRates);
	return;
}

void SolverOptions_run_program(SolverOptions * opt)
{
	if(opt->numDemes == 2)
		solve_D2(opt->fin, opt->numThetas, opt->thetas, opt->numMigRates, opt->migRates, opt->all, opt->ordered);
	else if(opt->numDemes == 1)
		solve_D1(opt->fin, opt->numThetas, opt->thetas, opt->ordered);
	return;
}

int main(int argc, char ** argv)
{
	SolverOptions opt;
	SolverOptions_parse_options(argc, argv, &opt);
	SolverOptions_run_program(&opt);
	return 0;
}
