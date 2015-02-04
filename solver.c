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
	double migRates[2];
	FILE * fin;
	FILE * fout;
} SolverOptions;

static struct option long_options[] =
{
	{"help", no_argument, 0, 'h'},
	{"input", required_argument, 0, 'i'},
	{"numdemes", required_argument, 0, 'D'},
	{"theta",  required_argument, 0, 't'},
	{"migrate", required_argument, 0, 'M'},
	{"theta-file", required_argument, 0, 1},
	{"migration-file", required_argument, 0, 2},
	{"output", required_argument, 0, 'o'},
	{0, 0, 0, 0}
};

void Solver_print_usage()
{
	fprintf(stderr, "\nUsage: solver [OPTIONS] -i datafile\n");
	return;
}

void SolverOptions_parse_options(int32_t argc, char ** argv, SolverOptions * opt)
{
	int c, i, optionIndex, success;
	opt->migRates[0] = opt->migRates[1] = -1.0;
	opt->numDemes = -1;
	strcpy(opt->filenameIn, "");
	strcpy(opt->filenameOut, "");
	//opt->theta = -1.0;
	while(1)
	{
		c = getopt_long(argc, argv, "hi:D:t:M:o:", long_options, &optionIndex);
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
			/*
			case 't':
				success = (int)sscanf(optarg, "%lf", &(opt->theta));
				if(success == 0)
					PERROR("Invalid input for theta (--theta, -t)");
				break;
			*/
			case 'M':
				success = (int)sscanf(optarg, "%lf %lf", &(opt->migRates[0]), &(opt->migRates[1]));
				if(success != 2)
					PERROR("Invalid input for migration rates (--migrates, -M)");
				break;
			case 'o':
				strncpy(opt->filenameOut, optarg, sizeof(opt->filenameOut));
				break;
			case 1:
				// deal with theta file
				break;
			case 2:
				// deal with migration-rate file
				break;
			default:
				abort();
				break;
		}
	}
	// now read thetas
	int32_t numThetas = argc - optind;
	double * thetas = (double *)calloc((size_t)numThetas, sizeof(double));
	CHECKPOINTER(thetas);
	opt->thetas = thetas;
	opt->numThetas = numThetas;
	int32_t thetaIdx = 0;
	for(i = optind; i < argc; i++)
	{
		//printf("argv[%i]: %s\n", i, argv[i]);
		success = (int)sscanf(argv[i], "%lf", &(opt->thetas[thetaIdx]));
		if(!success)
			PERROR("Could not read theta in parsing options.");
		thetaIdx++;
	}

	if(strlen(opt->filenameIn) == 0)
		PERROR("No input file specified");
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
	if(opt->numDemes == 2 && (opt->migRates[0] == -1 || opt->migRates[1] == -1))
		PERROR("Migration rate not specified (--migrates, -M)");
	//if(opt->theta == -1.0)
		//PERROR("Mutation rate not specified (--theta, -t)");
	//free(thetas);
	return;
}

int32_t check_mono_D1(FILE * inp)
{
	char line[DEFAULT_MAX_LINE_SIZE];
	if(fgets(line, DEFAULT_MAX_LINE_SIZE, inp) == NULL)
		return -1;
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
	return numLines+1;    // can't have extra lines at the end.
}

void solve_D1(FILE * fin, int32_t numThetas, double * thetas)
{
	int32_t i, k;
	double theta, prob;
	BMat bmat;
	int32_t mono = check_mono_D1(fin);
	if(!mono)
	{
		fseek(fin, 0, SEEK_SET);
		BMat_read_input(fin, &bmat);
		DataSet ds;
		DataSet_init(&ds, &bmat, numThetas, thetas);
		DataSet_free(&ds);
		BMat_free(&bmat);
	}
	// have to make this deal with multiple thetas as well.
	else
	{
		for(k = 0; k < numThetas; k++)
		{
			theta = thetas[k];
			prob = 1.0;
			for(i = mono-1; i > 0; i--)
				prob *= (double)i / ((double)i + theta);
			printf("%e\n", prob);
		}
	}
	return;
}

//solve_D2(opt->fin, opt->numThetas, opt->thetas, opt->migRates);
void solve_D2(FILE * fin, int32_t numThetas, double * thetas, double * migRates)
{
	BMat2d b2;
	DataSet2d ds;
	int32_t k;
	//double theta = 1.0;
	//double * migRates = (double *)malloc(sizeof(double) * 2);
	//double migRates[2];
	//migRates[0] = 0.5;
	//migRates[1] = 0.5;
	for(k = 0; k < numThetas; k++)
		thetas[k] /= 2.0; 		// to make probabilities maximally compatible with genetree, which defines theta as 4*N_{tot}*mu = 4*N*D*mu, which here is twice 4*N*mu.
	//BMat2d_read_input("testfile2d", &b2);
	BMat2d_read_input(fin, &b2);
	DataSet2d_init(&ds, &b2, numThetas, thetas, migRates);
	BMat2d_free(&b2);
	return;
}

void SolverOptions_run_program(SolverOptions * opt)
{
	if(opt->numDemes == 2)
		solve_D2(opt->fin, opt->numThetas, opt->thetas, opt->migRates);
	else if(opt->numDemes == 1)
		solve_D1(opt->fin, opt->numThetas, opt->thetas);
	return;
}

void SolverOptions_report_options(SolverOptions * opt)
{
	int32_t k;
	REPORTI(opt->numDemes);
	//REPORTF(opt->theta);
	for(k = 0; k < opt->numThetas; k++)
		REPORTF(opt->thetas[k]);
	REPORTF(opt->migRates[0]);
	REPORTF(opt->migRates[1]);
	printf("input: %s\n", opt->filenameIn);
	printf("output: %s\n", opt->filenameOut);
	return;
}

int main(int argc, char ** argv)
{
	SolverOptions opt;
	SolverOptions_parse_options(argc, argv, &opt);
	SolverOptions_run_program(&opt);
	return 0;
}
