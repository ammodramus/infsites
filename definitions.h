#ifndef DEFINITIONS_HEADER
#define DEFINITIONS_HEADER

#include <stdio.h>
#include <stdint.h>

#define PERROR(msg,...) do {fprintf(stderr,"\n\nProgram error:\n%s\n\n",msg); exit(-1);} while(0)
#define MAX(a,b) (a > b ? a : b)
#define MIN(a,b) (a > b ? b : a)
#define CHECKPOINTER(ptr) if(ptr == NULL) do {fprintf(stderr,"\nFailed memory allocation: %s\n",#ptr); exit(-1);} while(0)
#define NOTININTERVAL(x, a, b) (x < a || x > b)
#define REPORTI(x) printf(#x " = %i\n", x)
#define REPORTF(x) printf(#x " = %f\n", x)

#define INT32T_MAX 2147483647 

typedef int twoints[2];

#endif
