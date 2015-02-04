all: optimize

debug:
	gcc solver.c dataset2d.c datconfig.c dataset.c bmat.c bmat2d.c nodelist.c datconfig2d.c node.c hash.c matrix.c murmur3.c -lm -lumfpack -lamd -lcholmod -lblas -g -Wall -o solver_debug

optimize:
	gcc solver.c dataset2d.c datconfig.c dataset.c bmat.c bmat2d.c nodelist.c datconfig2d.c node.c hash.c matrix.c murmur3.c -lm -lumfpack -lamd -lcholmod -lblas -O3 -march=native -Wall -o solver

profile:
	gcc solver.c dataset2d.c dataset.c bmat.c bmat2d.c nodelist.c datconfig2d.c datconfig.c node.c hash.c matrix.c murmur3.c -lm -lumfpack -lamd -lcholmod -lblas -g -Wall -o solver_profile

static:
	gcc solver.c dataset2d.c dataset.c bmat.c bmat2d.c nodelist.c datconfig2d.c datconfig.c node.c hash.c matrix.c murmur3.c -lumfpack -lamd -lcholmod -lblas -lpthread -lm -O3 -msse3 --static -o solver_static
static2:
	gcc solver.c dataset2d.c dataset.c bmat.c bmat2d.c nodelist.c datconfig2d.c datconfig.c node.c hash.c matrix.c murmur3.c -lumfpack -lamd -lcholmod -lblas -lpthread -lm -O2 --static -o solver_static2
static3:
	gcc solver.c dataset2d.c dataset.c bmat.c bmat2d.c nodelist.c datconfig2d.c datconfig.c node.c hash.c matrix.c murmur3.c -lumfpack -lamd -lcholmod -lblas -lpthread -lm -g --static -o solver_static3

csparse:
	gcc -DCSPARSECOMPILE solver.c dataset2d.c datconfig.c dataset.c bmat.c bmat2d.c nodelist.c datconfig2d.c node.c hash.c matrix2.c murmur3.c -lm -g -Wall -o solver_csparse
csparseoptimize:
	gcc -DCSPARSECOMPILE solver.c dataset2d.c datconfig.c dataset.c bmat.c bmat2d.c nodelist.c datconfig2d.c node.c hash.c matrix2.c murmur3.c -lm -O3 -march=native -Wall -o solver_csparse
csparsestatic:
	gcc -DCSPARSECOMPILE solver.c dataset2d.c datconfig.c dataset.c bmat.c bmat2d.c nodelist.c datconfig2d.c node.c hash.c matrix2.c murmur3.c -lm -O3 -msse3 --static -Wall -o solver_csparse_static

