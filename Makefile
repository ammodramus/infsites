all: optimize

debug:
	gcc solver.c dataset2d.c datconfig.c dataset.c bmat.c bmat2d.c nodelist.c datconfig2d.c node.c hash.c matrix.c murmur3.c -lm -g -Wall -o solver_debug

optimize:
	gcc solver.c dataset2d.c datconfig.c dataset.c bmat.c bmat2d.c nodelist.c datconfig2d.c node.c hash.c matrix.c murmur3.c -lm -O3 -march=native -Wall -o solver

profile:
	gcc solver.c dataset2d.c dataset.c bmat.c bmat2d.c nodelist.c datconfig2d.c datconfig.c node.c hash.c matrix.c murmur3.c -lm -g -Wall -o solver_profile

static:
	gcc solver.c dataset2d.c dataset.c bmat.c bmat2d.c nodelist.c datconfig2d.c datconfig.c node.c hash.c matrix.c murmur3.c -lpthread -lm -O3 -march=native --static -o solver_static

dll:
	gcc -shared -o libsolver.so -fPIC solver.c dataset2d.c datconfig.c dataset.c bmat.c bmat2d.c nodelist.c datconfig2d.c node.c hash.c matrix.c murmur3.c -lm -lumfpack -lamd -lcholmod -lblas -O3 -march=native

dll_debug:
	gcc -shared -o libsolver.so -fPIC solver.c dataset2d.c datconfig.c dataset.c bmat.c bmat2d.c nodelist.c datconfig2d.c node.c hash.c matrix.c murmur3.c -lm -lumfpack -lamd -lcholmod -lblas -g
