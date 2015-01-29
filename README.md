Infinite-sites Likelihood Calculator
====================================

This program calculates sampling probabilities under the infinite-sites model
of mutation for one and two demes. The method used is that of
[Wu (2010)](http://ieeexplore.ieee.org/xpl/articleDetails.jsp?arnumber=5383348).

Installation
----------

Running `make` in this directory will create the `solver` binary. By default,
the program uses UMFPACK from the
[SuiteSparse](http://faculty.cse.tamu.edu/davis/suitesparse.html) library for
matrix factorization. I *think* it still supports the CSPARSE routines
(included) if you have trouble installing SuiteSparse. `make` can be run with
various different targets; see Makefile.

Usage
-----

Single-deme usage:
	
	solver --numdemes 1 --theta 0.5 -i example_seq_D1.txt

Two-deme usage:

	solver --numdemes 2 --theta 0.5 --migrates "1.0 1.0" -i example_seq_D2.txt

The definitions of theta = 4* N * mu and M = 4 * N * m should correspond to
those of Wu (2010). 
