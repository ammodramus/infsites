Infinite-sites Likelihood Calculator
====================================

This program calculates sampling probabilities under the infinite-sites
mutation model for one and two demes. The program uses the dynamic programming
approach of [Wu
(2010)](http://ieeexplore.ieee.org/xpl/articleDetails.jsp?arnumber=5383348),
essentially solving the back-in-time recursion of [Griffiths and Tavar\'e
(1994)](https://projecteuclid.org/euclid.ss/1177010378) forward in time,
starting from the base case.

Installation
----------

Running `make` in this directory will create the `solver` binary. `make` can be
run with various different targets; see Makefile. This software uses the
`[csparse](https://people.sc.fsu.edu/~jburkardt/c_src/csparse/csparse.html)`
library (included) for inverting matrices.

Usage
-----

Single-deme usage:

    solver -i example_seq_D1.txt -t example_t -D 1

Two-deme usage:

	solver -i example_seq_D2.txt -t example_t -M example_M -D 2

The definitions of theta = 4* N * mu and M = 4 * N * m should correspond to
those of Wu (2010). 
