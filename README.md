# plum
A sage package for analyzing plumbed manifolds.

To use:

0. You will need [sagemath](https://www.sagemath.org/) installed.

1. Download the file plum.sage.

2. Open sage and type: `load('replace this with the file path of plum.sage')`

For example, if the path is user/downloads/plum.sage, you would type: `load('user/downloads/plum.sage')`

3. Once you have loaded the script, you can create a new plumbing. For example, the plumbing consisting of 2 vertices joined by a single
    edge such that one vertex has weight -3 and the other has weight -10 can be created as follows: `P = Plumbing({0:-3,1:-10}, [(0,1)])`

4. There are several properties and methods of P. For example, if you want to:
- get the intersection form of P, type: `P.interesection_form`
- compute the first n levels of zhat for a spinc representative a (using GM notation and represented as a list), type: `P.zhat(a, n)`
- compute the first n levels of the bigraded root for a spinc representative k (using Nemethi notation and represented as a list), type: `P.bigraded_root(k, n)`
