# Plum
A sage module for analyzing plumbed manifolds.

To use:

0. You will need [sagemath](https://www.sagemath.org/) installed.

1. Download the file plum.sage from the latest
[release](https://github.com/peterkj1/plum/releases/tag/v1.0). You do not need
to download plum.py. That file is only used for generating sphinx documentation.

2. Open sage and type: `load('replace this with the file path of plum.sage')`

For example, if the path is user/downloads/plum.sage, you would
type: `load('user/downloads/plum.sage')`

3. Once you have loaded the script, you can create a new plumbing. For example,\  
`P = Plumbing({0:-1, 1:-2, 2:-3, 3:-7}, [(0,1), (0,2), (0,3)])`  

is the plumbing consisting of 4 vertices {0,1,2,3} with weights
-1, -2, -3, and -7 respectively. Vertices 0 and 1, 0 and 2, 0 and 3
are connected by edges.


For more details, see the [documentation](https://peterkj1.github.io/plum/index.html).
