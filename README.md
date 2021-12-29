## `TurboRDM-ED`

Minimal implementation of a fast-trace algorithm to compute reduced density matrices within quantum (cluster) impurity problems. Based on the [QcmPlab](https://github.com/QcmPlab) exact-diagonalization solvers, whose Fock-space structure is assumed. More specifically here we assume Sz-conservation, implemented in a generalized Lin-Gubernatis scheme. More info available [on the arXiv](https://arxiv.org/abs/2105.06806).

The algorithm is fully included on the master branch of [QcmPlab/CDMFT-LANC-ED](https://github.com/QcmPlab/CDMFT-LANC-ED) since commit [ed4ba2d](https://github.com/QcmPlab/CDMFT-LANC-ED/commit/ed4ba2d798d31588da7675565e2dd0c1ce821744).

#### Dependencies

The code relies on the [SciFortran](https://github.com/QcmPlab/SciFortran) library.  
The dependency list included therein applies (except MPI and Scalapack that here are totally irrelevant).

#### Code structure

- The `COMMON_VARS` module defines all the global types and methods, including the dynamical, sparse data structures needed to build a bath-impurity separated Fock-space map, suitable to boost bath-tracing speed.
- The `AUX_FUNX` module implements the Sz-sector building routine, i.e. the core of the present implementation, together with many auxiliary procedures (some of which imported verbatim from the QcmPlab ED libraries). A special mention goes to the subtracing procedures, that exploit a further bitwise factorizazion of the states, so to allow a site-by-site reduction of the impurity density matrix in the cluster case; the tracing therein does not implement the fast algorithm since the bit-dimension of the impurity states does not justify further loading of the data structures and the direct recipe has an acceptable performance.
- The `main` driver performs several tests on both the main fast-trace over the bath states and the site-by-site subtrace. It parses some input variables and flags so to build a suitable _fictitious_ eigenspace, by generating some random/uniform eigenvectors. Then it calls the map-building procedure and feeds everything to the slow- and fast-bath-tracing algorithms and to the subtracing routine. Finally it compares all the resulting matrices, together with the required cpu times.

#### Input file

The makefile compiles the program, putting the executable in the `test/` directory. There you can find an editable input-file, where you can define the dimension of the impurity problem and set some flags for the tests.

```
NLAT  : number of sites in the cluster; 1 corresponds to the standard Anderson Impurity Model (AIM).
NORB  : number of orbitals per site.
NBATH : number of bath levels per site, per orbital (*)
NUP   : the [Nup,Ndw] tuple defines the given sector by means of the Sz = Nup - Ndw identity;
NDW     do not forget the obvious constraint: Ns = Nup + Ndw = Nlat * Norb * (Nbath + 1).

VERBOSE : if T prints a lot (loads!) of useful information on screen.
FAST    : if T skips the slow algorithm for the bath trace.
RANDOM  : if T gives random eigenvectors; else all the sector will contribute with a uniform weight.

(*) Think of the bath as a set of noninteracting replicas of the cluster.
```


#### Copyright

© 2021 - Adriano Amaricci, Gabriele Bellomia.  
All rights reserved.

The software is provided with no license, as such it is protected by copyright. The software is provided as it is and can be read and copied, in agreement with the Terms of Service of GITHUB. Use of the code is constrained to authors agreement.
