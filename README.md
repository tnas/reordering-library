## Parallelizations of Irregular Algorithms for Sparse Matrices Reordering

This project aims the parallelization of some algorithms for the bandwidth and wavefront minimization problems. 
The related algorithms are:
* Reverse Cuthill McKee - RCM (Bandwidth Reduction)
* Sloan (Wavefront Reduction)

Three parallelized algorithms were builded from the serial RCM:
* Unordered RCM
* Bucket RCM
* Shrinked RCM

About the serial Sloan, two parallel algorithm were implemented:
* Bag Sloan
* Modified Sloan

All algorithms have been implemented using the OpenMP parallel platform. The used programming language was C.

[![DOI](https://www.zenodo.org/badge/82069776.svg)](https://www.zenodo.org/badge/latestdoi/82069776)
