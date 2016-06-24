## Parallelization of Reordering Algorithms for Bandwidth and Wavefront Reduction of Large Sparse Matrices

This project aims the parallelization of some algorithms for the bandwidth and wavefront reduction problems. 
The related algorithms are:
* Reverse Cuthill McKee - RCM (Bandwidth Reduction)
* Sloan (Wavefront Reduction)

Three parallelized algorithms were builded from the serial RCM:
* Unordered RCM
* Leveled RCM
* Bucket RCM

About the serial Sloan, just one parallel algorithm was implemented. It uses the concept of "Logical Bags" as main data structure.

All algorithms have been implemented using the OpenMP parallel platform. Moreover, the programming language used is C.
