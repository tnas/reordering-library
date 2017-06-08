## Non-Specultive Data-Driven Parallelizations of Irregular Algorithms for Sparse Matrices Reordering

This project aims the parallelization of some algorithms for the bandwidth and wavefront reduction problems. 
The related algorithms are:
* Reverse Cuthill McKee - RCM (Bandwidth Reduction)
* Sloan (Wavefront Reduction)

### Dependencies
On linux, execute the command to install libraries used by the program:

<code>
sudo apt-get install cmake libboost-all-dev gfortran libblas-dev
</code>

### Profiling
* For memory check, Valgrind has been used.
<code>
valgrind --leak-check=yes myprog arg1 arg2
</code>
<br/>
* For performance profiling, Callgrind/KCachegrind has been used.
<code>
valgrind --tool=callgrind program [program_options]
kcachegrind callgrind.out.XXX
</code>
