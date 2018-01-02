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

### Running
<code>
./reordering-library -m &lt;path of .mtx file&gt; -a &lt;algorithm&gt; -p &lt;number of threads&gt; -b &lt;percent of chunk&gt;
</code>

<br>

<code>
  &lt;algorithm&gt;
</code>

 * 	0: Serial RCM
 * 	1: Serial Sloan
 * 	2: HSL RCM
 * 	3: HSL Spectral
 * 	4: HSL Sloan
 * 	5: Unordered RCM
 * 	6: Leveled RCM
 * 	7: Bucket RCM
 * 	8: Relaxed Order Sloan
 * 	9: Boost RCM
 * 10: Boost Sloan
 * 11: Logical Bag Sloan
 * 12: Shrinked RCM

<code>
  &lt;percent of chunk&gt;
</code>

 * It is recommended the value of 0.5.
 
Example:

<code>
  ./reordering-library -m ./Matrices/rail_5177.mtx -a 5 -p 4 -b .5
</code>

* In this example, the matrix rail_5177 is processed by the Unordered RCM algorithm. It is executed with 4 threads.
 
### Profiling
<code>
valgrind --leak-check=yes myprog arg1 arg2
</code>

* For memory check, Valgrind has been used.

<code>
valgrind --tool=callgrind program [program_options]
</code>

* For performance profiling, Callgrind has been used.

<code>
kcachegrind callgrind.out.XXX
</code>

* For graphical performance visualization, KCachegrind has been used. The file callgrind.out.XXX is yielded by Callgrind, and XXX is the process identifier. 
