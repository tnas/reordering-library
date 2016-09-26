#!/bin/sh
set xrange [0:10]
set yrange [0:10] reverse
set size square
set terminal png
set output "original_matrix.png"
plot "original_matrix.mtx" using 1:2 notitle with dots linecolor rgb "#000000"
