#!/bin/sh
set xrange [0:17880]
set yrange [0:17880] reverse
set size square
set terminal png
set output "reordered_matrix.png"
plot "reordered_matrix.mtx" using 1:2 notitle with dots linecolor rgb "#000000"