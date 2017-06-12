#!/bin/bash

# Run this script following one of the two options:
# Option 1: . reordering.sh
# Option 2: source reordering.sh

# Disabling dynamic adjustment of the number of threads
# OPTIONS: true, false
export OMP_DYNAMIC="false" 

# Make waiting threads spin, consuming processor cycles while waiting
# OPTIONS: active, passive
export OMP_WAIT_POLICY="active"

# List of places. (<stride> default is 1)
# {<lower-bound>:<length>:<stride>},...,{<lower-bound>:<length>:<stride>}
export OMP_PLACES="{0:8},{8:8}"

# Defining threads affinity
# OPTIONS: true, false, master, spread, close
export OMP_PROC_BIND="close"

#./build/reordering-library -t 1

#nohup ./reordering-library -t 1 &
