#!/usr/bin/env bash

# Directory where we can find the data for this run.
DATA_DIR="../../metadata"

# File that contains metadata about this project.
METADATA="${DATA_DIR}/samples.txt"

# Computing how many samples we have for this run.
NUM_SAMPLES=`wc -l ${METADATA} | cut -d' ' -f1`
let "NUM_SAMPLES = NUM_SAMPLES - 1" # dummySample doesn't count as a valid input.

# Making a symbolic link for the final result.
for i in `seq 1 ${NUM_SAMPLES}`; do
   ln -s ../output/${i}-filteredSorted.bam final/${i}.bam
   ln -s ../output/${i}-filteredSorted.bam.bai final/${i}.bam.bai
done

exit 0
