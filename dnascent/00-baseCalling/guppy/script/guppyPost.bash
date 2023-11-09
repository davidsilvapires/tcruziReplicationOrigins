#!/usr/bin/env bash

# Making symbolic links for the final result.
for i in output/pass/*
do
  ln -s ../${i} final/
done

ln -s ../output/sequencing_summary.txt final/

exit 0

# Making a symbolic link for sequencing_summary.

