#!/usr/bin/env bash

# Verifing the total number of files
# ls output/ | wc -l

# Verifing if any file is empty (Probably due to a compression in step 02)
# find output/ -size 0 | wc -l

# Deleting empty files in output (in other directories the symbolic link has a size that I could not
# difference from not empty files -- check if its possible to delete from the symbolic links).

# find output/ -size 0 -delete

# Making symbolic links for the final results.

for i in output/*.brduDetect
do 
  ln -s ../${i} final/
done

for FILE in output/*.regions
do
  ln -s ../${FILE} final/
done

for FILE in output/*.forkSense
do
  ln -s ../${FILE} final/
done
