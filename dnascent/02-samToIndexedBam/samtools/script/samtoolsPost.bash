#!/usr/bin/env bash

# Making symbolic links for the final result.

for i in output/*sam.sortedBam
do
  NEWNAME=`basename ${i} | sed 's/sam.sortedBam/bam/'`
  ln -s ../${i} final/${NEWNAME}
done

for i in output/*sam.sortedBam.bai
do
  NEWNAME=`basename ${i} | sed 's/sam.sortedBam.bai/bam.bai/'`
  ln -s ../${i} final/${NEWNAME}
done

exit 0
