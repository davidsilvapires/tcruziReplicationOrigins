#!/usr/bin/env bash

# Setting up the initial directories structures.
mkdir final input job log output

# Making symbolic links for input data.
for FILE in `ls ../04-highMovingForkProbability/brdu/final/`
do
  ln -s ../../../../../04-highMovingForkProbability/awk/nonSinc/1/final/${FILE} input/${FILE:1}
done

# Writing the script.

for FILE in `ls input/`
do
  head -n 1 input/${FILE} > input/filtered.${FILE} 
  awk '$2>=0.7 || $3>=0.7 {print $0}' <(tail -n +2 input/${FILE}) >> input/filtered.${FILE}
done

exit 0

