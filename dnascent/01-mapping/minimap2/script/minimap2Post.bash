#!/usr/bin/env bash

# Making symbolic links for the final result.
for i in output/*
do
  ln -s ../${i} final/
done

exit 0

