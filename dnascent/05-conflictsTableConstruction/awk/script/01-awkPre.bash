#!/usr/bin/env bash

if [ $# -ne 1 ]
then
  echo "Usage: `basename $0` <FORK_PROBABILITY>"
  exit 1
fi

# Building conflicts table

numCompare() {
   awk -v n1="$1" -v n2="$2" 'BEGIN {printf (n1>=n2?"1":"0")}'
}

# FILE="output/filtered.0015aa63-71ef-45b4-8195-d26fca48029b_TcChr28-S_195683_224880_fwd.txt"
FILE="/tmp/2"

# Header definition.
HEADER=`echo "Read Name\tChrom\tMapped Read Start (bp)\tMapped Read End (bp)\tMapped Read Size (bp)\tStrand\tFork Start (bp)\tFork End (bp)\tFork Size (bp)\tFork direction (+/-)"`
echo -e ${HEADER}

# Read coordinates extraction.
READ_COORDINATES=`awk '{if(NR==1) print $0}' ${FILE} | cut -d' ' -f1-5 | awk -v OFS='\t' '{ print $1, $2, $3, $4, $4 - $3, $5 }'`

START=""
FORK_DIRECTION=""
PREVIOUS_COORDINATE=""

IFS='
'
for LINE in `tail -n +2 ${FILE}`; do
  COORDINATE=`echo ${LINE} | cut -f1`
  LEFT=`echo ${LINE} | cut -f2`
  if [ `numCompare ${LEFT} $1` -eq 1 ]; then
    if [ -z ${START} ]
    then
      START=${COORDINATE}
    fi
    FORK_DIRECTION="LEFT"
    PREVIOUS_COORDINATE=${COORDINATE}
    continue
  else
    if [ ! -z ${START} ]
    then
      END=${PREVIOUS_COORDINATE}
      let "DIFF = END - START"
      echo -e "${READ_COORDINATES}\t${START}\t${END}\t${DIFF}\t${FORK_DIRECTION}"
      START=""
      FORK_DIRECTION=""
    fi
  fi
done

# Dealing with the fucking last line.
if [ ! -z ${START} ]
then
  END=${PREVIOUS_COORDINATE}
  let "DIFF = END - START"
  echo -e "${READ_COORDINATES}\t${START}\t${END}\t${DIFF}\t${FORK_DIRECTION}"
fi


START=""
FORK_DIRECTION=""
PREVIOUS_COORDINATE=""

# Dealing with the fucking first line.

for LINE in `tail -n +2 ${FILE}`; do
  COORDINATE=`echo ${LINE} | cut -f1`
  RIGHT=`echo ${LINE} | cut -f3`
  if [ `numCompare ${RIGHT} $1` -eq 1 ]; then
    if [ -z ${START} ]
    then
      START=${COORDINATE}
    fi
    FORK_DIRECTION="RIGHT"
    PREVIOUS_COORDINATE=${COORDINATE}
    continue
  else
    if [ ! -z ${START} ]
    then
      END=${PREVIOUS_COORDINATE}
      let "DIFF = END - START"
      echo -e "${READ_COORDINATES}\t${START}\t${END}\t${DIFF}\t${FORK_DIRECTION}"
      START=""
      FORK_DIRECTION=""
    fi
  fi
done

# Dealing with the fucking last line.
if [ ! -z ${START} ]
then
  END=${PREVIOUS_COORDINATE}
  let "DIFF = END - START"
  echo -e "${READ_COORDINATES}\t${START}\t${END}\t${DIFF}\t${FORK_DIRECTION}"
fi

exit 0
