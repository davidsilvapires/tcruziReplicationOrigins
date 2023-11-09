#!/usr/bin/env bash

FILENAME=${1}
GFF_FILENAME=${2}

HEADER=`head -1 ${1}`
echo -e "${HEADER}\tFeature Start (bp)\tFeature End(bp)\tStrand\tFeature\tID\tDescription\tTranscription Direction (+/-)\tConflicts"

IFS='
'
# Pulga atrás da orelha: será que a linha abaixo realmente funciona? Testar em exemplos maiores.
for LINE in `tail -n +2 ${FILENAME} | sort -k2,2 -k3,3n -k4,4n -k6,6n -k7,7n -k8,8n -k10,10`
do
  CHROM=`echo ${LINE} | cut -f2`
  FORK_START=`echo ${LINE} | cut -f7`
  FORK_END=`echo ${LINE} | cut -f8`
  NEXT_CHROM="NO"
  # echo "${CHROM} ${FORK_START} ${FORK_END}"
  for GFF_LINE in `tail -n +2 ${GFF_FILENAME}`
  do
    # echo "GFF_LINE = ${GFF_LINE}"
    GFF_CHROM=`echo ${GFF_LINE} | cut -d, -f1`
    if [ "${NEXT_CHROM}" == "YES" ] && [ "${GFF_CHROM}" != "${CHROM}" ]
    then
      # echo "Ueba! Break"
      break
    fi
    if [ "${GFF_CHROM}" == "${CHROM}" ]
    then
      # echo "Chrom atingido!:${CHROM}"
      NEXT_CHROM="YES"
      FEATURE_START=`echo ${GFF_LINE} | cut -d, -f2`
      # echo "${FEATURE_START}"
      if [ "${FEATURE_START}" -gt "${FORK_END}" ]
      then
        continue
      fi
      FEATURE_END=`echo ${GFF_LINE} | cut -d, -f3` 
      if [ "${FEATURE_END}" -lt "${FORK_START}" ]
      then
        continue
      fi
      STRAND=`echo ${GFF_LINE} | cut -d, -f4`
      if [ "${STRAND}" == "+" ] 
      then
        echo -e -n "${LINE}\t$(echo ${GFF_LINE} | cut -d, -f2-7 | tr ',' '\t')\t+\t"
      elif [ "${STRAND}" == "-" ]
      then
        echo -e -n "${LINE}\t$(echo ${GFF_LINE} | cut -d, -f2-7 | tr ',' '\t')\t-\t"
      fi
      FORK_DIRECTION=`echo -e ${LINE} | cut -f10`
      # echo FORK_DIRECTION=${FORK_DIRECTION}
      if [ "${STRAND}" == "+" ] && [ "${FORK_DIRECTION}" == "LEFT" ]
      then
        echo "YES"
        continue
      fi
      if [ "${STRAND}" == "-" ] && [ "${FORK_DIRECTION}" == "RIGHT" ]
      then
        echo "YES"
        continue
      fi
      echo "NO"
    fi
  done
done

exit 0
