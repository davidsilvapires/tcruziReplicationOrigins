#!/usr/bin/env bash

###################################################################################################
#                                         awkMain                                            #
#                                                                                                 #
#                                                                                                 #
# This tutorial allows you to record the origins of DNA replication along the genome's compartment# 
#  We generated  and analyzed the compartment of the genome where the origins of DNA replication  #
# are located using  the "intersect" tool of the "bedtools" program. In summary, bedtools will    #
# compare the genomic coordinates of the DNA replication origins to the genomic coordinates of the#
# reference genome  and output a file containing the intersection regions. As a result, we will   #
# examine each DNA replication origins file to determine where genomic areas intersected.         #
# It writes the Slurm job script and submits all of them.                                         #
#                                                                                                 #
# Usage: saveCommand script/awkMain.bash                                                     #
#                                                                                                 #
# Copyleft (ɔ) 2023 by Thiago Andrade Franco <thiago.franco.esib@esib.butantan.gov.br>            #
#    12/07/2023: First version.                                                                   #
###################################################################################################

# Writing the job script to submits at Slurm
cat> job/awkMain.slurm <<EOI
#!/usr/bin/env bash

FILENAME=`basename \${BEDFILE}`
GFF_FILENAME=`basename \${GFFFILE}`

echo -e "Chrom\tStart (bp)\tEnd (bp)\tFeature\t Feature Start (bp)\tFeature End(bp)\tStrand\tID\tDescription\t"

IFS='
'
# Pulga atrás da orelha: será que a linha abaixo realmente funciona? Testar em exemplos maiores.
for LINE in `sort -k1,1 -k2,2n -k3,3n ${FILENAME}`
do
  #REP=`echo -e ${LINE} | cut -f11`
  #if [[ "${REP}" == "-" ]]
  #then
  #  continue
  #else
    CHROM=`echo -e ${LINE} | cut -f1`
    REP_START=`echo -e ${LINE} | cut -f2`
    REP_END=`echo -e ${LINE} | cut -f3`
    NEXT_CHROM="NO"
    # echo "${CHROM} ${FORK_START} ${FORK_END}"
    for GFF_LINE in `tail -n +2 ${GFF_FILENAME} | sort -k1,1 -k3,3n -k4,4n`
    do
      # echo "GFF_LINE = ${GFF_LINE}"
      GFF_CHROM=`echo -e ${GFF_LINE} | cut -f1`
      if [ "${NEXT_CHROM}" == "YES" ] && [ "${GFF_CHROM}" != "${CHROM}" ]
      then
        # echo "Ueba! Break"
        break
      fi
      if [ "${GFF_CHROM}" == "${CHROM}" ]
      then
        # echo "Chrom atingido!:${CHROM}"
        NEXT_CHROM="YES"
        FEATURE_START=`echo -e ${GFF_LINE} | cut -f3`
        # echo "${FEATURE_START}"
        if [ "${FEATURE_START}" -gt "${REP_END}" ]
        then
          continue
        fi
        FEATURE_END=`echo -e ${GFF_LINE} | cut -f4`
        if [ "${FEATURE_END}" -lt "${REP_START}" ]
        then
          continue
        fi
        GENE_STRAND=`echo -e ${GFF_LINE} | cut -f5`
        TABLE_PREFIX=`echo -e ${LINE} | cut -f1-6`
        TABLE_SUFFIX=`echo -e ${LINE} | cut -f11-14`
        if [ "${GENE_STRAND}" == "+" ]
        then
          echo -e "${TABLE_PREFIX}\t${TABLE_SUFFIX}\t$(echo -e ${GFF_LINE} | cut -f2-7)"
#          echo -e -n "${TABLE_PREFIX}\t${TABLE_SUFFIX}\t$(echo -e ${GFF_LINE} | cut -f2-7)\t+\t"
        fi
        if [ "${GENE_STRAND}" == "-" ]
        then
          echo -e "${TABLE_PREFIX}\t${TABLE_SUFFIX}\t$(echo -e ${GFF_LINE} | cut -f2-7)"
#          echo -e -n "${TABLE_PREFIX}\t${TABLE_SUFFIX}\t$(echo -e ${GFF_LINE} | cut -f2-7)\t-\t"
        fi
        if [ "${GENE_STRAND}" == "." ]
        then
          echo -e "${TABLE_PREFIX}\t${TABLE_SUFFIX}\t$(echo -e ${GFF_LINE} | cut -f2-7)"
#          echo -e -n "${TABLE_PREFIX}\t${TABLE_SUFFIX}\t$(echo -e ${GFF_LINE} | cut -f2-7)\t.\t"
        fi
#        FORK_DIRECTION=`echo -e ${LINE} | cut -f10`
#        # echo FORK_DIRECTION=${FORK_DIRECTION}
#        READ_STRAND=`echo -e ${LINE} | cut -f6`
#        # echo READ_STRAND=${READ_STRAND}
#        # GENE_STRAND = TRANSCRIPTION DIRECTION
#        if [ "${READ_STRAND}" == "fwd" ] && [ "${FORK_DIRECTION}" == "LEFT" ] && [ "${GENE_STRAND}" == "+" ]
#        then
#          echo "YES"
#          continue
#        fi
#        if [ "${READ_STRAND}" == "rev" ] && [ "${FORK_DIRECTION}" == "RIGHT" ] && [ "${GENE_STRAND}" == "+" ]
#        then
#          echo "YES"
#          continue
#        fi
#        if [ "${READ_STRAND}" == "fwd" ] && [ "${FORK_DIRECTION}" == "RIGHT" ] && [ "${GENE_STRAND}" == "-" ]
#        then
#          echo "YES"
#          continue
#        fi
#        if [ "${READ_STRAND}" == "rev" ] && [ "${FORK_DIRECTION}" == "LEFT" ] && [ "${GENE_STRAND}" == "-" ]
#        then
#          echo "YES"
#          continue
#        fi
#        echo "NO"
      fi
    done
  #fi
done

exit 0

EOI

# Submitting the job to Slurm.
/usr/bin/nice -n 19 /usr/bin/time --verbose --output=log/sbatch-awkMain.time memusg \
  --output-to-file log/sbatch-bedtoolsMain.memusg --time --shell "saveCommand sbatch --nodes 1 \
  --ntasks ${NUM_OF_CPUS} --mem ${MEMORY_SIZE} -o log/slurm-%A.out -J bedtoolsMain job/bedtoolsMain.slurm \
  2> log/sbatch-bedtoolsMain.err | tee log/sbatch-bedtoolsMain.out"

# Check how are going all the running awk processes. 
 echo "To check the runnning processes, execute one of the following commands:"
 echo "   $> watch -n 1 sequeue"

exit 0

