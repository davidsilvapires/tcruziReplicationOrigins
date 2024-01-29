#!/usr/bin/env bash

###################################################################################################
#                                         bedtoolsMain                                            #  
#                                                                                                 #
#   In this guide, we will analyze the region of DNA replication origins. We will check which     #
# features are around us and their distribution. To do this, we will use an improvised script and #
# the bedtools intersect tools.                                                                   #
#                                                                                                 #
# Usage: saveCommand script/bedtoolsMain.bash                                                     #
#                                                                                                 #
# Copyleft (ɔ) 2023 by  Thiago Andrade Franco and David da Silva Pires                            #
#                < david.pires@butantan.gov.br > and < thiago.franco.esib@esib.butantan.gov.br>   #
#    12/12/2023: First version.                                                                   #
###################################################################################################

# Writing the job script to submits at Slurm
cat > job/bedtoolsMain.slurm << EOI
#!/usr/bin/env bash

# Sorting bed files by chromosomes and by starting position.
for i in \`ls input/*.bed | xargs -i basename {} \`
do
  sort -k1,1 -k2,2n -k3,3n input/\${i} \\
  > output/\${i%.bed}-sorted.csv
done

# Creating bed3 file to sorted origins
for FILE in \`ls output/*-sorted.csv | xargs -i basename {}\`
do
  awk 'OFS="\t" {print \$1, \$2, \$3}' output/\${FILE} \\
  > output/\${FILE%-sorted.csv}.bed
done 

# Creating the gff file containing the features only with disjoint sets
  < input/genome.gff grep -v '^#' | 
  grep -v -P '.*\t.*\tmRNA\t.*\t.*\t.*\t.*\t.*\t.*' | 
  grep -v -P '.*\t.*\tprotein_coding_gene\t.*\t.*\t.*\t.*\t.*\t.*' | 
  grep -v -P '.*\t.*\texon\t.*\t.*\t.*\t.*\t.*\t.*' | 
  grep -v -P '.*\t.*\tncRNA_gene\t.*\t.*\t.*\t.*\t.*\t.*' |
  grep -v -P '.*\t.*\tpseudogenic_transcript\t.*\t.*\t.*\t.*\t.*\t.*' \\
  > output/genomeDisjointFeatures.gff

# Sorting the genomeDisjoinFeatures.
  sort  -k1,1 -k4,4n -k5,5n output/genomeDisjointFeatures.gff \\
  > output/genomeDisjointFeatures-sorted.gff

# Creating the bed file to up and downstream windows around origins and intersect with genome
# disjoint features.
for ORIGIN_TYPE in \${ORIGIN_TYPES}
do
   ORIGIN_ID=0
  # Creating a BED file with up- and downstream 3k windows around origins.
  IFS='
'
   for LINE in \`cat output/\${ORIGIN_TYPE}.bed\`
   do
      CHROM=\`echo \${LINE} | cut -f1\`
      CHROM_SIZE=\`grep -P "\${CHROM}\t" \${GENOME_CHROMSIZE} | cut -f2\`
      UPSTREAM_WINDOW_END=\`echo  \${LINE} | cut -f2\`
      DOWNSTREAM_WINDOW_START=\`echo \${LINE} | cut -f3\`
      let "UPSTREAM_WINDOW_START = UPSTREAM_WINDOW_END - 3000"
      if [ "\${UPSTREAM_WINDOW_START}" -lt "0" ]
      then
         UPSTREAM_WINDOW_START=0
      fi
      let "DOWNSTREAM_WINDOW_END = DOWNSTREAM_WINDOW_START + 3000"
      if [ "\${DOWNSTREAM_WINDOW_END}" -gt "\${CHROM_SIZE}" ]
      then 
            DOWNSTREAM_WINDOW_END=\${CHROM_SIZE}
      fi
      echo -e "\${CHROM}\t\${UPSTREAM_WINDOW_START}\t\${UPSTREAM_WINDOW_END}\t\${ORIGIN_ID}-u"
      echo -e "\${CHROM}\t\${DOWNSTREAM_WINDOW_START}\t\${DOWNSTREAM_WINDOW_END}\t\${ORIGIN_ID}-d"
      let "ORIGIN_ID = ORIGIN_ID + 1"
   done >> output/\${ORIGIN_TYPE}-3kWindow.bed4
   bedtools intersect -wo -a output/\${ORIGIN_TYPE}-3kWindow.bed4 \\
      -b output/genomeDisjointFeatures-sorted.gff \\
      > output/\${ORIGIN_TYPE}-intersect-genomeDisjointFeatures.tsv
done

# Dividing the intersection files by feature.
IFS=' '
for ORIGIN_TYPE in \${ORIGIN_TYPES}
do
  # Building a file that contains only one column with the features that compose the 7th column.
  cut -f7 output/\${ORIGIN_TYPE}-intersect-genomeDisjointFeatures.tsv | sort -u > \\
    output/\${ORIGIN_TYPE}-intersect-genomeDisjointFeatures-onlyFeatures.txt

  IFS='
'
  for LINE in \`cat output/\${ORIGIN_TYPE}-intersect-genomeDisjointFeatures-onlyFeatures.txt\`
  do
    FEATURE=\`echo \${LINE} | cut -f7\`
    grep -P ".*\t.*\t.*\t.*\t.*\t.*\t\${FEATURE}\t.*\t.*\t.*\t.*\t.*\t.*\t.*"  \\
    output/\${ORIGIN_TYPE}-intersect-genomeDisjointFeatures.tsv \\
    > output/\${ORIGIN_TYPE}-intersect-genomeDisjointFeatures-\${FEATURE}.tsv
  done
done

# Counting the number of features that participates in an intersection.
for FILE in output/*-intersect*-*.tsv
do
  FEATURE_NAME=\`echo "\${FILE}" | sed 's|output/.*-intersect-.*-\(.*\).tsv|\1|'\`
  # Counting the number of features that participates in an intersection.
  NUM_OF_FEATURES=\`cut -f13 \${FILE} | cut -d';' -f1 | cut -d= -f2 | sort | wc -l\`
  echo "Number of neighbouring \${FEATURE_NAME}: \${NUM_OF_FEATURES}" >> \\
    output/numberOfFeatures.txt

 for LINE in \`cat \${FILE}\`
    do 
      CHROM=\`echo \${LINE} | cut -f1\`
      WINDOW_START=\`echo \${LINE} | cut -f2\`
      WINDOW_END=\`echo \${LINE} | cut -f3\`
      FEATURE_START=\`echo \${LINE} | cut -f8\` 
      FEATURE_END=\`echo \${LINE} | cut -f9\`

      let "START_DIFFERENCE = FEATURE_START - WINDOW_START"
      let "END_DIFFERENCE = FEATURE_END - WINDOW_END"

      if [ "\${START_DIFFERENCE}" -ge "0" ] 
      then 
      INTERSECTION_START=\${FEATURE_START}
      else
      INTERSECTION_START=\${WINDOW_START}
      fi
      if [ "\${END_DIFFERENCE}" -ge "0" ]
      then
      INTERSECTION_END=\${WINDOW_END}
      else
      INTERSECTION_END=\${FEATURE_END}
      fi
      echo -e "\${CHROM}\t\${INTERSECTION_START}\t\${INTERSECTION_END}"
    done > \${FILE%.tsv}-justIntersection.bed3
done

# Counting the number of bases of features that participates in an intersection.
for FILE in output/*-justIntersection.bed3
do
  FEATURE_NAME=\`echo "\${FILE}" | sed 's|output/.*-intersect-.*-\(.*\).tsv|\1|'\`
  sort -k1,1 -k2,2n -k3,3n \${FILE} > \${FILE%.bed3}-sorted.bed  
  bedtools merge -i \${FILE%.bed3}-sorted.bed > \${FILE%-sorted.bed}-merged.bed3
  INTERSECTION_BASES=0
  for LINE in \`cat \${FILE%-sorted.bed}-merged.bed3\`
  do
    CHROM_START=\`echo \${LINE} | cut -f2\`
    CHROM_END=\`echo \${LINE} | cut -f3\`

    let "SIZE = CRHOM_END - CHROM_START"
    let "INTERSECTION_BASES = INTERSECTION_BASES + SIZE"
  done
  echo "Number of bases of neighbouring \${FEATURE_NAME}: \${INTERSECTION_BASES}" >> \\
    output/numberOfFeatureBases.txt
done

exit 0

EOI

# Submitting the job to Slurm.
  /usr/bin/nice -n 19 /usr/bin/time \
  --verbose \
  --output=log/sbatch-bedtoolsMain.time memusg \
  --output-to-file log/sbatch-bedtoolsMain.memusg  \
  --time --shell "saveCommand sbatch --nodes 1 \
  --ntasks ${NUM_OF_CPUS}  \
  --mem ${MEMORY_SIZE}  \
  -o log/slurm-%A.out   \
  -J bedtoolsMain job/bedtoolsMain.slurm \
  2> log/sbatch-bedtoolsMain.err | 
  tee log/sbatch-bedtoolsMain.out"

# Check how are going all the running awk processes. 
 echo "To check the runnning processes, execute one of the following commands:"
 echo "   $> watch -n 1 sequeue"

exit 0

