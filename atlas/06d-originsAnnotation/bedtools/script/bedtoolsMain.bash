#!/usr/bin/env bash

###################################################################################################
#                                         bedtoolsMain                                            #
#                                                                                                 #
#                                                                                                 #
#    This script allows you to record the origins of DNA replication in the genome. We generated  #
# and analyzed the sections of the genome where the origins of DNA replication are located using  #
# the "intersect" tool of the "bedtools" program. In summary, bedtools will compare the genomic   #
# coordinates of the DNA replication origins to the genomic coordinates of the reference genome   #
# and output a file containing the intersection regions. As a result, we will examine each DNA    #
# replication origins file to determine where genomic areas intersected. It writes the Slurm job  #
# script and submits all of them.                                                                 #
#                                                                                                 #
# Usage: saveCommand script/bedtoolsMain.bash                                                     #
#                                                                                                 #
# Copyleft (ɔ) 2023  Thiago Andrade Franco                                                        #
#                    thiago.franco.esib@esib.butantan.gov.br                                      #
#    08/07/2023: First version.                                                                   #
###################################################################################################

# Writing the job script to submits at Slurm
cat> job/bedtoolsMain.slurm <<EOI
#!/usr/bin/env bash

GFF_FILENAME=`basename \${GENOME}`

# preparing the input gff for analysis of genome features. In this step we will use the phyton
# program script developed by Dr. Alex Raniery Lima.
  python script/annotatePolycistron_SSR.py\\
  -g \${GENOME} \\
  -o output/\${GFF_FILENAME} 

# Sorting gff file from the genome
  sort -k1,1 -k4,4n -k5,5n output/\${GFF_FILENAME} \\
  > output/\${GFF_FILENAME%.gff}-sorted.gff

# Correcting the genomic coordinates of gff file  with polycistron
## greater than the end
< output/\${GFF_FILENAME%.gff}-sorted.gff \
  awk '\$4 > \$5 {printf("%s\t%s\t%s\t%d\t%d\t%s\t%s\t%s\t%s\n", \$1, \$2,\
  \$3, \$5, \$4, \$6, \$7, \$8, \$9)}' \\
  > output/\${GFF_FILENAME%.gff}-startGreaterThanEnd.gff   
## removing the header 
  grep -v "^#" output/\${GFF_FILENAME%.gff}-startGreaterThanEnd.gff \\
   > output/\${GFF_FILENAME%.gff}-startGreaterThanEnd-withoutHeader.gff

## less than the end
  < output/\${GFF_FILENAME%.gff}-sorted.gff \
  awk '\$4 < \$5 {print \$0}'  \
  > output/\${GFF_FILENAME%.gff}-startLessThanEnd.gff

## sort the files by chromosomes, start and end
  sort -k1,1 -k4,4n -k5,5n output/\${GFF_FILNAME%.gff}-startGreaterThanEnd-withoutHeader.gff \
  > output/\${GFF_FILNAME%.gff}-startGreaterThanEnd-withoutHeader-sorted.gff
  sort -k1,1 -k4,4n -k5,5n output/\${GFF_FILENAME%.gff}-startLessThanEnd.gff \
  > output/\${GFF_FILENAME%.gff}-startLessThanEnd-sorted.gff

## concatened
  cat output/\${GFF_FILNAME%.gff}-startGreaterThanEnd-withoutHeader-sorted.gff  \
      output/\${GFF_FILENAME%.gff}-startLessThanEnd-sorted.gff \
      > output/\${GFF_FILENAME%.gff}-polycistron-SSR.gff

      # Sorting bed files by chromosomes and by starting position.
for i in \`ls input/*.bed | xargs -i basename {}\`
do
  sort -k1,1 -k2,2n input/\${i} \\
  > output/\${i%.bed}-sorted.bed
done

# Intersect the bed files from DNA replication Origins with genome
for FILE in \`ls output/*-sorted.bed | xargs -i basename {}\`
do
bedtools intersect -wo -a output/\${FILE} \\
  -b output/\${GFF_FILENAME%.gff}-polycistron-SSR.gff \\
  > output/\${FILE%-sorted.bed}Overlapping\${GFF_FILENAME%.gff}.tsv \\
  2> log/\${FILE%-sorted.bed}Overlapping\${GFF_FILENAME%.gff}.err 

done

exit 0

EOI

# Submitting the job to Slurm.
  /usr/bin/nice -n 19 /usr/bin/time \
  --verbose \
  --output=log/sbatch-bedtoolsMain.time memusg \
  --output-to-file log/sbatch-bedtoolsMain.memusg \
    --time --shell "saveCommand sbatch --parsable \
    --nodes 1 \
  --ntasks ${NUM_OF_CPUS} \
  --mem ${MEMORY_SIZE} \
  -o log/slurm-%A.out \
  -J bedtoolsMain job/bedtoolsMain.slurm \
  2> log/sbatch-bedtoolsMain.err |
    tee log/sbatch-bedtoolsMain.out"

# Check how are going all the running awk processes.
echo "To check the runnning processes, execute one of the following commands:"
echo "   $> watch -n 1 sequeue"

exit 0
