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
# Copyleft (ɔ) 2023 by Marcela de Oliveira Vitarelli and Thiago Andrade Franco                    #
#                    <vitarelli.marcela@gmail.com> and thiago.franco.esib@esib.butantan.gov.br    #
#    12/07/2023: First version.                                                                   #
###################################################################################################

# Writing the job script to submits at Slurm
cat> job/bedtoolsMain.slurm <<EOI
#!/usr/bin/env bash

# Sorting bed files by chromosomes and by starting position.
for i in \`ls input/*.bed | xargs -i basename {}\`
do
 sort -k1,1 -k2,2n input/\${i} > output/\${i%.bed}-sorted.bed
done


# Sorting gff file from the genome
GFF_FILENAME=`basename \${GENOME}`
  sort -k1,1 -k4,4n -k5,5n input/\${GFF_FILENAME} > output/\${GFF_FILENAME%.gff}-sorted.gff

#Intersect the bed files from DNA replication Origins with genome 
for FILE in \`ls output/*-sorted.bed | xargs -i basename {}\`
do
bedtools intersect -wo -a output/\${FILE} \
  -b output/\${GFF_FILENAME%.gff}-sorted.gff \
  > output/\${FILE%-sorted.bed}Overlapping\${GFF_FILENAME%.gff}.tsv \
  2> log/\${FILE%-sorted.bed}Overlapping\${GFF_FILENAME%.gff}.err \

done
    
   exit 0 

EOI

# Submitting the job to Slurm.
/usr/bin/nice -n 19 /usr/bin/time --verbose --output=log/sbatch-bedtoolsMain.time memusg \
  --output-to-file log/sbatch-bedtoolsMain.memusg --time --shell "saveCommand sbatch --nodes 1 \
  --ntasks ${NUM_OF_CPUS} --mem ${MEMORY_SIZE} -o log/slurm-%A.out -J bedtoolsMain job/bedtoolsMain.slurm \
  2> log/sbatch-bedtoolsMain.err | tee log/sbatch-bedtoolsMain.out"

# Check how are going all the running awk processes. 
 echo "To check the runnning processes, execute one of the following commands:"
 echo "   $> watch -n 1 sequeue"

exit 0

