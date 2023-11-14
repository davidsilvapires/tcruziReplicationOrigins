#!/usr/bin/env bash

###################################################################################################
#                                         bedtoolsMain                                            #
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
# Usage: saveCommand script/bedtoolsMain.bash                                                     #
#                                                                                                 #
# Copyleft (ɔ) 2023 by Thiago Andrade Franco <thiago.franco.esib@esib.butantan.gov.br>            #
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

# Sorting gff files by chromosomes and by starting position.
for i in \`ls input/*.gff | xargs -i basename {}\`
do
 sort -k1,1 -k2,2n input/\${i} > output/\${i%.gff}-sorted.gff
done

# Intersect the bed files from DNA replication Origins with core compartment
for FILE in \`ls output/*-sorted.bed | xargs -i basename {}\`
do
bedtools intersect -a output/\${FILE} \
  -b output/core-sorted.gff \
  > output/\${FILE%-sorted.bed}OverlappingCore.tsv \
  2> log/\${FILE%-sorted.bed}OverlappingCore.err \


done
 
# Intersect the bed files from DNA replication Origins with disruptive compartment
for FILE in \`ls output/*-sorted.bed | xargs -i basename {}\`
do
bedtools intersect -wo -a output/\${FILE} \
  -b output/disruptive-sorted.gff \
  > output/\${FILE%-sorted.bed}OverlappingDisruptive.tsv \
  2> log/\${FILE%-sorted.bed}OverlappingDisruptive.err \


done   

# Intersect the bed files from DNA replication Origins with both compartment
for FILE in \`ls output/*-sorted.bed | xargs -i basename {}\`
do
bedtools intersect -wo -a output/\${FILE} \
  -b output/both-sorted.gff \
  > output/\${FILE%-sorted.bed}OverlappingBoth.tsv \
  2> log/\${FILE%-sorted.bed}OverlappingBoth.err \


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

