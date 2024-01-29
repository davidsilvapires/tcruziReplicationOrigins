#!/usr/bin/env bash

###################################################################################################
#                                         bedtoolsMain                                            #  
#                                                                                                 #
#  This script allows the creation pf Predominants origins database. In this manner, we utilized  #
# the bedtools intersect to generate the predominant origins database. In summary, we found which #
# MFA-Seq DNA replication origins had genomic coordinates overlapping with Orc1cdc6 chip-seq peaks#
#  As a result, the MFA-seq origins that overlap with Orc1Cdc6 will be considered the predominant #
# origins. It writes the Slurm job script and submits all of them.                                #
#                                                                                                 #
# Usage: saveCommand script/bedtoolsMain.bash                                                     #
#                                                                                                 #
# Copyleft (ɔ) 2023 by Thiago Andrade Franco                                                      #
#                    < thiago.franco.esib@esib.butantan.gov.br                                    #
#    12/14/2023: First version.                                                                   #
###################################################################################################

# Writing the job script to submits at Slurm
cat> job/bedtoolsMain.slurm <<EOI
#!/usr/bin/env bash

# Building the bed files.
for FILE in \`ls input/*.bed | xargs -i basename {}\`
do 
  cut -f1-3 input/\${FILE} | 
  sort -u  \\
  > output/\${FILE}
done

# Sorting bed files by chromosomes and by starting position.
for i in \`ls output/*.bed | xargs -i basename {}\`
do
 sort -k1,1 -k2,2n -k3,3n input/\${i} | 
 sort -u \\
 > output/\${i%.bed}-sorted.bed
done

# Get fasta sequence of bed file to process the meme analisys. 
for FILE in \`ls output/*-sorted.bed | xargs -i basename {}\`
do
  bedtools getfasta -fi \${GENOME_FASTA} \\
    -bed output/\${FILE} \\
    > output/\${FILE%-sorted.bed}.fasta
done

exit 0 

EOI

# Submitting the job to Slurm.
  /usr/bin/nice -n 19 /usr/bin/time \
  --verbose \
  --output=log/sbatch-bedtoolsMain.time memusg \
  --output-to-file log/sbatch-bedtoolsMain.memusg \
  --time --shell "saveCommand sbatch --nodes 1 \
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

