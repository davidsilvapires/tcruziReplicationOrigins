#!/usr/bin/env bash

###################################################################################################
#                                        deeptoolsMain                                            #  
#                                                                                                 #
#  This script allows verifying the co-location of the orc1cdc6 peaks in the genes of the         #
# multigenic families, in the features of the genome and also in the origins of replication       #
# identified in this  work. For this analysis we will use the deeptools compute matrix and        #
# plotheatmap tools. Script  also writes another script to run on Slurm.                          #
#                                                                                                 #
# Usage: saveCommand script/bedtoolsMain.bash                                                     #
#                                                                                                 #
# Copyleft (ɔ) 2023 by Thiago Andrade Franco                                                      #
#                   < thiago.franco.esib@esib.butantan.gov.br>                                    #
#    12/07/2023: First version.                                                                   #
###################################################################################################

# Writing the job script to submits at Slurm
cat> job/deeptoolsMain5.slurm <<EOI
#!/usr/bin/env bash

# Creating a dataset of random sequences from T. cruzi genome with 6kb
  seqkit sliding -s 6000 -W 6000 \\
  input/tryCru-clb1.fasta \\
  -o output/slicedGenome6kb.fasta

# Count the number of each nitrogenous (N) base
  faCount output/slicedGenome6kb.fasta \\
  > output/countbase.csv

#  Select only sequences that N < 500 
  awk -F "\t" '\$7 < 500 {print \$0}' \\
  output/countbase.csv |
  sort -k7,7n \\
  > output/sequencesLessThan500N.csv 

# Creating radom sequence to control
IFS='
'
for LINE in \`cat output/sequencesLessThan500N.csv | cut -f1\`
do
  CHROM=\`echo \${LINE} | cut -d "_" -f1\`
  CHROM_START=\`echo \${LINE} | cut -d ":" -f2 | cut -d "-" -f1\`
  CHROM_END=\`echo \${LINE} | cut -d ":" -f2 | cut -d "-" -f2\`
  
  echo -e "\${CHROM}\t\${CHROM_START}\t\${CHROM_END}"
done > output/control.bed

# Generating bed files with 6kb random sequences for Histone heatmap control.
## To predominant
  shuf -n83 output/control.bed | 
  sort -k1,1 -k2,2n -k3,3n \\
  > output/radon_sequencesPredominants.bed 

## To flexible
 shuf -n662 output/control.bed | 
  sort -k1,1 -k2,2n -k3,3n \\
  > output/radon_sequencesFlexible.bed 

## To Dormant
 shuf -n451 output/control.bed | 
  sort -k1,1 -k2,2n -k3,3n \\
  > output/radon_sequencesDormant.bed 

exit 0
EOI

# Submitting the job to Slurm.
  /usr/bin/nice -n 19 /usr/bin/time \
  --verbose \
  --output=log/sbatch-deeptoolsMain.time memusg \
  --output-to-file log/sbatch-deeptoolsMain.memusg \
  --time --shell "saveCommand sbatch --parsable \
  --nodes 1 \
  --ntasks ${NUM_OF_CPUS} \
  --mem ${MEMORY_SIZE} \
  -o log/slurm-%A.out \
  -J deeptoolsMain5 job/deeptoolsMain5.slurm \
  2> log/sbatch-deeptoolsMain5.err | 
  tee log/sbatch-deeptoolsMain5.out"

# Check how are going all the running awk processes. 
echo "To check the runnning processes, execute one of the following commands:"
echo "   $> watch -n 1 sequeue"

exit 0



