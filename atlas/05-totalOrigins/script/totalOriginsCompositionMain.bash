#!/usr/bin/env bash

###################################################################################################
#                                      totalOriginsCompositionMain.bash                           #  
#                                                                                                 #
#  This script allows the creation of Total origins database. In this manner, We utilized the     #
# linux command cat to concatenate the dataset of predominants, flexible and dormant origins to   #
# create the total origins. It writes the Slurm job script and submits all of them.               #
#                                                                                                 #
# Usage: saveCommand script/totalOriginsCompositionMain.bash                                      #
#                                                                                                 #
# Copyleft (ɔ) 2023 by  Thiago Andrade Franco                                                     #
#                       thiago.franco.esib@esib.butantan.gov.br                                   #
#    11/10/2023: First version.                                                                   #
###################################################################################################

# Writing the job script to submits at Slurm
cat> job/totalOriginsCompositionMain.slurm <<EOI
#!/usr/bin/env bash

# Classify the origins of each bank
for FILE in \`ls input/*.bed | xargs -i basename {} \`
do
  awk 'OFS="\t" {print \$1, \$2, \$3 }' input/\${FILE} \\
  > output/\${FILE%.bed}.csv
done


# Sorting bed files by chromosomes and by starting position.
for FILE in \`ls output/*.csv | xargs -i basename {} \`
do
  sort -k1,1 -k2,2n -k3,3n output/\${FILE} \\
  > output/\${FILE%.csv}-sorted.tsv
done


# Classify the origins of each bank
for FILE in \`ls output/*.tsv | xargs -i basename {} \`
do
  awk 'OFS="\t" {print \$1, \$2, \$3 , output/"\${FILE%-sorted.tsv}"}' output/\${FILE} \\
  > output/\${FILE%-sorted.tsv}.bed
done

# To create the total origns  dataset, we will concatenate the predominant, flexible, and dormant
# origins. 
  cat output/*.bed > output/totalOrigins-notSorted.bed

# Sorting the bed file by chromosomes and by starting position of the total origns file
  sort -k1,1 -k2,2n -k3,3n output/totalOrigins-notSorted.bed \\
  > output/totalOrigins.bed

exit 0 

EOI

# Submitting the job to Slurm.
  /usr/bin/nice -n 19 /usr/bin/time  --verbose \
    --output=log/sbatch-totalOriginsCompositionMain.time memusg \
    --output-to-file log/sbatch-totalOriginsCompositionMain.memusg \
    --time \
    --shell "saveCommand sbatch --nodes 1 \
    --ntasks ${NUM_OF_CPUS} \
    --mem ${MEMORY_SIZE} \
    -o log/slurm-%A.out \
    -J totalOriginsComposition job/totalOriginsCompositionMain.slurm \
    2> log/sbatch-totalOriginsCompositionMain.err |
    tee log/sbatch-totalOriginsCompositionMain.out"

# Check how are going all the running awk processes. 
 echo "To check the runnning processes, execute one of the following commands:"
 echo "   $> watch -n 1 sequeue"

exit 0

