#!/usr/bin/env bash

###################################################################################################
#                                         bedtoolsMain                                            #
#                                                                                                 #
#  This script calculates the distance between the replication origins' center point and the      #
# center point of orc1cdc6. We utilize the closest tool from Bedtools for this. This tool will    #
# calculate the distance between the DNA replication origins and the Chip-Seq Orc1cdc6 peaks. It  #
# writes the Slurm job script and submits all of them.                                            #
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

# defining the center point of genomic coordinates of the DNA replication origins
for FILE in \`ls output/*-sorted.bed | xargs -i basename {}\`
do
   awk 'function abs(v) { return v < 0 ? -v : v } function min(a, b) { return a > b ? b : a } \
     OFS = "\t" {print \$1, int(abs(\$3 - \$2) / 2) + min(\$3, \$2), int(abs(\$3 - \$2) / 2) + \
     min(\$3, \$2) }' output/\${FILE}  > output/\${FILE%-sorted.bed}-CenterPoint.bed
done

for FILE in \`ls output/*CenterPoint.bed | xargs -i basename {}\`
do
bedtools closest -d -t first -a output/\${FILE} \
-b output/orc1cdc6-CenterPoint.bed \
> output/distanceBetweenOrc1cdc6-\${FILE%-CenterPoint.bed}.tsv \
2> log/distanceBetweenOrc1cdc6-\${FILE%-CenterPoint.bed}.err 

done

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

