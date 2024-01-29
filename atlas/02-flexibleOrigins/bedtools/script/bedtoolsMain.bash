#!/usr/bin/env bash

###################################################################################################
#                                         bedtoolsMain                                            #  
#                                                                                                 #
#  This script allows the creation of Flexible origins database. In summary, the Flexible origins #
#  database contains DNA replication origins discovered by D-Nascent that have intersection of    #
#  their genomic coordinates with the coordinates of the Orc1Cdc6 peaks but not with the          #
# coordinates of the predominant origins.  In this manner, we utilized  the bedtools intersect and#
# bedtools subtract to generate the flexible origins database. It writes the Slurm job script     #
# and submits all of them.                                                                        #
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
for i in \`ls input/\`
do
 sort -k1,1 -k2,2n input/\${i} > output/\${i%.bed}-sorted.bed
done

# Subtract predominant origins from dnascent. If any overlap is found the entire feature is removed.
bedtools subtract -A -a output/dnascent-sorted.bed \
  -b output/predominantOrigins-sorted.bed \
  > output/subtractDnascent-PredominantOrigins.bed \
  2>log/subtractDnascent-PredominantOrigins.err 

# Intersect function to identify overlapping peaks between dnascent-predominantes origins with a
# +/-3kb window and orc1cdc6 => Flexible origins.
bedtools intersect -wo -a output/subtractDnascent-PredominantOrigins.bed \
  -b output/orc1cdc6-sorted.bed \
  > output/intersectSubtractedDnascent-PredominantsWithOrc1cdc6.bed \
  2>log/intersectSubtractedDnascent-PredominantsWithOrc1cdc6.err

# Generating flexible origins file.
cat output/intersectSubtractedDnascent-PredominantsWithOrc1cdc6.bed| cut -f1-3 | sort -u \
  > output/flexibleOrigins.bed

# Generating flexible origins file. -NotUniq
cat output/intersectSubtractedDnascent-PredominantsWithOrc1cdc6.bed| cut -f1-3 \
  > output/flexibleOrigins-NotUniq.bed

# Generating orc1cdc6 coordinates that overlap with mfa origins.
cat output/intersectSubtractedDnascent-PredominantsWithOrc1cdc6.bed | cut -f5-7 | sort -u \
  > output/orc1Cdc6ThatMakeUpFlexibleOrigins.bed

# Generating orc1cdc6 coordinates that overlap with mfa origins.-NotUniq
cat output/intersectSubtractedDnascent-PredominantsWithOrc1cdc6.bed | cut -f5-7 \
  > output/orc1Cdc6ThatMakeUpFlexibleOrigins-NotUniq.bed

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

