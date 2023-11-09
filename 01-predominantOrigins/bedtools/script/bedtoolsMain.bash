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

# Intersect function to identify overlapping peaks between MFA-origins with a +/-3kb window and
# Orc1Cdc6 peaks => Predominant origins.
 bedtools intersect -wo -a output/mfa-sorted.bed \
   -b output/orc1cdc6-sorted.bed \
   > output/intersectOriginsMfaWithOrc1cdc6.bed \
   2> log/intersectOriginsMfaWithOrc1cdc6.err

# We will create two origin databases, one with the coordinates of the Orc1Cdc6 that intersect with
# the origins from the MFA-Seq and the other with the coordinates of the MFA-Seq that intersect 
# with the coordinates of the Orc1Cdc6.

# As a result, we will use the steps below to select only the coordinates of the MFA-Seq origins
# that intersect with the Orc1cdc6 binding sites.
cat output/intersectOriginsMfaWithOrc1cdc6.bed| cut -f1-3 | sort -u \
  > output/predominantsOrigins.bed

# As a result, we will use the steps below to select only the coordinates of the MFA-Seq origins
# that intersect with the Orc1cdc6 binding sites.
cat output/intersectOriginsMfaWithOrc1cdc6.bed| cut -f1-3 \
  > output/predominantsOrigins-NotUniq.bed


# Now we use the step below to select only the coordinates of the Orc1Cdc6 that overlapping MFA
# origins.
cat output/intersectOriginsMfaWithOrc1cdc6.bed | cut -f4-6 | sort -u \
  > output/orc1Cdc6ThatMakeupPredominantOrigins.bed

# Now we use the step below to select only the coordinates of the Orc1Cdc6 that overlapping MFA
# origins. NotUniq
cat output/intersectOriginsMfaWithOrc1cdc6.bed | cut -f4-6 \
  > output/orc1Cdc6ThatMakeupPredominantOrigins-NotUniq.bed

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

