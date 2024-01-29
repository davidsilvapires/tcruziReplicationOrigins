#!/usr/bin/env bash

###################################################################################################
#                                         bedtoolsMain                                            #  
#                                                                                                 #
#   This script allows the creation of Dormant origins database. In summary, the Dormant origins  #
# database is made up of the genomic coordinates of Orc1Cdc6 ChIP-Seq peaks that do not overlap   #
# with the genomic coordinates of the predominant and flexible origins.In this manner, we utilized#
# the bedtools subtract to generate the dormant origins database. It writes the Slurm job script  #
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

#Sorting bed files by chromosomes and by starting position.
for i in \`ls input/\`
do 
 sort -k1,1 -k2,2n input/\${i} > output/\${i%.bed}-sorted.bed
done

# In the first step, we will eliminate all Orc1Cdc6 genomic coordinates that overlap with the
# genomic coordinates of the predominant origins.
bedtools subtract -A -a output/orc1cdc6-sorted.bed \
  -b output/orc1Cdc6ThatMakeUpPredominantOrigins-sorted.bed \
  >  output/subtractOrc1cdc6-Orc1Cdc6ThatMakeUpPrediminantOrigins.bed \
  2> log/subtractOrc1cdc6-Orc1Cdc6ThatMakeUpPrediminantOrigins.err

# Then, we will remove all Orc1Cdc6 genomic coordinates that overlap with the genomic coordinates of
# the flexible origins
bedtools subtract -A -a output/subtractOrc1cdc6-Orc1Cdc6ThatMakeUpPrediminantOrigins.bed \
  -b output/orc1Cdc6ThatMakeUpFlexibleOrigins-sorted.bed \
  >  output/subtractOrc1Cdc6-Orc1Cdc6ThatMakeUP-PredominantAndFlexibleOrigins.bed \
  2>  log/subtractOrc1Cdc6-Orc1Cdc6ThatMakeUP-PredominantAndFlexibleOrigins.err

  # Generating dormant origins file.
  cp  output/subtractOrc1Cdc6-Orc1Cdc6ThatMakeUP-PredominantAndFlexibleOrigins.bed \
    output/dormantOrigins.bed

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

