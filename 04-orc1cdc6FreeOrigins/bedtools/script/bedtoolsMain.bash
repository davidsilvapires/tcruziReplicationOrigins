#!/usr/bin/env bash

###################################################################################################
#                                         bedtoolsMain                                            #  
#                                                                                                 #
#  This script allows the creation of Orc1Cdc6-free origins database.In summary, the Orc1Cdc6-free#
# origins database is made up of the genomic coordinates of DNA replication origins discovered by #
# D-nascent that do not overlap with the coordinates of Orc1Cdc6 peaks.In this manner, we utilized#
# the bedtools subtract to generate the Orc1Cdc6-free origins database.It writes the Slurm job    #
# script and submits all of them.                                                                 #
#                                                                                                 #
# Usage: saveCommand script/bedtoolsMain.bash                                                     #
#                                                                                                 #
# Copyleft (ɔ) 2023 by Thiago Andrade Franco                                                      #
#                   < thiago.franco.esib@esib.butantan.gov.br>                                    #
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
-b output/orc1cdc6-sorted.bed \
> output/subtractDnascent-orc1cdc6.bed \
2>log/subtractDnascent-orc1cdc6.err 

# Generating flexible origins file.
cp output/subtractDnascent-orc1cdc6.bed output/orc1Cdc6-free.bed

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
