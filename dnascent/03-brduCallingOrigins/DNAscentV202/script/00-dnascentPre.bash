# Setting up the initial directories structures.
mkdir final input job log output

# Indexing the raw reads and mapped bam files.
RAW_READS="../../../../../data/miniON/brdu/tryCru/brdu2/final/"
SEQ_SUMMARY="../../../00-baseCalling/guppy/brdu2/final/sequencing_summary.txt"

# Making symbolic links for reference genome and input data.

ln -s ../../../../../../orgn/tryCru-clb/tryCru-clb1/genome/final/tryCru-clb1.fasta input/

for FILE in ../../../02-samToIndexedBam/samtools/brdu2/final/*.bam*
do
  ln -s ../${FILE} input/
done

# Writing the job script.
cat > job/dnascent1.job << EOI
#!/usr/bin/env bash

# Making DNAsecent index

/usr/bin/nice -19 /usr/bin/time --verbose --output=log/index.time saveCommand DNAscent index -f ${RAW_READS} -s ${SEQ_SUMMARY} -o output/index.dnascent > log/index.out 2> log/index.err &
wait

# Detection of BrdU analog

for FILE in \`ls input | grep '.bam$'\`
do
  /usr/bin/nice -n 19 /usr/bin/time --verbose --output=log/\${FILE}.brduDetect.time saveCommand DNAscent detect --threads 80 --quality 20 --length 1000 --bam input/\${FILE} --reference input/tryCru-clb1.fasta --index output/index.dnascent --output output/\${FILE}.brduDetect > log/\${FILE}.brduDetect.out 2> log/\${FILE}.brduDetect.err 
done

# Call regions of analogue incorporation

#for FILE in \`ls output | grep '.brduDetect$'\`
#do
#  /usr/bin/nice -n 19 /usr/bin/time --verbose --output=log/\${FILE}.time saveCommand DNAscent regions -d output/\${FILE} -o output/\${FILE}.regions > log/\${FILE}.regions.out 2> log/\${FILE}.regions.err &
#done

# Estimating fork direction

#for FILE in \`ls output | grep '.brduDetect$'\`
#do
#  /usr/bin/nice -n 19 /usr/bin/time --verbose --output=log/\${FILE}.time saveCommand DNAscent forkSense --markOrigins --markTerminations -t 10 -d output/\${FILE} -o output/\${FILE}.forkSense > log/\${FILE}.forkSense.out 2> log/\${FILE}.forkSense.err &
#done

wait

exit 0
EOI

# Composing the sbatch command.

saveCommand sbatch -n80 --mem=200G -w "vital" -o log/sbatch-brduDetect.log job/dnascent1.job 2> log/sbatch-brduDetect.err

