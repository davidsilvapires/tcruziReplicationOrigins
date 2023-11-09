# Setting up the initial directories structures.
mkdir final input job log output

# Making symbolic links for input data file.

for FILE in ../../../01-mapping/minimap2/brdu2/final/*
do
  ln -s ../${FILE} input/
done

# Writing the job script.
cat > job/samtools.job << EOI
#!/usr/bin/env bash 

for FILE in \`ls input/\`
do

/usr/bin/nice -n 19 /usr/bin/time --verbose --output=log/\${FILE}.bam.time saveCommand samtools view -@ 10 -hSb input/\${FILE} -o output/\${FILE}.bam > log/\${FILE}.bam.out 2> log/\${FILE}.bam.err && \\
/usr/bin/nice -n 19 /usr/bin/time --verbose --output=log/\${FILE}.sortedBam.time saveCommand samtools sort -@ 10 output/\${FILE}.bam -o output/\${FILE}.sortedBam > log/\${FILE}.sortedBam.out 2> log/\${FILE}.sortedBam.err && \\
/usr/bin/nice -n 19 /usr/bin/time --verbose --output=log/\${FILE}.indexedSortedBam.time saveCommand samtools index -@ 10 output/\${FILE}.sortedBam >log/\${FILE}.indexedSortedBam.out 2> log/\${FILE}.indexedSortedBam.err

done
wait

exit 0
EOI

# Composing the sbatch command.

saveCommand sbatch -n10 --mem=10G -w "vital" -o log/sbatch.log job/samtools.job 2> log/sbatch.err
