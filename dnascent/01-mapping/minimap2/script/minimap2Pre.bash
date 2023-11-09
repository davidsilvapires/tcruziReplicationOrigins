# minimap2 -ax map-ont -t 24 referenceGenome.fasta baseCalledReads.fastq (or .fasta) -o output.sam
# -a : output in SAM format
# -x : tuned for optimal performance and accuracy
# map-ont : uses ordinary minimizers as seeds. Use this seeting for Oxfore Nanopore reads.
# -t : number of threads
# Note: there is no need to index the genome, unless it's the human one.

# Setting up the initial directories structures.
mkdir final input job log output

# Making symbolic links for input data files.

for FILE in ../../../00-baseCalling/guppy/brdu2/final/*.fastq
do
  ln -s ../${FILE} input/
done

# Writing the job script.
cat > job/minimap2.job << EOI
#!/usr/bin/env bash

for FILE in \`ls input/\`
do

/usr/bin/nice -n 19 /usr/bin/time --verbose --output=log/minimap2.time saveCommand minimap2 -ax map-ont -t 80 /project/carol/dnascent/orgn/tryCru-clb/tryCru-clb1/genome/final/tryCru-clb1.fasta input/\${FILE} -o output/\${FILE}aln.sam > log/\${FILE}.out 2> log/\${FILE}.err &

done
wait

exit 0
EOI

# Composing the sbatch command.

saveCommand sbatch -n80 --mem=300G -w "vital" -o log/sbatch.log job/minimap2.job 2> log/sbatch.err
