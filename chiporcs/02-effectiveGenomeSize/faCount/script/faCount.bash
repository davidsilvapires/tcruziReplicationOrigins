# Calculating the effective genome size according with deepTools webpage tutorial
# (https://deeptools.readthedocs.io/en/develop/content/feature/effectiveGenomeSize.html). The EGS is
# defined as the length of the "mappable" genome and one way to calculate it, if multimapping reads
# are included in the analysis, is to remove the number of N bases in the genome.

# Setting up initial directories structure.
mkdir -p final input job log output

# Making symbolic links for input data files.

ln -s /project/carol/chiporcs/orgn/tryCru-clb/tryCru-clb1/genome/final/tryCru-clb1.fasta input/tryCru-clb1.fa

# Writing job script.
cat > job/faCount.job << EOI
#!/usr/bin/env bash

/usr/bin/nice -n 19 /usr/bin/time --verbose --output=log/faCount.time faCount input/tryCru-clb1.fa > output/tryCru-clb1-egs.txt 2> log/faCount.err

exit 0
EOI

# Composing sbatch command.
saveCommand sbatch -n1 --mem=300G -w "vital" -o log/sbatch.log job/faCount.job 2> log/sbatch.err | tee log/sbatch.out





