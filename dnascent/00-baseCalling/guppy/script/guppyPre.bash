# Setting up the initial directories structures.
mkdir final input job log output

# Making symbolic links for input data files.

for FILE in ../../../../../data/miniON/brdu/tryCru/brdu2/final/*.fast5
do 
   ln -s ../${FILE} input/
done

# Writing the job script.
cat > job/guppy.job << EOI
#!/usr/bin/env bash

   /usr/bin/nice -n 0 /usr/bin/time --verbose --output=log/guppy.time saveCommand guppy_basecaller -i input/ -s output/ --flowcell FLO-MIN106 --kit SQK-LSK109 --cpu_threads_per_caller 240 --num_callers 6 --qscore_filtering --trim_strategy dna > log/guppy.out 2> log/guppy.err

exit 0
EOI

# Composing the sbatch command.

saveCommand sbatch -n80 --mem=200G -w "vital" -o log/sbatch.log job/guppy.job 2> log/sbatch.err

