# Writing the job script.
cat > job/dnascent3.job << EOI
#!/usr/bin/env bash

# Estimating fork direction

for FILE in \`ls output | grep '.brduDetect$'\`
do
  /usr/bin/nice -n 19 /usr/bin/time --verbose --output=log/\${FILE}.forkSense.time saveCommand DNAscent forkSense --threads 40 --markOrigins --markTerminations --markForks -d output/\${FILE} -o output/\${FILE}.forkSense > log/\${FILE}.forkSense.out 2> log/\${FILE}.forkSense.err 
done

wait

exit 0
EOI

# Composing the sbatch command.

saveCommand sbatch -n40 --mem=200G -w "vital" -o log/sbatch-forkSense.log job/dnascent3.job 2> log/sbatch-forkSense.err

