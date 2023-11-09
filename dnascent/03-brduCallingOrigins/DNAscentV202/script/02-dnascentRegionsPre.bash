# Writing the job script.
cat > job/dnascent2.job << EOI
#!/usr/bin/env bash

# Call regions of analogue incorporation

for FILE in \`ls output | grep '.brduDetect$'\`
do
  /usr/bin/nice -n 19 /usr/bin/time --verbose --output=log/\${FILE}.regions.time saveCommand DNAscent regions -d output/\${FILE} -o output/\${FILE}.regions > log/\${FILE}.regions.out 2> log/\${FILE}.regions.err 
done

wait

exit 0
EOI

# Composing the sbatch command.

saveCommand sbatch -n40 --mem=200G -w "vital" -o log/sbatch-regions.log job/dnascent2.job 2> log/sbatch-regions.err

