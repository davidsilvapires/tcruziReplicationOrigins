# Setting up the initial directories structures.
mkdir final input job log output

# Making symbolic links for input data.
for FILE in `ls ../03-brduCallingOrigins/DNAscent/final/*.regions`
do
  ln -s ../${FILE} input/
done

# Concatenating files.
#cat input/fastq_*.regions >> input/allReads.regions

# Writing the job script.
cat > job/awk.bash << EOI
#!/usr/bin/env bash

for FILE in \`ls input/*.regions\`
do

   # To identify the beginning and name of each read and redirecting it to a new file. Changig
   # spaces to underscor in reads names and making a new file with start of the read and read name.
   < \${FILE} grep -n '^>' | cut -d':' -f1 > \${FILE}.readStart.txt
   < \${FILE} grep -n '^>' | cut -d':' -f2 > \${FILE}.readNames.txt
   sed -i 's/ /_/g' \${FILE}.readNames.txt
   paste \${FILE}.readStart.txt \${FILE}.readNames.txt > \${FILE}.readStartNames.txt

   # To indentify the last line of the file and adding in the end of the readStartNames file. Adding
   # one extra line to get the last read and adding 0 to identify the end of the file.
   NUM_OF_LINES=\`cat \${FILE} | wc -l\`
   let "NUM_OF_LINES = NUM_OF_LINES + 1"
   echo \${NUM_OF_LINES} >> \${FILE}.readStartNames.txt
   echo "0" >> \${FILE}.readStartNames.txt

   LINE=\`head -n 1 \${FILE}.readStart.txt\`
   READ_NAME=\`head -n 1 \${FILE}.readNames.txt\`
   STARTS=\`tail -n +2 \${FILE}.readStartNames.txt\`

   # To show read by read.
   IFS='
'
   for i in \${STARTS} 
   do
      echo "i = \${i}"
      read
      NEXT_LINE=\`echo \${i} | cut -f1\`
      if [[ NEXT_LINE == 0 ]]
      then
         break
      fi
      < \${FILE} awk -v LINE="\${LINE}" -v NEXT_LINE="\${NEXT_LINE}" 'NR==LINE, NR==NEXT_LINE-1 {print}' > output/"\${READ_NAME}.txt"
      let "LINE = NEXT_LINE"
      READ_NAME=\`echo \${i} | cut -f2\`
   done
done

exit 0
EOI

# Composing the sbatch command.

saveCommand sbatch -n30 --mem=30G -w "vital" -o log/sbatch.log job/awk.bash 2> log/sbatch.err 

