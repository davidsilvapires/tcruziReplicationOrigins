# Setting up the initial directories structures.
mkdir final input job log output 

LAST_FILE_NUMBER=`ls ../../03-brduCallingOrigins/DNAscent/brdu/final/*_*_*.forkSense | cut -d_ -f4 | sort -V | tail -n1`

for i in `seq 0 ${LAST_FILE_NUMBER}`
do
  mkdir input/${i}

  for FILE in ../../03-brduCallingOrigins/DNAscent/brdu/final/*_${i}_*.forkSense
  do
    ln -s ../../${FILE} input/${i}
  done

  # Concatenating files.
  cat input/${i}/* > input/${i}/allReads.forkSense
done

# Writing the job script.
cat > job/awk.bash << EOI
#!/usr/bin/env bash

LAST_FILE_NUMBER=`ls ../../03-brduCallingOrigins/DNAscent/brdu/final/*_*_*.forkSense | cut -d_ -f4 | sort -V | tail -n1`

for INDEX in \`seq 0 ${LAST_FILE_NUMBER}\`
do

  for FILE in \`ls input/\${INDEX}/*.forkSense\`
  do
  
     # To identify the beginning and name of each read and redirecting it to a new file. Changing
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
done &
wait

exit 0
EOI

# Composing the sbatch command.
saveCommand sbatch -n34 --mem=340G -w "vital" -o log/sbatch.log job/awk.bash 2> log/sbatch.err | tee log/sbatch.out
