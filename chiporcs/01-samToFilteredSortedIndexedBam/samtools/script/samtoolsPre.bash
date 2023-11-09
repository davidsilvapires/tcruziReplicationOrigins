#!/usr/bin/env bash

# samToSortedFilteredIndexedBam pre-processing pipeline.

# Example of usage:
#    $> ./samToSortedFilteredIndexedBam.bash

# Directory where we can find the data for this run.
DATA_DIR="../../metadata"

# File that contains metadata about this project.
METADATA="${DATA_DIR}/samples.txt"

# Computing how many samples we have for this run.
NUM_SAMPLES=`wc -l ${METADATA} | cut -d' ' -f1`
let "NUM_SAMPLES = NUM_SAMPLES - 1" # dummySample doesn't count as a valid input.

# Computing how many CPUs each process can use.
NUM_CPUS=120
let "CPU_PER_PROC = ${NUM_CPUS} / ${NUM_SAMPLES}"

# Minimum accepted mapping quality.
QUALITY=10

# Checking if there is a new version.
echo "Your current samtools version is:"
samtools --version | head -1 | cut -d' ' -f2
echo
echo "Please, verify if there is a new version available for samtools at the following URL:"
echo "https://sourceforge.net/projects/samtools/files/"
echo "If a new version is available, please, install it and read the changelog at the README.txt"
echo "file in the Summary tab at the following URL:"
echo "https://sourceforge.net/projects/samtools"
echo
echo "Press <ENTER> to continue or <CTRL + C> to stop."
read

# Setting up the initial directories structure.
mkdir final input job log output

# Making symbolic links for input data files.
# The input field separator has to be '\n' such that we can parse each line instead of each word.
IFS='
'
for ID in `seq 1 ${NUM_SAMPLES}`; do
  # The ".." is necessary because we are creating a link inside a subdirectory.
  # Remember that the files are numbered according to metadata samples file unique ID at first
  # column. This favors the use of job arrays through the use of environment variables.
  ln -s ../../../00-mapping/hisat2/final/${ID}.sam input/
done

# Writing the job script.
cat > job/samToFilteredSortedIndexedBam.job << EOI
#!/usr/bin/env bash

for FILE in \`seq 1 ${NUM_SAMPLES}\`; do
   /usr/bin/nice -n 19 /usr/bin/time --verbose --output=log/\${FILE}-filtered.time saveCommand samtools view -@ ${CPU_PER_PROC} -hSb -q ${QUALITY} input/\${FILE}.sam -o output/\${FILE}-filtered.bam > log/\${FILE}-filtered.out 2> log/\${FILE}-filtered.err && \\
   /usr/bin/nice -n 19 /usr/bin/time --verbose --output=log/\${FILE}-filteredSorted.time saveCommand samtools sort -@ ${CPU_PER_PROC} -o output/\${FILE}-filteredSorted.bam output/\${FILE}-filtered.bam 2> log/\${FILE}-filteredSorted.err > log/\${FILE}-filteredSorted.out && \\
   /usr/bin/nice -n 19 /usr/bin/time --verbose --output=log/\${FILE}-filteredSortedIndexed.time saveCommand samtools index -@ ${CPU_PER_PROC} output/\${FILE}-filteredSorted.bam > log/\${FILE}-filteredSortedIndexed.out 2> log/\${FILE}-filteredSortedIndexed.err &
done
wait

exit 0
EOI

# Editing the job script to confirm that everything is OK.
vim job/samToFilteredSortedIndexedBam.job

# Composing the sbatch and watch commands.
QSUB_COMMAND="saveCommand sbatch -n120 --mem=200G -w \"vital\" -o log/sbatch.log job/samToFilteredSortedIndexedBam.job 2> log/sbatch.err | tee log/sbatch.out"
WATCH_COMMAND="watch -n 1 squeue"

# Explaining how to run the job script.
echo "If you are happy about the contents of the job script that you just edited, you can run it by doing:"
echo "${QSUB_COMMAND}"
echo
echo "To see the status of your job array, do:"
echo "${WATCH_COMMAND}"
echo
echo "Please, after all the jobs are done, consider also running the script samToFilteredSortedIndexedBamPost.bash."

# Composing an initial cmdLine.bash file.
GROUP_FOLDER=`echo $(pwd) | sed -e 's|/project/chIPSeq.*||'`
{
  echo "# Copying the template for samToFilteredSortedIndexedBam tool (pre-processing step)."
  echo "cp ${GROUP_FOLDER}/pipeline/template/chIPSeq/01-samToFilteredSortedIndexedBam/samtools/samToFilteredSortedIndexedBamPre.bash script/"
  echo
  echo "# Executing the pre-processing template."
  echo "saveCommand ${0}"
  echo
  echo "# Running the job script."
  echo "${QSUB_COMMAND}"
  echo
  echo "# Seeing the status of the job array."
  echo "${WATCH_COMMAND}"
} > cmdLine.bash

exit 0
