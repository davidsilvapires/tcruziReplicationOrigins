#!/usr/bin/env bash

# hisat2 pre-processing pipeline.

# Example of usage:
#    $> ./hisat2Pre.bash

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

# Checking if there is a new version.
#echo "Your current bowtie2 version is:"
#bowtie2 --version | head -1 | cut -d' ' -f3
#echo
#echo "Please, verify if there is a new version available for bowtie2 at the following URL:"
#echo "https://sourceforge.net/projects/bowtie-bio/files/bowtie2/"
#echo "If a new version is available, please, install it and read the changelog at the following URL:"
#echo "http://bowtie-bio.sourceforge.net/bowtie2/index.shtml"
#echo
#echo "Press <ENTER> to continue or <CTRL + C> to stop."
#read

# Setting up the initial directories structure.
mkdir final input job log output

# Making symbolic links for input data files.
# The input field separator has to be '\n' such that we can parse each line instead of each word.
IFS='
'
for i in `head -n -1 ${METADATA}`; do
  # We take fields 1, 2 and 3, corresponding to ID, R1 mate and R2 mate, respectively.
  ID=`echo ${i} | cut -d' ' -f1`
  R1_MATE=`echo ${i} | cut -d' ' -f2`
  R2_MATE=`echo ${i} | cut -d' ' -f3`
  
	# The ".." is necessary because we are creating a link inside a subdirectory.
	# From now on, the files are numbered according to metadata samples file unique ID at first
  # column. This favors the use of job arrays through the use of environment variables.
	ln -s ../../../../../../../../00-readFiltering/leta/marcelaVitarelli/01-adapterTrimming/trimmomatic/final/${R1_MATE} input/${ID}_R1.fastq
	ln -s ../../../../../../../../00-readFiltering/leta/marcelaVitarelli/01-adapterTrimming/trimmomatic/final/${R2_MATE} input/${ID}_R2.fastq
done

# Making a symbolic link for bowtie2 index for T. cruzi CLBrener Esmeraldo like genome version 32.
ln -s /project/carol/chiporcs/orgn/tryCru-clb/tryCru-clb1/genome/index/hisat2/final input/index

# Writing the job script.
# -q reads are FASTQ files.


cat > job/hisat2.job << EOI
#!/usr/bin/env bash

for FILE in \`seq 1 ${NUM_SAMPLES}\`
do
  /usr/bin/nice -n 19 /usr/bin/time --verbose --output=log/\${FILE}.time saveCommand hisat2 --threads ${CPU_PER_PROC} -q --phred33 --n-ceil L,0,0.15 --score-min L,-0.6,-0.6 --no-spliced-alignment --no-unal -x input/index/tryCru-clb1 -1 input/\${FILE}_R1.fastq -2 input/\${FILE}_R2.fastq -S output/\${FILE}.sam > log/\${FILE}.out 2> log/\${FILE}.err &
done
wait

exit 0
EOI

# Editing the job script to confirm that everything is OK.
vim job/hisat2.job

# Composing the sbatch and watch commands.
QSUB_COMMAND="saveCommand sbatch -n96 --mem=400G -w \"vital\" -o log/sbatch.log job/hisat2.job 2> log/sbatch.err | tee log/sbatch.out"
WATCH_COMMAND="watch -n 1 squeue"

# Explaining how to run the job script.
echo "If you are happy about the contents of the job script that you just edited, you can run it by doing:"
echo "${QSUB_COMMAND}"
echo
echo "To see the status of your job array, do:"
echo "${WATCH_COMMAND}"
echo
echo "Please, after all the jobs are done, consider also running the script bowtie2Post.bash."

# Composing an initial cmdLine.bash file.
GROUP_FOLDER=`echo $(pwd) | sed -e 's|/project/chIPSeq.*||'`
{
  echo "# Copying the template for bowtie2 tool (pre-processing step)."
  echo "cp ${GROUP_FOLDER}/pipeline/template/chIPSeq/00-mapping/bowtie2/bowtie2Pre.bash script/"
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
