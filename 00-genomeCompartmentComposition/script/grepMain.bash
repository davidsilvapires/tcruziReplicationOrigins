#!/usr/bin/env bash

###################################################################################################
#                                         grepMain                                                #  
#                                                                                                 #
#  This script allows the creation of databases of the compartments of the Trypanosoma cruzi      #
# genome according to the information found in the literature. This script has been modified      #
# to work  with the genome version 32. Minor code changes are required for the most recent        #
# versions. It writes the Slurm job script and submits all of them.                               #
#                                                                                                 #
# Usage: saveCommand script/grepMain.bash                                                         #
#                                                                                                 #
# Copyleft (ɔ) 2023 by Thiago Andrade Franco <thiago.franco.esib@esib.butantan.gov.br>            #
#    12/07/2023: First version.                                                                   #
###################################################################################################

# Writing the job script to submits at Slurm
cat> job/grepMain.slurm <<EOI
#!/usr/bin/env bash

# In the first step, we will create a gff file containing the coding sequence for all analyses.
# To isolate only the coding sequence region, we will need to select only lines with "gene" in the third column.
grep -P '.+\t.*\tgene\t.*\t.*\t.*\t.*\t.*\t.*' \
  input/TriTrypDB-32_TcruziCLBrenerEsmeraldo-like.gff \
  > output/TriTrypDB-32_TcruziCLBrenerEsmeraldo-like-gene.gff

# In this step, we will create the files of each gene of the multigene family
for FEATURE in  'DGF|dispersed' 'GP63|protease' 'RHS|retrotransposon' 'TS|trans-sialidase'
do
  /usr/bin/nice -n 19 grep -Pi $'.+\t.*\t.*\t.*\t.*\t.*\t.*\t.*\t.*'\${FEATURE}$'.*' \ 
    output/TriTrypDB-32_TcruziCLBrenerEsmeraldo-like-gene.gff \ 
    > output/TriTrypDB-32_TcruziCLBrenerEsmeraldo-like-gene-\${FEATURE_TAG}.gff 
done

for FEATURE in 'Mucin-associated'  'mucin-associated' 
do
  /usr/bin/nice -n 19 grep -P $'.+\t.*\t.*\t.*\t.*\t.*\t.*\t.*\t.*'\${FEATURE}$'.*' \ 
    output/TriTrypDB-32_TcruziCLBrenerEsmeraldo-like-gene.gff \ 
    > output/TriTrypDB-32_TcruziCLBrenerEsmeraldo-like-gene-\${FEATURE}.gff 
done

grep -P 'mucin$' output/TriTrypDB-32_TcruziCLBrenerEsmeraldo-like-gene.gff \
> output/TriTrypDB-32_TcruziCLBrenerEsmeraldo-like-gene-mucin.gff

wait

# To create the database for the Core compartment, we must first extract all gene sequences from multigene families.
# As a result, we remove all the genes from the multigene families from the gff file assembled in the previous step,
# which only contains the coding region (CDS) .
  grep -v $'.*\t.*\t.*\t.*\t.*\t.*\t.*\t.*\t.*'Mucin-associated$'.*' \
    output/TriTrypDB-32_TcruziCLBrenerEsmeraldo-like-gene.gff|
  grep -v $'.*\t.*\t.*\t.*\t.*\t.*\t.*\t.*\t.*'mucin$'.*' |
  grep -v $'.*\t.*\t.*\t.*\t.*\t.*\t.*\t.*\t.*'Mucin$'.*' |
  grep -v $'.*\t.*\t.*\t.*\t.*\t.*\t.*\t.*\t.*'dispersed$'.*' |
  grep -v $'.*\t.*\t.*\t.*\t.*\t.*\t.*\t.*\t.*'protease$'.*' |
  grep -v $'.*\t.*\t.*\t.*\t.*\t.*\t.*\t.*\t.*'retrotransposon$'.*' |
  grep -v $'.*\t.*\t.*\t.*\t.*\t.*\t.*\t.*\t.*'trans-sialidase$'.*' \
  > output/TriTrypDB-32_TcruziCLBrenerEsmeraldo-like-coreCompartment.gff

# The Disruptive compartment, according to the literature, is made up of mucin, transialidases, and mucin-associated genes.
# As a result, we combined the files from these genes into a single file to create the disruptive compartment database.
cat output/TriTrypDB-32_TcruziCLBrenerEsmeraldo-like-gene-TS.gff \
    output/TriTrypDB-32_TcruziCLBrenerEsmeraldo-like-gene-mucin-associated.gff \
    output/TriTrypDB-32_TcruziCLBrenerEsmeraldo-like-gene-mucin.gff \
    output/TriTrypDB-32_TcruziCLBrenerEsmeraldo-like-gene-Mucin-associated.gff \
    > output/TriTrypDB-32_TcruziCLBrenerEsmeraldo-like-disruptiveCompartment.gff

# According to the literature, DGF-1, GP63, and retrotransposon genes can be found in both core and disruptive compartments.
# We created a third compartment called both that contained the sequences of these three distinct classes of genes for the purpose of our research.
cat output/TriTrypDB-32_TcruziCLBrenerEsmeraldo-like-gene-GP63.gff \
    output/TriTrypDB-32_TcruziCLBrenerEsmeraldo-like-gene-DGF.gff \
    output/TriTrypDB-32_TcruziCLBrenerEsmeraldo-like-gene-RHS.gff \
    > output/TriTrypDB-32_TcruziCLBrenerEsmeraldo-like-bothCompartment.gff

exit 0

EOI

# Submitting the job to Slurm.
/usr/bin/nice -n 19 /usr/bin/time --verbose --output=log/sbatch-grepMain.time memusg \
  --output-to-file log/sbatch-bedtoolsMain.memusg --time --shell "saveCommand sbatch --nodes 1 \
  --ntasks ${NUM_OF_CPUS} --mem ${MEMORY_SIZE} -o log/slurm-%A.out -J grepMain job/grepMain.slurm \
  2> log/sbatch-grepMain.err | tee log/sbatch-grepMain.out"

# Check how are going all the running awk processes. 
 echo "To check the runnning processes, execute one of the following commands:"
 echo "   $> watch -n 1 sequeue"

 exit 0
