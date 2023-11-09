###################################################################################################
#                                          bedtoolsVariables                                      #
#                                                                                                 #
# This script contains all the variables needed in this computational analysis.                   #
#                                                                                                 #
# Usage: . script/deeptoolsVariables.bash                                                          #
#                                                                                                 #
# Copyleft (ɔ) 2023 by Thiago Andrade Franco                                                      #
#                     <thiago.franco.esib@esib.butantan.gov.br>                                   #    
#    11/07/2023: First version.                                                                   #
###################################################################################################

# Number of Memery necessary to slurm
export MEMORY_SIZE="30G"

# Number of CPUS necessary to rum job at slurm
export NUM_OF_CPUS=30

# Genome fasta file
export GENOME="input/tryCru-clb1.fa"

# Genomic coordinates (bed file ) that will be converted into bigwig
export BEDFILE="input/totalNucleosome.bed"
