###################################################################################################
#                                          awkVariables                                          #
#                                                                                                 #
# This script contains all the variables needed in this computational analysis.                   #
#                                                                                                 #
# Usage: . script/awkVariables.bash                                                          #
#                                                                                                 #
# Copyleft (ɔ) 2023 by  Thiago Andrade Franco     <thiago.franco.esib@esib.butantan.gov.br>       #    
#    11/07/2023: First version.                                                                   #
###################################################################################################

# Number of Memery necessary to slurm
export MEMORY_SIZE="10G"

# Number of CPUS necessary to rum job at slurm
export NUM_OF_CPUS=10

# The file contains the genomic coordinates of the DNA replication origins that were used in the analysis.
export FILENAME="input/predominantOrigins.bed"

# The file contains the genomic coordinates of the reference genome  that were used in the analysis.
export GFF_FILENAME="input/gff.gff"










