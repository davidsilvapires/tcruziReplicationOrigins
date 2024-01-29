###################################################################################################
#                                          grepVariables                                          #
#                                                                                                 #
# This script contains all the variables needed in this computational analysis.                   #
#                                                                                                 #
# Usage: . script/bedtoolsVariables.bash                                                          #
#                                                                                                 #
# Copyleft (ɔ) 2023 by Thiago Andrade Franco     and David da Silva Pires                         #
#                     <thiago.franco.esib@esib.butantan.gov.br> and <david.pires@butantan.gov.br  #    
#    12/12/2023: First version.                                                                   #
###################################################################################################

# Number of Memery necessary to slurm
export MEMORY_SIZE="1G"

# Number of CPUS necessary to rum job at slurm
export NUM_OF_CPUS=1

# genome chromsize
export GENOME_CHROMSIZE="/project/carol/dnascent/orgn/tryCru-clb/tryCru-clb5/genome/final/tryCru-clb5-chromSizes.txt"

# DNA replication origins analized
export ORIGIN_TYPES="RHS orc1cdc6"





