###################################################################################################
#                                          grepVariables                                          #
#                                                                                                 #
# This script contains all the variables needed in this computational analysis.                   #
#                                                                                                 #
# Usage: . script/grepVariables.bash                                                              #
#                                                                                                 #
# Copyleft (ɔ) 2023 by Thiago Andrade Franco   <thiago.franco.esib@esib.butantan.gov.br>          #    
#    11/07/2023: First version.                                                                   #
###################################################################################################

# Number of Memery necessary to slurm
export MEMORY_SIZE="20G"

# Number of CPUS necessary to rum job at slurm
export NUM_OF_CPUS=10

# tag to file articles
export FEATURE_TAG=`echo ${FEATURE} | cut -d '|' -f1`






