###################################################################################################
#                                          grepVariables                                          #
#                                                                                                 #
# This script contains all the variables needed in this computational analysis.                   #
#                                                                                                 #
# Usage: . script/bedtoolsVariables.bash                                                          #
#                                                                                                 #
# Copyleft (ɔ) 2023 by Marcela de Oliveira Vitarelli and Thiago Andrade Franco                    #
#                     <vitarelli.marcela@gmail.com> and <thiago.franco.esib@esib.butantan.gov.br> #    
#    11/07/2023: First version.                                                                   #
###################################################################################################

# Number of Memery necessary to slurm
export MEMORY_SIZE="10G"

# Number of CPUS necessary to rum job at slurm
export NUM_OF_CPUS=10

# Genome file that will be used in the annotation
export GENOME="input/genome.gff"








