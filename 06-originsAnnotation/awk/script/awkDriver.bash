#!/usr/bin/env bash

###################################################################################################
#                                        awkDriver                                           #
#                                                                                                 #
# This script will direct the other pre-processing, processing, and post-processing scripts.      #
# Usage: saveCommand script/awkDriver.bash                                                   #
#                                                                                                 #
# Copyleft (ɔ) 2023 by Marcela de Oliveira Vitarelli and Thiago Andrade Franco                    #
#                    <vitarellimarcela@gmail.com>  and  <thiago.franco.esib@esib.butantan.gov.br> #
#    12/07/2023: First version.                                                                   #
###################################################################################################

# Loading environment variables about this experiment.
. script/awkVariables.bash

# Running genomeComparmentComposition pre-processing, main and post-processing scripts.
saveCommand script/awkPre.bash
saveCommand script/awkMain.bash

exit 0
