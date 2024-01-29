#!/usr/bin/env bash

###################################################################################################
#                                 totalOriginsCompositionDriver                                   #
#                                                                                                 #
# This script will direct the other pre-processing, processing, and post-processing scripts.      #
# Usage: saveCommand script/totalOriginsCompositionDriver.bash                                    #
#                                                                                                 #
# Copyleft (ɔ) 2023 by Thiago Andrade Franco                                                      #
#                     <thiago.franco.esib@esib.butantan.gov.br>                                   #
#    11/10/2023: First version.                                                                   #
###################################################################################################

# Loading environment variables about this experiment.
. script/totalOriginsCompositionVariables.bash

# Running genomeComparmentComposition pre-processing, main and post-processing scripts.
saveCommand script/totalOriginsCompositionPre.bash
saveCommand script/totalOriginsCompositionMain.bash

exit 0
