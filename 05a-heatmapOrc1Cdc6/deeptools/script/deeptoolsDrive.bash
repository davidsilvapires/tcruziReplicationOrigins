#!/usr/bin/env bash

###################################################################################################
#                                       deeptoolsDriver                                           #
#                                                                                                 #
# This script will direct the other pre-processing, processing, and post-processing scripts.      #
# Usage: saveCommand script/deeptoolsDriver.bash                                                  #
#                                                                                                 #
# Copyleft (ɔ) 2023 by Thiago Andrade Franco                                                      #
#                     <thiago.franco.esib@esib.butantan.gov.br>                                   #
#    12/07/2023: First version.                                                                   #
###################################################################################################

# Loading environment variables about this experiment.
. script/deeptoolsVariables.bash

# Running genomeComparmentComposition pre-processing, main and post-processing scripts.
saveCommand script/deeptoolsPre.bash
saveCommand script/deeptoolsMain.bash

exit 0
