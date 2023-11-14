#!/usr/bin/env bash

###################################################################################################
#                                           grepDriver                                            #
#                                                                                                 #
# This script will direct the other pre-processing, processing, and post-processing scripts.      #
# Usage: saveCommand script/grepDriver.bash                                                       #
#                                                                                                 #
# Copyleft (ɔ) 2023 by Thiago Andrade Franco  <thiago.franco.esib@esib.butantan.gov.br>           #
#    12/07/2023: First version.                                                                   #
###################################################################################################

# Loading environment variables about this experiment.
. script/grepVariables.bash

# Running genomeComparmentComposition pre-processing, main and post-processing scripts.
saveCommand script/grepPre.bash
saveCommand script/grepMain.bash

exit 0
