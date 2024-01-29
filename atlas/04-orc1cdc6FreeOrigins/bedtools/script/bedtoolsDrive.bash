#!/usr/bin/env bash

###################################################################################################
#                                        bedtoolsDriver                                           #
#                                                                                                 #
# This script will direct the other pre-processing, processing, and post-processing scripts.      #
# Usage: saveCommand script/bedtoolsDriver.bash                                                   #
#                                                                                                 #
# Copyleft (ɔ) 2023 by Thiago Andrade Franco                                                      #
#                     <thiago.franco.esib@esib.butantan.gov.br>                                   #
#    12/07/2023: First version.                                                                   #
###################################################################################################

# Loading environment variables about this experiment.
. script/bedtoolsVariables.bash

# Running genomeComparmentComposition pre-processing, main and post-processing scripts.
saveCommand script/bedtoolsPre.bash
saveCommand script/bedtoolsMain.bash

exit 0
