#!/usr/bin/env bash

###################################################################################################
#                                  grepPre                                                        #
#                                                                                                 #
# This script is a required pre-processing step to composing the Genome Compartment database.     #
# It sets up the directory structure and creates symbolic links for input data.                   #                                                                                                 #
#                                                                                                 #
# Usage: saveCommand script/grepPre.bash                                                          #
#                                                                                                 #
# Copyleft (ɔ) 2023 by Thiago Andrade Franco   <thiago.franco.esib@esib.butantan.gov.br>          #
#    12/07/2023: First version.                                                                   #
###################################################################################################

# Setting up the initial directories structures.
mkdir final input job output 

# Storing initial processing date and time.
date > log/processingStart.date

# Making a symbolic link for input data 
      ln -s /project/carol/dnascent/orgn/tryCru-clb/tryCru-clb1/gff/final/tryCru-clb1.gff \
        input/TriTrypDB-32_TcruziCLBrenerEsmeraldo-like.gff

# Writing the command lines necessary to run this pipeline.
cat > cmdLine.bash << EOI
# Copying the template scripts.
cp -a /project/carol/dnascent/pipeline/atlas/00-genomeCompartmentComposition/* .

# Running grepdriver script.
saveCommand script/grepDriver.bash\
  2>&1 | tee log/grepDriver.out

# Checking result.
saveCommand script/grepCheck.bash\
  2>&1 | tee log/grepCheck.out

# Removing intermediate files.
saveCommand script/grepClean.bash\
  2>&1 | tee log/grepClean.out
EOI

exit 0
