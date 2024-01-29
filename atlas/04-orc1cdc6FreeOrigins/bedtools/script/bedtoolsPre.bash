#!/usr/bin/env bash

###################################################################################################
#                                  bedtoolsPre                                                    #
#                                                                                                 #
# This script is a required pre-processing step to composing the Orc1Cdc6-free origins database.  #
# It sets up the directory structure and creates symbolic links for input data.                   #                                                                                                 #
#                                                                                                 #
# Usage: saveCommand script/bedtoolsPre.bash                                                      #
#                                                                                                 #
# Copyleft (ɔ) 2023 by Thiago Andrade Franco                                                      #
#                      <thiago.franco.esib@esib.butantan.gov.br>                                  #
#    12/07/2023: First version.                                                                   #
###################################################################################################

# Setting up the initial directories structures.
mkdir final input job output 

# Storing initial processing date and time.
date > log/processingStart.date

# Making a symbolic link for input data 
ln -s /project/carol/dnascent/project/atlas/metaData/orc1cdc6.bed\
  input/orc1cdc6.bed
ln -s /project/carol/dnascent/project/atlas/metaData/dnascentNonSinc3kb.bed\
  input/dnascent.bed

# Writing the command lines necessary to run this pipeline.
cat > cmdLine.bash << EOI
# Copying the template scripts.
cp -a /project/carol/dnascent/pipeline/atlas/04-orc1cdc6FreeOrigins/bedtools/* .

# Running bedtoolsDriver script.
saveCommand script/bedtoolsDriver.bash\
  2>&1 | tee log/bedtoolsDriver.out

# Checking result.
saveCommand script/bedtoolsCheck.bash\
  2>&1 | tee log/bedtoolsCheck.out

# Removing intermediate files.
saveCommand script/bedtoolsClean.bash\
  2>&1 | tee log/bedtoolsClean.out
EOI

exit 0
