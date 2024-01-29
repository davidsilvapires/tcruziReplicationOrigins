#!/usr/bin/env bash

###################################################################################################
#                                  awkPre                                                        #
#                                                                                                 #
# This script is a required pre-processing step to composing the Predominants origins  database.  #
# It sets up the directory structure and creates symbolic links for input data.                   #                                                                                                 
#                                                                                                 #
# Usage: saveCommand script/awkPre.bash                                                      #
#                                                                                                 #
# Copyleft (ɔ) 2023 by Thiago Andrade Franco  <thiago.franco.esib@esib.butantan.gov.br>           #
#    12/07/2023: First version.                                                                   #
###################################################################################################

# Setting up the initial directories structures.
mkdir final input job output 

# Storing initial processing date and time.
date > log/processingStart.date

# Making a symbolic link for input data 
     ln -s /project/carol/dnascent/project/atlas/metaData/orc1cdc6.bed \
        input/orc1cdc6.bed
     ln -s /project/carol/dnascent/project/atlas/metaData/origens_mfa_3kb.bed \
        input/mfa.bed
     ln -s /project/carol/dnascent/project/atlas/metaData/dnascentNonSinc3kb.bed \
       input/dnascent.bed
     ln -s /project/carol/dnascent/project/atlas/01-predominantOrigins/bedtools/\
final/predominantsOrigins.bed input/predominantsOrigins.bed
     ln -s /project/carol/dnascent/project/atlas/02-flexibleOrigins/bedtools/\
final/flexibleOrigins.bed input/flexibleOrigins.bed
     ln -s /project/carol/dnascent/project/atlas/03-dormantOrigins/bedtools/\
final/dormantOrigins.bed input/dormantOrigins.bed
     ln -s /project/carol/dnascent/project/atlas/04-orc1cdc6FreeOrigins/bedtools/\

# Writing the command lines necessary to run this pipeline.
cat > cmdLine.bash << EOI
# Copying the template scripts.
cp -a /project/carol/dnascent/pipeline/atlas/06-originsAnnotation/awk/* .

# Running awkDriver script.
saveCommand script/awkDriver.bash\
  2>&1 | tee log/awkDriver.out

# Checking result.
saveCommand script/awkCheck.bash\
  2>&1 | tee log/awkCheck.out

# Removing intermediate files.
saveCommand script/awkClean.bash\
  2>&1 | tee log/awkClean.out
EOI

exit 0
