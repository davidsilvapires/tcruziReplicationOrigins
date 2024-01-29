#!/usr/bin/env bash

###################################################################################################
#                                  totalOriginsComposition                                        #
#                                                                                                 #
# This script is a required pre-processing step to composing the Total origins  dataset           #   
# It sets up the directory structure and creates symbolic links for input data.                   #
#                                                                                                 #
# Usage: saveCommand script/PretotalOrigins.bash                                                  #
#                                                                                                 #
# Copyleft (ɔ) 2023 by  Thiago Andrade Franco                                                     #
#                       <thiago.franco.esib@esib.butantan.gov.br>                                 # 
#    11/10/2023: First version.                                                                   #
###################################################################################################

# Setting up the initial directories structures.
mkdir final input job output 

# Storing initial processing date and time.
date > log/processingStart.date

# Making a symbolic link for input data 
  ln -s /project/carol/dnascent/project/atlas/01-predominantOrigins/bedtools/\
final/predominantsOrigins.bed input/
  ln -s /project/carol/dnascent/project/atlas/02-flexibleOrigins/bedtools/\
final/flexibleOrigins.bed input/
  ln -s /project/carol/dnascent/project/atlas/03-dormantOrigins/bedtools/\
final/dormantOrigins.bed input/

# Writing the command lines necessary to run this pipeline.
cat > cmdLine.bash << EOI
# Copying the template scripts.
cp -a /project/carol/dnascent/pipeline/atlas/05-totalOrigins/* .

# Running bedtoolsDriver script.
saveCommand script/totalOriginsCompositionDriver.bash\
  2>&1 | tee log/totalOriginsCompositionDriver.out

# Checking result.
saveCommand script/totalOriginsCompositionCheck.bash\
  2>&1 | tee log/totalOriginsCompositionCheck.out

# Removing intermediate files.
saveCommand script/totalOriginsCompositionClean.bash\
  2>&1 | tee log/totalOriginsCompositionClean.out

EOI

exit 0
