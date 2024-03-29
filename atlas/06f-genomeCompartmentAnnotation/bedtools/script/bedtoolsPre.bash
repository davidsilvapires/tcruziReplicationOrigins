#!/usr/bin/env bash

###################################################################################################
#                                  grepPre                                                        #
#                                                                                                 #
# This script is a required pre-processing step to composing the Predominants origins  database.  #
# It sets up the directory structure and creates symbolic links for input data.                   #                                                                                                 
#                                                                                                 #
# Usage: saveCommand script/bedtoolsPre.bash                                                      #
#                                                                                                 #
# Copyleft (ɔ) 2023 by Thiago Andrade Franco  <thiago.franco.esib@esib.butantan.gov.br>           #
#    12/07/2023: First version.                                                                   #
###################################################################################################

# Setting up the initial directories structures.
mkdir final input job output 

# Storing initial processing date and time.
date > log/processingStart.date

# Making a symbolic link for input data 

# to origins from ChIP-seq , MFA-seq and D-NAscent     
  ln -s /project/carol/dnascent/project/atlas/metaData/orc1cdc6.bed \
        input/orc1cdc6.bed
  ln -s /project/carol/dnascent/project/atlas/metaData/origens_mfa_3kb.bed \
        input/mfa.bed
  ln -s /project/carol/dnascent/project/atlas/metaData/dnascentNonSinc3kb.bed \
       input/dnascent.bed
 
# to atlas origins  
  ln -s /project/carol/dnascent/project/atlas/01-predominantOrigins/bedtools/\
final/predominantsOrigins.bed input/predominantsOrigins.bed
  ln -s /project/carol/dnascent/project/atlas/02-flexibleOrigins/bedtools/\
final/flexibleOrigins.bed input/flexibleOrigins.bed
  ln -s /project/carol/dnascent/project/atlas/03-dormantOrigins/bedtools/\
final/dormantOrigins.bed input/dormantOrigins.bed
  ln -s /project/carol/dnascent/project/atlas/04-orc1cdc6FreeOrigins/bedtools/\
final/orc1Cdc6-free.bed   input/orc1Cdc6-free.bed
  ln -s /project/carol/dnascent/project/atlas/05-totalOrigins/\
final/totalOrigins.bed input/atlasOrigins.bed

# to genome compartment      
  ln -s /project/carol/dnascent/project/atlas/00-genomeCompartmentComposition/\
release-32/final/TriTrypDB-32_TcruziCLBrenerEsmeraldo-like-coreCompartment.gff \
input/core.gff
  ln -s /project/carol/dnascent/project/atlas/00-genomeCompartmentComposition/\
release-32/final/TriTrypDB-32_TcruziCLBrenerEsmeraldo-like-disruptiveCompartment.gff \
input/disruptive.gff
  ln -s /project/carol/dnascent/project/atlas/00-genomeCompartmentComposition/\
release-32/final/TriTrypDB-32_TcruziCLBrenerEsmeraldo-like-bothCompartment.gff \
input/both.gff

# Writing the command lines necessary to run this pipeline.
cat > cmdLine.bash << EOI
# Copying the template scripts.
cp -a /project/carol/dnascent/pipeline/atlas/01-predominantOrigins/bedtools/* .

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
