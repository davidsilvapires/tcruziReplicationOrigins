#!/usr/bin/env bash

###################################################################################################
#                                  grepPre                                                        #
#                                                                                                 #
# This script is a required pre-processing step to composing the Predominants origins  database.  #
# It sets up the directory structure and creates symbolic links for input data.                   #                                                                                                 #
#                                                                                                 #
# Usage: saveCommand script/bedtoolsPre.bash                                                      #
#                                                                                                 #
# Copyleft (ɔ) 2023 by Marcela de Oliveira Vitarelli and  Thiago Andrade Franco                   #
#                      <vitarelli.marcela@gmail.com> and <thiago.franco.esib@esib.butantan.gov.br>#
#    12/07/2023: First version.                                                                   #
###################################################################################################

# Setting up the initial directories structures.
mkdir final input job output 

# Storing initial processing date and time.
date > log/processingStart.date

# Making a symbolic link for input data 
  ## Genome used in this annotated 
 ln -s /project/carol/dnascent/project/atlas/metaData/tryCru-clb1.gff \
       input/genome.gff
 
  ## Origins identifed by mfa, ChIP-Seq and D-NAscent
     ln -s /project/carol/dnascent/project/atlas/metaData/orc1cdc6.bed \
        input/orc1cdc6.bed
     ln -s /project/carol/dnascent/project/atlas/metaData/origens_mfa_3kb.bed \
        input/mfa.bed
     ln -s /project/carol/dnascent/project/atlas/metaData/dnascentNonSinc3kb.bed \
       input/dnascent.bed

 ## Origins from atlas
 
     ln -s /project/carol/dnascent/project/atlas/01-predominantOrigins/bedtools/\
final/predominantsOrigins.bed input/predominantsOrigins.bed
     ln -s /project/carol/dnascent/project/atlas/02-flexibleOrigins/bedtools/\
final/flexibleOrigins.bed input/flexibleOrigins.bed
     ln -s /project/carol/dnascent/project/atlas/03-dormantOrigins/bedtools/\
final/dormantOrigins.bed input/dormantOrigins.bed
     ln -s /project/carol/dnascent/project/atlas/04-orc1cdc6FreeOrigins/bedtools/\
final/orc1Cdc6-free.bed   input/orc1Cdc6-free.bed
     ln -s /project/carol/dnascent/project/atlas/05-totalOrigins/\
final/totalOrigins.bed  input/totalOrigins.bed

# Writing the command lines necessary to run this pipeline.
cat > cmdLine.bash << EOI

# Copying the template scripts.
cp -a /project/carol/dnascent/pipeline/atlas/06d-originsAnnotation/bedtools/* .

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
