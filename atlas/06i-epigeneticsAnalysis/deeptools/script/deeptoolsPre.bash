#!/usr/bin/env bash

###################################################################################################
#                                  deeptoolsPre                                                   #
#                                                                                                 #
# This script is a required pre-processing step to plotHeatmap to Orc1cdc6 colocalization with    #
# histones.                                                                                       #
# It sets up the directory structure and creates symbolic links for input data.                   # 
#                                                                                                 #
# Usage: saveCommand script/deeptoolsPre.bash                                                     #
#                                                                                                 #
# Copyleft (ɔ) 2023 by Thiago Andrade Franco                                                      #
#                      <thiago.franco.esib@esib.butantan.gov.br>                                  #
#    12/15/2023: First version.                                                                   #
###################################################################################################

# Setting up the initial directories structures.
mkdir final input job output 

# Storing initial processing date and time.
date > log/processingStart.date

# Making a symbolic link for input data from ChIP-Seq peaks
ln -s /project/carol/dnascent/project/atlas/metaData/orc1cdc6.bed\
  input/orc1cdc6.bed

# Making a symbolic link for input data from Origins' atlas  
ln -s /project/carol/dnascent/project/atlas/01-predominantOrigins/bedtools/\
final/orc1Cdc6ThatMakeupPredominantOrigins-NotUniq.bed input/predominantOrigins.bed
ln -s /project/carol/dnascent/project/atlas/02-flexibleOrigins/bedtools/\
final/orc1Cdc6ThatMakeUpFlexibleOrigins-NotUniq.bed input/flexibleOrigins.bed
ln -s /project/carol/dnascent/project/atlas/03-dormantOrigins/bedtools/\
final/dormantOrigins.bed  input/dormantOrigins.bed

# Making symbolic links to the variant histones  data obtained in collaboration with the researcher Roson et al.
ln -s /project/carol/dnascent/project/atlas/metaData/Epi1-h2bv.bw input/
ln -s /project/carol/dnascent/project/atlas/metaData/Epi1-h3.bw input/
ln -s /project/carol/dnascent/project/atlas/metaData/Epi2-h2bv.bw input/
ln -s /project/carol/dnascent/project/atlas/metaData/Epi2-h3.bw input/
ln -s /project/carol/dnascent/project/atlas/metaData/Epi-h4v_1.bw input/
ln -s /project/carol/dnascent/project/atlas/metaData/Epi-h4v_2.bw input/
ln -s /project/carol/dnascent/project/atlas/metaData/Epi-h4v_3.bw input/
ln -s /project/carol/dnascent/project/atlas/metaData/tryCru-clb1.fasta input/

# Writing the command lines necessary to run this pipeline.
cat > cmdLine.bash << EOI
# Copying the template scripts.
cp -a /project/carol/dnascent/pipeline/atlas/05a-heatmapOrc1Cdc6/deeptools/* .

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
