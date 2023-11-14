#!/usr/bin/env bash

###################################################################################################
#                                  deeptoolsPre                                                   #
#                                                                                                 #
# This script is a required pre-processing step to plotHeatmap to Orc1cdc6.                       #
# It sets up the directory structure and creates symbolic links for input data.                   # 
#                                                                                                 #
# Usage: saveCommand script/deeptoolsPre.bash                                                      #
#                                                                                                 #
# Copyleft (ɔ) 2023 by Thiago Andrade Franco                                                      #
#                      <thiago.franco.esib@esib.butantan.gov.br>                                  #
#    12/07/2023: First version.                                                                   #
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
final/predominantsOrigins-NotUniq.bed input/predominantOrigins.bed
ln -s /project/carol/dnascent/project/atlas/02-flexibleOrigins/bedtools/\
final/flexibleOrigins-NotUniq.bed input/flexibleOrigins.bed

# Making a symbolic link for input data from multigenic families
ln -s /project/carol/dnascent/project/atlas/00-genomeCompartmentComposition/release-32/\
final/TriTrypDB-32_TcruziCLBrenerEsmeraldo-like-gene-DGF.bed input/DGF-1-Gene.bed
ln -s /project/carol/dnascent/project/atlas/00-genomeCompartmentComposition/release-32/\
final/TriTrypDB-32_TcruziCLBrenerEsmeraldo-like-gene-MASP.bed input/MASP-Gene.bed
ln -s /project/carol/dnascent/project/atlas/00-genomeCompartmentComposition/release-32/\
final/TriTrypDB-32_TcruziCLBrenerEsmeraldo-like-gene-RHS.bed input/RHS-Gene.bed
ln -s /project/carol/dnascent/project/atlas/00-genomeCompartmentComposition/\
release-32/final/TriTrypDB-32_TcruziCLBrenerEsmeraldo-like-gene-TS.bed input/TS-Gene.bed
ln -s /project/carol/dnascent/project/atlas/00-genomeCompartmentComposition/\
release-32/final/TriTrypDB-32_TcruziCLBrenerEsmeraldo-like-gene-mucin.bed input/mucin-Gene.bed

# Making a symbolic link for input data from genome feature
ln -s /project/carol/chiporcs/rawData/labData/bedFeatures/output_bedFiles\
/Intergenic_Polycistron_Region.bed input/intergenicPolycistronRegions-Feature.bed
ln -s /project/carol/chiporcs/rawData/labData/bedFeatures/output_bedFiles\
/snoRNA.bed input/snoRNA-Feature.bed
ln -s /project/carol/dnascent/project/atlas/metaData/tryCru-clb1.fasta\
  input/tryCru-clb1.fa

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
