###################################################################################################
#                             orc1Cdc6-genomeFeatureDistribution.R                                #
#                                                                                                 #
#  In this script, we will graphically present the annotation of peaks from Orc1cdc6 ChIP-seq on  #
# the genomic features of Trypanosoma cruzi. The data were generated in previous processes using  #
# bedtools' intersect tools, as described in the earlier stages of the pipeline.                  #
# We suggest executing this script in Rstudio for enhanced visualization, facilitating specific   #
# adjustments.                                                                                    #
# Usage: saveCommand orc1Cdc6-genomeFeatureDistribution.R                                         #
#                                                                                                 #
# Copyleft (ɔ) 2023 by Thiago Andrade Franco                                                      #
#                    <thiago.franco.esib@esib.butantan.gov.br                                     #
#    01/09/2023: First version.                                                                   #
###################################################################################################

# First step, load the required libraries; if not already installed, install and then load.
packs <- c("ggplot2", "dplyr","rio","psych","plyr")
if(sum(as.numeric(!packs %in% installed.packages())) != 0){
  installer <- packs[!packs %in% installed.packages()]
  for(i in 1:length(installer)) {
    install.packages(installer, dependencies = T)
    break()}
  sapply(packs, require, character = T)
} else {
  sapply(packs, require, character = T)
}

# Second step, adjust the working directory
setwd("~/Documents/myScripts/RScript/BarPerGraph")

# Import the processed data on the server to plot the graph.
# Control data
totalbp                     <- import("controlFeature.txt")
colnames(totalbp)           <- c('Feature', 'Peaks')
totalbp

# Orc1Cdc6 peaks
orc1cdc6                    <- import("orc1cdc6Feature.txt")
colnames(orc1cdc6)          <- c('Feature', 'Peaks')
orc1cdc6

# Defining the factors for the data frame (database)
C_factor                    <- c(rep('Total bp', dim(totalbp)[1]))
o_factor                    <- c(rep('Orc1Cdc6', dim(orc1cdc6)[1]))

# Column composition to factors
totalbp$V3                   <- C_factor
orc1cdc6$V3                  <-o_factor

# Creating dataframe
df                           <- rbind(totalbp, orc1cdc6)
colnames(df)                 <-c('Feature','Peaks', 'Groups')
df

# Confirm your matrix is a data frame
is.data.frame(df)
df
df$Feature=factor(df$Feature, levels = c("CDS","cSSR", "dSSR", 
                                         "Inter-PTU", "snoRNA",
                                         "rRNA", "tRNA"))
df

# ploting bar grapgh  with percentage
feature <- ggplot(df, aes(x = Groups,y = Peaks, fill = Feature)) +
  geom_bar(stat="identity", position ="fill") +
  theme_minimal() +
  theme(axis.text = element_text(size = 8 , color ="black"),
        legend.position = "bottom") +
  labs(x =NULL, y = NULL) +
  scale_y_continuous(labels = scales::percent) +
  coord_flip() +
  scale_fill_manual(values=c( "#f0027f","#ffff99","#386cb0","#fdc086","#beaed4","#bf5b17", "#7fc97f"))

feature

# Exporting the graphic in pdf
ggsave("Atlas/orc1cdc6GenomeFeature.pdf",feature, height =2 , width= 4, dpi= 1200)

