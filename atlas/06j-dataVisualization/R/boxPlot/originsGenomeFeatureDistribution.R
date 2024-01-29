###################################################################################################
#                             origins-genomeFeatureDistribution.R                                 #
#                                                                                                 #
#    In this script, we will graphically present the annotation of DNA replication origins in     #
# the genomic features of Trypanosoma cruzi. The data were generated in previous processes        #
# using bedtools' intersect tools, as described in the earlier stages of the pipeline.            #
# We suggest executing this script in Rstudio for enhanced visualization, facilitating specific   #
# adjustments.                                                                                    #
# Usage: saveCommand origins-genomeFeatureDistribution.R                                          #
#                                                                                                 #
# Copyleft (ɔ) 2023 by Thiago Andrade Franco                                                      #
#                    <thiago.franco.esib@esib.butantan.gov.br                                     #
#    01/09/2023: First version.                                                                   #
###################################################################################################
# First step, load the required libraries; if not already installed, install and then load.
packs <- c("ggplot2", "dplyr","rio","scales","cowplot")
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
setwd("~/Documents/MyScripts/RScript/BarPerGraph")

# Import the processed data on the server to plot the graph.
# Control data
totalbp<- import("controlFeature.txt")
colnames(totalbp)<- c('Feature', 'Peaks')
totalbp

# Experiment data predominant
predominant <- import("predominantFeature.txt")
colnames(predominant)<- c('Feature', 'Peaks')
predominant

# Experiment data flexible
# Experiment data predominant
flexible <- import("flexibleFeature.txt")
colnames(flexible)<- c('Feature', 'Peaks')
flexible

dormant <- import("dormantFeature.txt")
colnames(dormant)<- c('Feature', 'Peaks')
dormant

# Defining the factors for the data frame (database)
C_factor <- c(rep('Total bp', dim(totalbp)[1]))
P_factor <- c(rep('Predominant', dim(predominant)[1]))
F_factor <-c(rep('Flexible', dim(flexible)[1]))
D_factor <-c(rep('Dormant', dim(dormant)[1]))

# Column composition to factors
totalbp$V3 <- C_factor
predominant$V3 <-P_factor
flexible$V3 <- F_factor
dormant$V3 <- D_factor

# Creating dataframe
df <- rbind(totalbp, predominant, flexible, dormant)
colnames(df)<-c('Feature','Peaks', 'Groups')
df

# Confirm your matrix is a data frame
is.data.frame(df)
df
df$Groups=factor(df$Groups, levels = c("Dormant","Flexible","Predominant","Total bp"))
df$Feature=factor(df$Feature, levels = c("CDS","cSSR", "dSSR", 
                                         "Inter-PTU", "snoRNA",
                                         "rRNA", "tRNA"))
df

# ploting bar grapgh  with percentage
feature <- ggplot(df, aes(x = Groups,y = Peaks, fill = Feature)) +
  geom_bar(stat="identity", position ="fill") +
  theme_minimal() +
  theme(axis.text = element_text(size = 12 , color ="black"),
        legend.position ="right") +
  labs(x =NULL, y = NULL) +
  scale_y_continuous(labels = scales::percent) +
  coord_flip() +
  scale_fill_manual(values=c( "#f0027f","#ffff99","#386cb0","#fdc086","#beaed4","#bf5b17",
                                       "#7fc97f"))

feature
# Exporting the graphic in pdf
ggsave("Atlas/originsGenomeFeature.pdf",feature, height =7 , width=10.5, dpi= 1200)
