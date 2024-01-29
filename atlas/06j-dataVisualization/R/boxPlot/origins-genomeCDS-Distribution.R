###################################################################################################
#                             origins-genomeFeatureDistribution.R                                #
#                                                                                                 #
#    In this script, we will graphically present the annotation of DNA replication origins in     #
# the genomic CDS of Trypanosoma cruzi. The data were generated in previous processes using       #
# bedtools' intersect tools, as described in the earlier stages of the pipeline.                  #
# We suggest executing this script in Rstudio for enhanced visualization, facilitating specific   #
# adjustments.                                                                                    #
# Usage: saveCommand origins-genomeFeatureDistribution.R                                         #
#                                                                                                 #
# Copyleft (ɔ) 2023 by Thiago Andrade Franco                                                      #
#                    <thiago.franco.esib@esib.butantan.gov.br                                     #
#    01/09/2023: First version.                                                                   #
###################################################################################################

# This script allows you to make a bar plot with the percentage
# First step download the necessary libraries
packs <- c("ggplot2", "dplyr","rio","scales")
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

# Import data from the contingency table (Control)
totalbp               <- import("controlGene.txt")
colnames(totalbp)     <- c('CDS', 'Peaks')
totalbp

# Import data from the contingency table (Predominant)
Predominant              <- import("predominantGene.txt")
colnames(Predominant )   <- c('CDS', 'Peaks')
Predominant

# Import data from the contingency table (Flexible)
Flexible              <- import("flexibleGene.txt")
colnames(Flexible)    <- c('CDS', 'Peaks')
Flexible

# Import data from the contingency table (Dormant)
Dormant                <- import("dormantGene.txt")
colnames(Dormant)      <- c('CDS', 'Peaks')
Dormant

# Defining the factors for the data frame (database)
C_factor <- c(rep('Total bp', dim(totalbp)[1]))
P_factor <- c(rep('Predominant', dim(Predominant)[1]))
F_factor <- c(rep('Flexible', dim(Flexible)[1]))
D_factor <- c(rep('Dormant',dim(Dormant)[1]))

# Column composition Control
totalbp$V3      <- C_factor
Predominant$V3  <- P_factor 
Flexible$V3     <- F_factor
Dormant$V3      <- D_factor

# Creating dataframe
df               <- rbind(totalbp, Predominant,Flexible,Dormant)
colnames(df)     <-c('CDS','Peaks', 'Groups')
df$Groups=factor(df$Groups, levels = c("Dormant","Flexible","Predominant","Total bp"))
df$CDS=factor(df$CDS, levels = c( "ATP-dependent", "C/D","Cysteine","DGF-1",
                                  "Elongation", "Flagellar","GP63","Hypothetical",
                                  "Histone", "MASP", "Mucin", "Receptor-type", 
                                  "Retrotransposon", "Trans-sialidase",
                                  "UDP-Gal", "Others"))
df

# Confirm your matrix is a data frame
is.data.frame(df)
df

# ploting bar grapgh  with percentage
CDS <- ggplot(df, aes(x = Groups,y = Peaks, fill = CDS)) +
  geom_bar(stat="identity", position ="fill") +
  theme_minimal() +
  theme(axis.text = element_text(size = 8 , color ="black")) +
  labs(x =NULL, y = NULL) +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values=c( "#a6611a","#9e0142","#d53e4f","#f46d43","#fdae61","#fee08b","#ffffbf","#e6f598","#abdda4","#66c2a5","#3288bd","#5e4fa2", "#d7191c","#fdae61","#abdda4","#2b83ba")) +
  coord_flip()
CDS

# Exporting the graphic in pdf
ggsave("Atlas/OriginsGenomeCDS.pdf",CDS, height =4 , width= 6, dpi= 1200)

