###################################################################################################
#                                orc1Cdc6-genomeCDSDistribution.R                                 #
#                                                                                                 #
#  In this script, we will graphically present the annotation of peaks from Orc1cdc6 ChIP-seq on  #
# the genomic CDS of Trypanosoma cruzi. The data were generated in previous processes             #
# using  bedtools' intersect tools, as described in the earlier stages of the pipeline.           #
# We suggest executing this script in Rstudio for enhanced visualization, facilitating specific   #
# adjustments.                                                                                    #
# Usage: saveCommand orc1Cdc6-genomeCDSDistribution.R                                             #
#                                                                                                 #
# Copyleft (ɔ) 2023 by Thiago Andrade Franco                                                      #
#                    <thiago.franco.esib@esib.butantan.gov.br                                     #
#    09/01/2024: First version.                                                                   #
###################################################################################################
# First step, load the required libraries; if not already installed, install and then load.
packs <- c("ggplot2", "dplyr","rio")
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
totalbp                                  <- import("controlGene.txt")
colnames(totalbp)                        <- c('CDS', 'Peaks')
totalbp

# Import data from the contingency table (Predominant)
Orc1cdc6                                 <- import("orc1cdc6Gene.txt")
colnames(Orc1cdc6)                       <- c('CDS', 'Peaks')
Orc1cdc6

# Defining the factors for the data frame (database)
C_factor                                 <- c(rep('Total bp', dim(totalbp)[1]))
O_factor                                 <- c(rep('Orc1cdc6', dim(Orc1cdc6)[1]))

# Column composition Control
totalbp$V3                               <- C_factor
Orc1cdc6$V3                              <- O_factor

# Creating dataframe
df                                        <- rbind(totalbp,Orc1cdc6)
colnames(df)                              <-c('CDS','Peaks', 'Groups')
df

### Confirm your matrix is a data frame
is.data.frame(df)
df
df$CDS=factor(df$CDS, levels = c( "ATP-dependent", "C/D","Cysteine","DGF-1",
                                  "Elongation", "Flagellar","GP63","Hypothetical",
                                  "Histone", "MASP", "Mucin", "Receptor-type", 
                                  "Retrotransposon", "Trans-sialidase",
                                  "UDP-Gal", "Others"))
df
# ploting bar grapgh  with percentage
### For this graph we don't need creating percentade column just use scaly_y_continuos(labels=scales::percent)
CDS <- ggplot(df, aes(x = Groups,y = Peaks, fill = CDS)) +
  geom_bar(stat="identity", position ="fill") +
  theme_minimal() +
  theme(axis.text = element_text(size = 8 , color ="black"),
        legend.position = "bottom") +
  labs(x =NULL, y = NULL) +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values=c( "#a6611a","#9e0142","#d53e4f","#f46d43","#fdae61","#fee08b","#ffffbf","#e6f598","#abdda4","#66c2a5","#3288bd","#5e4fa2", "#d7191c","#fdae61","#abdda4","#2b83ba")) +
  coord_flip()


CDS
# Exporting the graphic in pdf
ggsave("Atlas/orc1cdc6GenomeCDS.pdf",CDS, height =3 , width= 4.5, dpi= 1200)
