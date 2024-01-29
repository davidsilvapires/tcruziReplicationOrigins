###################################################################################################
#                                orc1Cdc6Free-genomeCompartment.R                                 #
#                                                                                                 #
#  In this script, we will graphically present the annotation of DNA replication origins without  #
# Orc1Cdc6 in the genomic compartment of Trypanosoma cruzi.The data were generated in previous    #
# processes using bedtools' subtract tools, as described in the earlier stages of the pipeline.   #
# We suggest executing this script in Rstudio for enhanced visualization, facilitating specific   #
# adjustments.                                                                                    #
#                                                                                                 #
# Usage: saveCommand orc1Cdc6Free-genomeCompartment.R                                             #
#                                                                                                 #
# Copyleft (ɔ) 2023 by Thiago Andrade Franco                                                      #
#                    <thiago.franco.esib@esib.butantan.gov.br                                     #
#    09/01/2024: First version.                                                                   #
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
setwd("~/Documents/myScripts/RScript/BarPerGraph")

# Import the processed data on the server to plot the graph.
# Control data
totalbp<- import("controlCompartment.txt")
colnames(totalbp)<- c('Compartment', 'Peaks')
totalbp

# Experiment data dormant
Orc1cdc6Free <- import("Orc1Cdc6-Free-Compartment.txt")
colnames(Orc1cdc6Free)<- c('Compartment', 'Peaks')
Orc1cdc6Free

# Defining the factors for the data frame (database)
C_factor <- c(rep('Total bp', dim(totalbp)[1]))
o_factor <- c(rep('Orc1cdc6-Free', dim(Orc1cdc6Free)[1]))

# Column composition to factors
totalbp$V3        <- C_factor
Orc1cdc6Free$V3       <- o_factor

# Creating dataframe
df <- rbind(totalbp,Orc1cdc6Free)
colnames(df)<-c('Compartment','Peaks', 'Groups')
df$Compartment=factor(df$Compartment, levels = c("Core", "Disruptive", "Both" ) )
df

# Confirm your matrix is a data frame
is.data.frame(df)
df

# ploting bar grapgh  with percentage
f <- ggplot(df, aes(x = Groups,y = Peaks, fill = Compartment)) +
  geom_bar(stat="identity", position ="fill") +
  theme_minimal() +
  theme(axis.text = element_text(size = 8 , color ="black")) +
  labs(x =NULL, y = NULL) +
  scale_y_continuous(labels = scales::percent) +
  coord_flip() +
  scale_fill_manual(values=c("#1b9e77","#d95f02","#7570b3"))
f

# Exporting the graphic in pdf
ggsave("Atlas/Orc1Cdc6Free-GenomeCompartment.pdf",f, height =5 , width=7.5 , dpi= 1200)
