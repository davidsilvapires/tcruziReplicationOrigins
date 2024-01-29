###################################################################################################
#                             origins-originsGenomeCompartmentDistribution.R                      #
#                                                                                                 #
#    In this script, we will graphically present the annotation of DNA replication origins in     #
# the genomic Corpmartment of Trypanosoma cruzi. The data were generated in previous processes    #
# using bedtools' intersect tools, as described in the earlier stages of the pipeline.            #
# We suggest executing this script in Rstudio for enhanced visualization, facilitating specific   #
# adjustments.                                                                                    #
# Usage: saveCommand origins-originsGenomeCompartmentDistribution.R                               #
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
setwd("~/Documents/myScripts/RScript/BarPerGraph")

# Import the processed data on the server to plot the graph.
# Control data
totalbp                     <- import("controlCompartment.txt")
colnames(totalbp)           <- c('Compartment', 'Peaks')
totalbp

# Experiment data predominant
predominant                 <- import("Predominant-Compartment.txt")
colnames(predominant)       <- c('Compartment', 'Peaks')
predominant

# Experiment data flexible
flexible                   <- import("Flexible-Compartment.txt")
colnames(flexible)         <- c('Compartment', 'Peaks')
flexible

# Experiment data dormant
dormant                     <- import("Dormant-Compartment.txt")
colnames(dormant)           <- c('Compartment', 'Peaks')
dormant


# Defining the factors for the data frame (database)
C_factor <- c(rep('Total bp', dim(totalbp)[1]))
p_factor <- c(rep('Predominant', dim(predominant)[1]))
f_factor <- c(rep('Flexible', dim(flexible)[1]))
d_factor <- c(rep('Dormant', dim(dormant)[1]))

# Column composition to factors
totalbp$V3        <- C_factor
predominant$V3    <- p_factor
flexible$V3       <- f_factor
dormant$V3        <- d_factor

# Creating dataframe
df               <- rbind(totalbp,predominant, flexible, dormant)
colnames(df)     <-c('Compartment','Peaks', 'Groups')
df$Compartment=factor(df$Compartment, levels = c("Core", "Disruptive", "Both" ) )
df$Groups=factor(df$Groups, levels = c("Dormant","Flexible","Predominant","Total bp"))
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
ggsave("Atlas/originsGenomeCompartment.pdf",f, height =4 , width=6 , dpi= 1200)
