###################################################################################################
#                            orc1Cdc6-genomeFeatureDistribution.R                                 #
#                                                                                                 #
# In this script, we will graphically present the distances between the peaks of Orc1cdc6 ChIP-seq#
# and DNA replication origins. The data were generated in previous processes using bedtools       #
# tools, as described in the preceding steps of the pipeline.                                     #
#                                                                                                 #
# We suggest executing this script in Rstudio for enhanced visualization, facilitating specific   #
# adjustments.                                                                                    #
#                                                                                                 #
# Usage: saveCommand orc1Cdc6-genomeFeatureDistribution.R                                         #
#                                                                                                 #
# Copyleft (ɔ) 2023 by Thiago Andrade Franco                                                      #
#                    <thiago.franco.esib@esib.butantan.gov.br                                     #
#    01/09/2023: First version.                                                                   #
###################################################################################################

# First step download the necessary libraries
packs <- c("ggplot2", "rio", "dplyr","plyr", "ggsignif", 
           "rstatix", "RVAideMemoire", "car") 

# Checking the packages installed
if(sum(as.numeric(!packs %in% installed.packages())) != 0){
  installer <- packs[!packs %in% installed.packages()]
  for(i in 1:length(installer)) {
    install.packages(installer, dependencies = T)
    break()}
  sapply(packs, require, character = T) 
} else {
  sapply(packs, require, character = T) 
}

# Adjust the working directory
setwd("~/Documents/myScripts/RScript/BoxplotGraphic")

# Insert the data for the data frame composition
mfa              <-import(file="closetmfa.txt")
dnascent         <-import(file="closetdnascent.txt")

# Factor composition
mfa$Origins                 <- rep('MFA-Seq',dim(mfa)[1])
colnames(mfa)               <- c('Length', 'Origins')
dnascent$Origins            <- rep('D-NAscent',dim(dnascent)[1])
colnames(dnascent)          <- c('Length', 'Origins')

# In this step, we separate two different approaches. 
# First, we look at the length of genes (or replicated origins).
# In this case, we use the Mann-Whithey test (non-parametric) following the Bonferroni post hoc test. 
# And then, we look at the frequency of genes. 

# dataframe composition (database)
df2<- rbind(mfa[,1:2], dnascent[,1:2])
df2
df2$Origins <- factor(df2$Origins, levels=c("MFA-Seq","D-NAscent"))
df2

### Confirm your matrix is a data frame
is.data.frame(df2)
df2

#view data frame
glimpse(df2)

# Next step, statistical analysis
mu <- ddply(df2, "Origins", summarise, grp.mean=mean(Length), grp.median=median(Length), grp.moda= mode(Length))
head(mu)

# Check data distribution
byf.shapiro(Length ~ Origins, df2)
shapiro.test(df2$Length~Origins)

# Verification of the homogeneity of variances
leveneTest(Length~Origins, df2, center=median)

# Mann-whitney  (non-parametric) = to length 
wilcox.test(Length~Origins,data=df2)

# post-hoc test To length
dunn_test(Length ~ Origins, data = df2, p.adjust.method = "bonferroni")

# Descriptive data analysis
df2 %>% group_by(Origins) %>% 
  get_summary_stats(Length, type = "median_iqr")

# Data visualization
Box <- ggplot(df2, aes(x=Origins, y =Length)) +
  geom_boxplot(aes(color = Origins, fill = Origins),
               alpha = 0.3, position = "identity") +
  theme_classic() +
  labs(y="Distance (bp)", x=" Origins") +
  theme(axis.text = element_text(size = 10 , color ="black"), 
        axis.title = element_text(size=12,color = "black")) +
  geom_signif(comparisons = list(c("MFA-Seq","D-NAscent"))),
              map_signif_level = TRUE,y_position = c(200000 ), textsize=2, test = "wilcox.test") +
  coord_cartesian(ylim = c(0, 250000))
Box

# Export graph
ggsave("Atlas/DistancepeakOrc1Cdc6ToMFAandDNAscentOrigins.pdf",Box, height = 4, width= 6, dpi= 3000)
