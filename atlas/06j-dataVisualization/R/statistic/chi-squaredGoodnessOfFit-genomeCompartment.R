###################################################################################################
#                            chi-squaredGoodnessOfFit.R                                           #
#                                                                                                 #
#  In this script, we will present the statistical analysis employed for annotating peaks of      #
# Orc1Cdc6 and replication origins within the genomic Compartment of Trypanosoma cruzi. As the    #
# data exhibit a non-parametric nature and follow a standard distribution of features, we         #
# utilized the chi-square goodness-of-fit test for analysis.                                      #
#                                                                                                 #
# We suggest executing this script in Rstudio for enhanced visualization, facilitating specific   #
# adjustments.                                                                                    #
#                                                                                                 #
# Usage: saveCommand chi-squaredGoodnessOfFit.R                                                   #
#                                                                                                 #
# Copyleft (ɔ) 2023 by Thiago Andrade Franco                                                      #
#                    <thiago.franco.esib@esib.butantan.gov.br                                     #
#    01/09/2023: First version.                                                                   #
###################################################################################################
# This experiment will now be statistically analyzed.
# Initially, a contingency table containing the data to be analyzed and the
# reference data (Control) was created.
# First step download the necessary libraries
packs <- c("ggplot2", "dplyr","rstatix", "psych", "rio", "plyr","ggrepel", "lsr")
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
setwd("~/Documents/myScripts/RScript/Statistic")

# Now, we will import the contingency table to Chi-Square
table <- read.csv2('contigencyTableCDS.csv',row.names = 1)
table

# Afterwards, we will perform the chi-square test
# Given that our data had frequencies lower than 5, we used the parameter
# "simulate.p.value = TRUE,B = 1000" to simulate a frequency with 1000 replications.
chisq <- chisq.test(table, p= c(0.83,0.131,0.039), simulate.p.value = TRUE,B = 1000)
chisq

chisq <- chisq.test(table, p= c(0.83,0.131,0.039))
chisq

# The following command extends the number of lines of printing your results in the console
options(scipen=999)

#  Analysis of expected frequencies
chisq$expected
chisq$estimate
chisq$observed

# Analysis of adjusted standardized residuals
# Standardized Residual (SPSS) - Pearson Residuals:
chisq$residuals

# Adjusted standardized residual (SPSS):
chisq$stdres

# Adjusted standardized residual > 1.96 or < -1.96 -- 5% alpha
# Calculation of the cut-off point for the standardized residues
# Calculate the new alpha:
# Where "l" is the number of rows and "c" is the number of columns
# We will divide the4 0.05 by the product c*l (number of cells)
newalpha <- 0.05/length(table)
newalpha

# Calculate the cut point, based on the new alpha:
# The division by two is because it is two-tailed
qnorm(newalpha/2)  # this number will dictate the range
# Significant residuals: > 2.39 or < -2.39 -- new alpha: 0.007142857
# Calculation of p for residuals
round(2*(1-pnorm(abs(chisq$stdres))),20)

# Effect Size - Cramer's V
cramersV(table,  p= c(0.83,0.131,0.039))

# The interpretation depends on the degrees of freedom (df):
# gl = columns-1
# In this case: df =  and Cramer's V corresponds to a small effect size (Cohen, 1988)
gl <- (ncol(table)-1)
gl

#  p value View
chisq$p.value
