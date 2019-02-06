# I love snail mucus

# Transpose data from excel sheet to data.frame, add columns for snail and coral

data <- t(Updated_Abundance_02_05_19)
df <- as.data.frame(data[,1:n])
colnames(df) <- c(df[1,])
df <- df[-1,]

coral_status <- t(corralivore_table)
coral_status <- as.data.frame(coral_status[,1:n])
colnames(coral_status) <- c(coral_status[1,])
coral_status <- coral_status[-1,]

n <- dim(data)[2] - 3

# Now let's make our full combined dataset.
fulldf <- cbind(df, coral_status)

# Sanity checking that things are there and that we can filter correctly
unique(fulldf$Feeding)
corals <- fulldf[fulldf$Feeding == "Corallivore",]
noncorals <- fulldf[fulldf$Feeding == "Non-corallivore",]
wilcox.test(corals$"1", noncorals$"1")

# Trying Hotellings T^2 just for kicks -- not part of the data analysis 
# install.packages("ICSNP")
# library("ICSNP")
# HotellingsT2(corals[,1:1252], noncorals[,1:1252])

# cols <- 1:1252
# 
# library(rrcov)
# T2.test(corals[,cols], noncorals[,cols])

# Testing how to run Mann-Whitney on 1252 columns
# column_number <- 4
# corals <- fulldf[fulldf$Feeding == "Corallivore", column_number]
# noncorals <- fulldf[fulldf$Feeding == "Non-corallivore", column_number]
# corals; noncorals
# names(wilcox.test(corals, noncorals))
# wilcox.test(corals, noncorals)$p.value

# function run_mann() takes a column of fulldf and runs two-sided Mann Whitney U test using function wilcox.test() on it, returns the p-value
run_mann <- function(column, data) {
  corals <- data[data$Feeding == "Corallivore", column]
  noncorals <- data[data$Feeding == "Non-corallivore", column]
  
  pval <- wilcox.test(corals, noncorals)$p.value
  return(pval)
}

# Tried using map() but kept getting errors -- instead resorted to using for loop
# library(tidyverse)
# cols <- 1:1252
# columns <- cols %>% map(run_mann(data = fulldf)) %>% unlist()
# str(columns)

pvals <- rep(NA, 1252)
columns <- 1:1252

# For loop iterates over all of the columns and runs run_mann() every time
for (i in columns) {
  pval <- suppressWarnings(run_mann(column = i, data = fulldf))
  pvals[i] <- pval
}

# Ties produce NaN (aka when we have nothing but 0s)
length(pvals[pvals == "NaN"])

# allp is a dataframe with the OTUs and corresponding Mann Whitney p-value
allp <- as.data.frame(cbind(columns, pvals)) %>% filter(pvals != 'NaN')
colnames(allp) <- c("taxa", "pvalues")
head(allp)

which(pvals <= .05)

# Bonferroni correction

p.sorted <- sort(pvals)
m <- length(pvals)

p.bonferroni <- pmin(m*p.sorted, 1)
p.bonferroni
which(p.bonferroni <= .05)

# Benjamini-Hochberg Procedure

p.hochberg <- p.adjust(p.sorted, method = "hochberg")
p.hochberg
which(p.hochberg <= .05)

# Q-values

install.packages("BiocManager")
BiocManager::install("qvalue", version = "3.8")
browseVignettes("qvalue")

library(qvalue)
qobj <- qvalue(p = pvals)
qobj
