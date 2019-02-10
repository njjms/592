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
fulldf[,dim(fulldf)[2] - 1]
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
length(pvals[pvals != "NaN"])
hist(pvals[pvals != "NaN"],
     main = "Histogram of Uncorrected p-values from 1094 Mann-Whitney U tests",
     xlab = "p-values")
remove_i <- which(pvals == "NaN")
remove_i
pvals <- pvals[-remove_i]
otus <- columns[-remove_i]

length(pvals); length(otus)

# allp is a dataframe with the OTUs and corresponding Mann Whitney p-value
allp <- as.data.frame(cbind(otus, pvals))
colnames(allp) <- c("taxa", "pvalues")
library(ggplot2)
ggplot(aes(x = pvalues), data = allp) + geom_histogram(fill = "dodgerblue1", bins = 50) +
  ggtitle("Histogram of Uncorrected p-values")
head(allp)
str(allp$pvalues)

which(allp$pvalues <= .05)
length(which(allp$pvalues <= .05))

# Bonferroni correction
length(order(allp$pvalues))
allp <- allp[order(allp$pvalues),]
m <- length(allp$pvalues); m

allp$p.bonferroni <- pmin(m*allp$pvalues, 1)
head(allp)
ggplot(aes(x = p.bonferroni), data = allp) + geom_histogram(fill = "dodgerblue1", bins = 50) +
  ggtitle("Histogram of Benjamini-Hochberg Corrected P-values")
which(allp$p.bonferroni <= .05)

# Benjamini-Hochberg Procedure

allp$p.hochberg <- p.adjust(allp$pvalues, method = "BH")
ggplot(aes(x = p.hochberg), data = allp) + geom_histogram(fill = "dodgerblue1", bins = 50) +
  ggtitle("Histogram of Benjamini-Hochberg Corrected P-values")
summary(allp)
which(allp$p.hochberg <= .05)
head(allp)

# Q-values
install.packages("BiocManager")
BiocManager::install("qvalue", version = "3.8")
browseVignettes("qvalue")

library(qvalue)
qobj <- qvalue(p = allp$pvalues)
allp$qvalues <- (qobj$qvalues)
head(allp, n = 50)
