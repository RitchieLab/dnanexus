#!/usr/bin/env Rscript

library(MASS)

args<-commandArgs(TRUE)

# get the arguments\
ancfn <- args[1]
evfn <- args[2]
thresh <- as.numeric(args[3])
numev <- as.integer(args[4])

anc <- read.table(ancfn, header=F, sep='\t')
evec <- read.table(evfn, header=F)
anc_known <- anc[(anc$V2 != "Unknown" & anc$V2 != "Other"), ]

for(i in 2:ncol(evec)){ colnames(evec)[i] <- paste("PC", as.character(i-1), sep='')}
colnames(evec)[1] = "ID"
colnames(anc_known) <- c("ID", "race")
evec_race <- merge(evec, anc_known)
evec_race$race <- factor(evec_race$race)
qres <- qda(evec_race[,2:(numev + 1)], evec_race$race)
p_anc <- predict(qres, evec[,2:(numev + 1)])

# now, get the ancestry + max probability
race_df <- data.frame(ID=evec$ID, race=apply(pqd$posterior, 1, which.max), prob=apply(pqd$posterior, 1, max))
race_labels <- levels(evec_race$race)
race_df$race <- race_labels[race_df$race]

# and compare the max posterior probability with the threshold and set to "NA/Other"
race_df$race[race_df$prob < thresh] <- "Other"

write.table(race_df, row.names=F, sep='\t', quote=F)
