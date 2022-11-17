## ---------------------------
##
## metagene.cx.track_maker.r
##
## Builds a methylation bed file for downprocessing
## Author: Juan Santos, SLU, 2021 December
##
## ---------------------------

## Cleans the workspace and disable scientific notation 
rm(list = ls())
options(scipen=999)

## sets working directory
setwd("/Users/jsantos/metagene_express/")

## reads the input
read.table("sample.cx_report.CG.txt", header = FALSE) -> data1

data1[data1$V1!= "chloroplast" & data1$V1!= "mitochondria", ] -> data1
data.frame(data1[, 1:6]) -> data1
colnames(data1) <- c("chr", "pos", "strand", "met", "unmet", "con")

data1["name"] <- "."
data1["start"] <- data1$pos - 1
data1["methyl"] <- data1$met/(data1$met + data1$unmet)

head(data1)
str(data1)

###  selects the columns with the values in the right order
track <- data.frame(data1$chr, data1$start, data1$pos, data1$name, data1$methyl, data1$strand)
track <- data1[ , c("chr", "start", "pos", "methyl") ]
head(track)

track <- na.omit(track)

## writes a first bed table
write.table(track, "track.CG.bed", quote =  FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")

system("sort -k1,1 -k2,2n track.CG.bed> track.sorted.CG.bed; rm track.CG.bed", intern=TRUE)
