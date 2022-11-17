## ---------------------------
##
## metagene.matrixmaker.r
##
## Makes a data matrix genes/TES x bins  with coverage values
## Author: Juan Santos, SLU, 2021 December
##
## ---------------------------

## Cleans the workspace and disable scientific notation 
rm(list = ls())
options(scipen=999)

## Sets working directory
setwd("/Users/jsantos/Documents/github_ingefast/metagene_express/example")

## Defines the input file  - a bed file with normalised coverage/methylation values
Sys.setenv(sample_file = "/Users/jsantos/Documents/github_ingefast/metagene_express/example/track.sorted.bed")

## Defines the dir with the genomic reference files
Sys.setenv(genomic_reference_files="/Users/jsantos/Documents/github_ingefast/metagene_express/example/genomic_reference_files")


##################################################################################
#################### maps coverage files into defined bins   #####################
##################################################################################

system("bedtools map -a $genomic_reference_files/tss_plus.sorted.bed -b $sample_file -c 4 -o mean >track.tss_plus.sorted.txt;", intern=TRUE)

system("bedtools map -a $genomic_reference_files//tss_minus.sorted.bed -b $sample_file -c 4 -o mean >track.tss_minus.sorted.txt;", intern=TRUE)

system("bedtools map -a $genomic_reference_files/body_plus.sorted.bed -b $sample_file -c 4 -o mean >track.body_plus.sorted.txt;", intern=TRUE)

system("bedtools map -a $genomic_reference_files/body_minus.sorted.bed -b $sample_file -c 4 -o mean >track.body_minus.sorted.txt;", intern=TRUE)

system("bedtools map -a $genomic_reference_files/tes_plus.sorted.bed -b $sample_file -c 4 -o mean >track.tes_plus.sorted.txt;", intern=TRUE)

system("bedtools map -a $genomic_reference_files/tes_minus.sorted.bed -b $sample_file -c 4 -o mean >track.tes_minus.sorted.txt;", intern=TRUE)


######################################################################################
##  This works out the final data frame for the PLUS strand
######################################################################################

## This reshapes the TSS plus

tss <-read.table("track.tss_plus.sorted.txt", na.strings='.', header=FALSE)
colnames(tss) <- c("chr", "start", "end", "id", "win", "strand", "count")

tss_matrix <- reshape(tss[c(4, 5, 7)], idvar ="id", timevar= "win", direction = "wide")

## This reshapes the BODY plus

body <- read.table("track.body_plus.sorted.txt", na.strings='.', header=FALSE)
colnames(body) <- c("chr", "start", "end", "id", "win", "strand", "count")

body_matrix <- reshape(body[c(4, 5, 7)], idvar ="id", timevar= "win", direction = "wide")

## This reshapes the TES plus

tes <- read.table("track.tes_plus.sorted.txt", na.strings='.', header=FALSE)
colnames(tes) <- c("chr", "start", "end", "id", "win", "strand", "count")

tes_matrix <- reshape(tes[c(4, 5, 7)], idvar ="id", timevar= "win", direction = "wide")

## This puts the TSS, BODY and TES together

profile_plus <- cbind(tss_matrix, body_matrix[ ,2:41], tes_matrix[ ,2:21])

colnames(profile_plus) <- c("id", paste("tss", seq(from=2000, to=100, by=-100), sep="_"),  paste("body", seq(from=1, to=40, by=1), sep="_"), paste("tes", seq(from=100, to=2000, by=100), sep="_"))


######################################################################################
##  This works out the final data frame for the MINUS strand
######################################################################################

##  This reshapes the TSS minus

tss <- read.table("track.tss_minus.sorted.txt", na.strings='.', header=FALSE)
colnames(tss) <- c("chr", "start", "end", "id", "win", "strand", "count")

tss_matrix <- reshape(tss[c(4, 5, 7)], idvar ="id", timevar= "win", direction = "wide")

## This reshapes the BODY minus

body <- read.table("track.body_minus.sorted.txt", na.strings='.', header=FALSE)
colnames(body) <- c("chr", "start", "end", "id", "win", "strand", "count")

body_matrix <- reshape(body[c(4, 5, 7)], idvar ="id", timevar= "win", direction = "wide")

## This reshapes the TES minus

tes <- read.table("track.tes_minus.sorted.txt", na.strings='.', header=FALSE)
colnames(tes) <- c("chr", "start", "end", "id", "win", "strand", "count")

tes_matrix <- reshape(tes[c(4, 5, 7)], idvar ="id", timevar= "win", direction = "wide")

## This puts the TSS, BODY and TES together

profile_minus <- cbind(tss_matrix[ ,1], rev(tss_matrix[ , 2:21]), rev(body_matrix[ , 2:41]),  rev(tes_matrix[ , 2:21]))

colnames(profile_minus) <- c("id", paste("tss", seq(from=2000, to=100, by=-100), sep="_"),  paste("body", seq(from=1, to=40, by=1), sep="_"), paste("tes", seq(from=100, to=2000, by=100), sep="_"))


######################################################################################
## This binds MINUS and PLUS together and saves the FINAL data frame and saves output
######################################################################################

profile_total <- rbind(profile_plus, profile_minus)
profile_total <- profile_total[order(profile_total$id), ]

write.table(profile_total, "track.matrix.txt", quote =  FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")

system("rm track.tes* track.tss* track.input*", intern=TRUE)

q("no");
