## ---------------------------
##
## metagene.frame_builder.r
##
## Builds a set of frame file to define the binning in genes and TEs
## Author: Juan Santos, SLU, 2021 December
##
## ---------------------------

## Cleans the workspace and disable scientific notation 
rm(list = ls())
options(scipen=999)

## Sets working directory
setwd("/Users/jsantos/metagene_express/")

########################################################################################
## Inputs the annotation file in bed format
########################################################################################

master <- read.table("TAIR10_genes_transposons.bed", header=FALSE)

master_bed <- master[ , c(1:6)]
colnames(master_bed) <- c("chr", "start", "end", "id", "score", "strand")

master_bed_plus <- master_bed[master_bed$strand == "+", ]
master_bed_minus <- master_bed[master_bed$strand == "-", ]

write.table(master_bed_plus, "body_plus.bed", col.names=F, row.names=F, sep="\t", quote=F)
write.table(master_bed_minus, "body_minus.bed", col.names=F, row.names=F, sep="\t", quote=F)


## defines the 2Kb flanks on genomic features

system("bedtools flank -i body_plus.bed -g /Users/juan/Documents/genomes/Arabidopsis_thaliana/TAIR10/TAIR10.chrom.notPt.notMt.sizes -l 2000 -r 0  >tss_plus.bed;")

system("bedtools flank -i body_plus.bed -g /Users/juan/Documents/genomes/Arabidopsis_thaliana/TAIR10/TAIR10.chrom.notPt.notMt.sizes -l 0 -r 2000  >tes_plus.bed;")

system("bedtools flank -i body_minus.bed -g /Users/juan/Documents/genomes/Arabidopsis_thaliana/TAIR10/TAIR10.chrom.notPt.notMt.sizes -l 0 -r 2000 >tss_minus.bed;")

system("bedtools flank -i body_minus.bed -g /Users/juan/Documents/genomes/Arabidopsis_thaliana/TAIR10/TAIR10.chrom.notPt.notMt.sizes -l 2000 -r 0 >tes_minus.bed;")


## creates the binning

system("bedtools makewindows -b body_plus.bed -n 40 -i srcwinnum > body_plus.window.txt")

system("bedtools makewindows -b body_minus.bed -n 40 -i srcwinnum > body_minus.window.txt")

system("bedtools makewindows -b tss_plus.bed -w 100 -i srcwinnum > tss_plus.window.txt")

system("bedtools makewindows -b tes_plus.bed -w 100 -i srcwinnum > tes_plus.window.txt")

system("bedtools makewindows -b tss_minus.bed -w 100 -i srcwinnum > tss_minus.window.txt")

system("bedtools makewindows -b tes_minus.bed -w 100 -i srcwinnum > tes_minus.window.txt")


########################################################################################
## works out tss, tes and body per strand
########################################################################################

## tss_plus
tss_plus <-read.table("tss_plus.window.txt", header=F)
colnames(tss_plus) <- c("chr", "feature", "start", "idwin")

tss_plus$idwin <- as.character(tss_plus$idwin)
unlist(strsplit(tss_plus$idwin,"_")) -> pairs
matrix(pairs, ncol=2, byrow=TRUE) -> newpairs
tss_plus["id"] <- newpairs[,1]
tss_plus["win"] <- as.numeric(newpairs[,2])
tss_plus[order(tss_plus$chr, tss_plus$start),] -> tss_plus
tss_plus["strand"] <- "+"
tss_plus <- tss_plus[-4]
head(tss_plus)

write.table(tss_plus, "tss_plus.frame.bed", col.names=F, row.names=F, sep="\t", quote=F) 

## tss_minus
tss_minus <-read.table("tss_minus.window.txt", header=F)
colnames(tss_minus) <- c("chr", "feature", "start", "idwin")

tss_minus$idwin <- as.character(tss_minus$idwin)
unlist(strsplit(tss_minus$idwin,"_")) -> pairs
matrix(pairs, ncol=2, byrow=TRUE) -> newpairs
tss_minus["id"] <- newpairs[,1]
tss_minus["win"] <- as.numeric(newpairs[,2])
tss_minus[order(tss_minus$chr, tss_minus$start),] -> tss_minus
tss_minus["strand"] <- "-"
tss_minus <- tss_minus[-4]
head(tss_minus)

write.table(tss_minus, "tss_minus.frame.bed", col.names=F, row.names=F, sep="\t", quote=F) 

## body_plus
body_plus <-read.table("body_plus.window.txt", header=F)
colnames(body_plus) <- c("chr", "feature", "start", "idwin")

body_plus$idwin <- as.character(body_plus$idwin)
unlist(strsplit(body_plus$idwin,"_")) -> pairs
matrix(pairs, ncol=2, byrow=TRUE) -> newpairs
body_plus["id"] <- newpairs[,1]
body_plus["win"] <- as.numeric(newpairs[,2])
body_plus[order(body_plus$chr, body_plus$start),] -> body_plus
body_plus["strand"] <- "+"
body_plus <- body_plus[-4]
head(body_plus)

write.table(body_plus, "body_plus.frame.bed", col.names=F, row.names=F, sep="\t", quote=F) 

## body_minus
body_minus <-read.table("body_minus.window.txt", header=F)
colnames(body_minus) <- c("chr", "feature", "start", "idwin")

body_minus$idwin <- as.character(body_minus$idwin)
unlist(strsplit(body_minus$idwin,"_")) -> pairs
matrix(pairs, ncol=2, byrow=TRUE) -> newpairs
body_minus["id"] <- newpairs[,1]
body_minus["win"] <- as.numeric(newpairs[,2])
body_minus[order(body_minus$chr, body_minus$start),] -> body_minus
body_minus["strand"] <- "-"
body_minus <- body_minus[-4]
head(body_minus)

write.table(body_minus, "body_minus.frame.bed", col.names=F, row.names=F, sep="\t", quote=F) 

## tes_plus
tes_plus <-read.table("tes_plus.window.txt", header=F)
colnames(tes_plus) <- c("chr", "feature", "start", "idwin")

tes_plus$idwin <- as.character(tes_plus$idwin)
unlist(strsplit(tes_plus$idwin,"_")) -> pairs
matrix(pairs, ncol=2, byrow=TRUE) -> newpairs
tes_plus["id"] <- newpairs[,1]
tes_plus["win"] <- as.numeric(newpairs[,2])
tes_plus[order(tes_plus$chr, tes_plus$start),] -> tes_plus
tes_plus["strand"] <- "+"
tes_plus <- tes_plus[-4]
head(tes_plus)

write.table(tes_plus, "tes_plus.frame.bed", col.names=F, row.names=F, sep="\t", quote=F) 

## tes_minus
tes_minus <-read.table("tes_minus.window.txt", header=F)
colnames(tes_minus) <- c("chr", "feature", "start", "idwin")

tes_minus$idwin <- as.character(tes_minus$idwin)
unlist(strsplit(tes_minus$idwin,"_")) -> pairs
matrix(pairs, ncol=2, byrow=TRUE) -> newpairs
tes_minus["id"] <- newpairs[,1]
tes_minus["win"] <- as.numeric(newpairs[,2])
tes_minus[order(tes_minus$chr, tes_minus$start),] -> tes_minus
tes_minus["strand"] <- "-"
tes_minus <- tes_minus[-4]
head(tes_minus)

write.table(tes_minus, "tes_minus.frame.bed", col.names=F, row.names=F, sep="\t", quote=F) 


########################################################################################
## Sorts all the final bed frame files
########################################################################################

system("sort -k1,1 -k2,2n tss_plus.frame.bed > tss_plus.sorted.bed;")
system("sort -k1,1 -k2,2n tss_minus.frame.bed > tss_minus.sorted.bed;")

system("sort -k1,1 -k2,2n body_plus.frame.bed > body_plus.sorted.bed;")
system("sort -k1,1 -k2,2n body_minus.frame.bed > body_minus.sorted.bed;")

system("sort -k1,1 -k2,2n tes_plus.frame.bed > tes_plus.sorted.bed;")
system("sort -k1,1 -k2,2n tes_minus.frame.bed > tes_minus.sorted.bed;")

system("rm *frame.bed *s.bed *window.txt")
