## ---------------------------
##
## metagene.heatmap.r
##
## Plots heatmap of single genomic coverage and DNA methylation data 
## Author: Juan Santos, SLU, 2021 December
##
## ---------------------------

## Cleans the workspace and disable scientific notation 
rm(list = ls())
options(scipen=999)

## Sets working directory
setwd("/Users/jsantos/metagene_express/")

## loads libraries
library(gplots)
library(pheatmap)  

## reads the data and split after genomic feature

track1 <- c("track.matrix.txt")

TRACK_matrix1 <- read.table(track1, check.names = FALSE, header=TRUE)
TRACK_matrix_gene1 <- TRACK_matrix1[ grep("G", TRACK_matrix1$id), ]
TRACK_matrix_te1 <- TRACK_matrix1[ grep("TE", TRACK_matrix1$id), ]


########################################################################################
### this plots the heatmap for genes
########################################################################################

##############  This can make some rescaling of max/min  ###############################

## this rescale the maximum values in color #gradient
TRACK_matrix_gene1[TRACK_matrix_gene1 >= 0.5] <- 1

## this rescale the minimum values in color gradient
#TRACK_matrix_gene1[TRACK_matrix_gene1 <= -1.5] <- -1.5 

##############  This bring a certain ordering of rows #################################

## Ordering criteria after row name ie genomic position
#TRACK_matrix_gene1 <- TRACK_matrix_gene1[sort(TRACK_matrix_gene1$id)), ]  

## Ordering criteria after decreasing/increasing  levels
TRACK_matrix_gene1 <- TRACK_matrix_gene1[order(apply(TRACK_matrix_gene1[, 2:81], 1, mean, na.rm=TRUE), decreasing = TRUE), ] 


## 
main_title_genes <- c("DNA methylation ratio (CG) in genes")

#png(paste("Heatmap.", main_title_genes,".png"), width = 3, height = 9, units = "in", res = 200, type = c("quartz"))

pheatmap( TRACK_matrix_gene1[ , 2:81], main = main_title_genes, fontsize = 8, cluster_rows = FALSE, cluster_cols = FALSE, show_rownames = F, show_colnames = F, annotation_names_col=F, border_color = NA, color = colorRampPalette(colors=c("green", "black", "red"))(100))

#dev.off()


########################################################################################
### this plots the heatmap for TEs
########################################################################################

##############  This can make some rescaling of max/min  ###############################

## this rescale the maximum values in color #gradient
TRACK_matrix_te1[TRACK_matrix_te1 >= 0.5] <- 1

## this rescale the minimum values in color gradient
#TRACK_matrix_te1[TRACK_matrix_te1 <= -1.5] <- -1.5 

##############  This bring a certain ordering of rows #################################

## Ordering criteria after row name ie genomic position
#TRACK_matrix_te1 <- TRACK_matrix_te1[sort(TRACK_matrix_te1$id)), ]  

## Ordering criteria after decreasing/increasing  levels
TRACK_matrix_te1 <- TRACK_matrix_te1[order(apply(TRACK_matrix_te1[, 2:81], 1, mean, na.rm=TRUE), decreasing = TRUE), ] 


## 
main_title_tes <- c("DNA methylation ratio (CG) in TEs")

#png(paste("Heatmap.", main_title_tes,".png"), width = 3, height = 9, units = "in", res = 200, type = c("quartz"))

pheatmap( TRACK_matrix_te1[ , 2:81], main = main_title_tes, fontsize = 8, cluster_rows = FALSE, cluster_cols = FALSE, show_rownames = F, show_colnames = F, annotation_names_col=F, border_color = NA, color = colorRampPalette(colors=c("green", "black", "red"))(100))

#dev.off()