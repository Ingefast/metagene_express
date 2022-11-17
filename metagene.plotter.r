## ---------------------------
##
## metagene.plotter.r
##
## Plots line graph and boxplots of genomic coverage and DNA methylation data 
## Author: Juan Santos, SLU, 2021 December
##
## ---------------------------

## Cleans the workspace and disable scientific notation 
rm(list = ls())
options(scipen=999)

## Sets working directory
setwd("/Users/jsantos/Documents/github_ingefast/metagene_express/example")


## Specifies input with genomic coverage matrix files

track1 <- c("sample1_track.matrix.txt")

track2 <- c("sample2_track.matrix.txt")

track3 <- c("sample3_track.matrix.txt")

track4 <- c("sample4_track.matrix.txt")


## reads the data and split after genomic feature

read.table(track1, check.names = FALSE, header=TRUE) -> TRACK_matrix1
TRACK_matrix_gene1 <- TRACK_matrix1[ grep("G", TRACK_matrix1$id), ]
TRACK_matrix_te1 <- TRACK_matrix1[ grep("TE", TRACK_matrix1$id), ]

read.table(track2, check.names = FALSE, header=TRUE) -> TRACK_matrix2
TRACK_matrix_gene2 <- TRACK_matrix2[ grep("G", TRACK_matrix2$id), ]
TRACK_matrix_te2 <- TRACK_matrix2[ grep("TE", TRACK_matrix2$id), ]

read.table(track3, check.names = FALSE, header=TRUE) -> TRACK_matrix3
TRACK_matrix_gene3 <- TRACK_matrix3[ grep("G", TRACK_matrix3$id), ]
TRACK_matrix_te3 <- TRACK_matrix3[ grep("TE", TRACK_matrix3$id), ]

read.table(track4, check.names = FALSE, header=TRUE) -> TRACK_matrix4
TRACK_matrix_gene4 <- TRACK_matrix4[ grep("G", TRACK_matrix4$id), ]
TRACK_matrix_te4 <- TRACK_matrix4[ grep("TE", TRACK_matrix4$id), ]


## Titles for the two main types of genomic features
main_title_gene <- c("Genes in CG context")
main_title_te <- c("TEs in CG context")


########################################################################################
#########################   Plotting of Genes ##########################################
########################################################################################

op <- par(mfrow = c(1, 2), mar = c(10, 7, 2, 2))


## Plots metagene plot as combined line graph

plot(colMeans(TRACK_matrix_gene1[, 2:81], na.rm = TRUE), col = "red4", lwd = 4, type = "l", ylim = c(0, 0.35), main = main_title_gene, xlab = "", ylab ="", xaxs = "i", yaxs = "i", cex.axis =2, axes = FALSE)
lines(colMeans(TRACK_matrix_gene2[, 2:81], na.rm = TRUE), col = "red1", lwd = 4, type = "l")
lines(colMeans(TRACK_matrix_gene3[, 2:81], na.rm = TRUE), col = "steelblue4", lwd = 4, type = "l")
lines(colMeans(TRACK_matrix_gene4[, 2:81], na.rm = TRUE), col = "steelblue1", lwd = 4, type = "l")

mtext("DNA methylation ratio (CG)", 2, line = 1, cex = 1)
axis(1, lwd = 3, at = c (1, 20, 60, 80), labels = c("-2kb", "TSS", "TES", "+2kb"), cex.axis = 1.3)
axis(2, at = c(0, 0.35), las = 2, lwd = 3, cex.axis = 1.3)
abline(v = c(20, 60), lty = 'longdash', lwd = 1.4)

legend("top", c( "wt_rep1", "wt_rep2", "mutantA_rep1", "mutantA_rep2" ), title.col = "black", col = c("red4", "red1", "steelblue4", "steelblue1"), inset=0.01, cex = 0.6, text.col = c("red4", "red1", "steelblue4", "steelblue1"), lty = c(1), lwd=c(4), pch = c(-1), bty="n")


## Plots a boxplot with means in as asterisks

boxplot(rowMeans(TRACK_matrix_gene1[, 21:61], na.rm = TRUE), rowMeans(TRACK_matrix_gene2[, 21:61], na.rm = TRUE), rowMeans(TRACK_matrix_gene3[, 21:61], na.rm = TRUE), rowMeans(TRACK_matrix_gene4[, 21:61], na.rm = TRUE), main = main_title_gene, ylab = "DNA methylation ratio (CG)", names = c("wt_rep1", "wt_rep2", "mutantA_rep1", "mutantA_rep2"),  col = c("red4", "red1", "steelblue4", "steelblue1"), outline = FALSE, cex.lab = 1, notch = T, las = 2)

mean.value <- c(  mean(rowMeans(TRACK_matrix_gene1[, 21:61], na.rm = TRUE), na.rm = TRUE), mean(rowMeans(TRACK_matrix_gene2[, 21:61], na.rm = TRUE), na.rm = TRUE), mean(rowMeans(TRACK_matrix_gene3[, 21:61], na.rm = TRUE), na.rm = TRUE), mean(rowMeans(TRACK_matrix_gene4[, 21:61], na.rm = TRUE), na.rm = TRUE)  )

points(seq(1:4), mean.value, pch = "*", col = "gold", cex = 2) 


########################################################################################
#########################   Plotting of TEs ##########################################
########################################################################################

op <- par(mfrow = c(1, 2), mar = c(10, 7, 2, 2))


## Plots metagene plot as combined line graph

plot(colMeans(TRACK_matrix_te1[, 2:81], na.rm = TRUE), col = "red4", lwd = 4, type = "l", ylim = c(0.2, 0.90), main = main_title_te, xlab = "", ylab ="", xaxs = "i", yaxs = "i", cex.axis =2, axes = FALSE)
lines(colMeans(TRACK_matrix_te2[, 2:81], na.rm = TRUE), col = "red1", lwd = 4, type = "l")
lines(colMeans(TRACK_matrix_te3[, 2:81], na.rm = TRUE), col = "steelblue4", lwd = 4, type = "l")
lines(colMeans(TRACK_matrix_te4[, 2:81], na.rm = TRUE), col = "steelblue1", lwd = 4, type = "l")

mtext("DNA methylation ratio (CG)", 2, line = 1, cex = 1)
axis(1, lwd = 3, at = c (1, 20, 60, 80), labels = c("-2kb", "TSS", "TES", "+2kb"), cex.axis = 1.3)
axis(2, at = c(0.2, 0.90), las = 2, lwd = 3, cex.axis = 1.3)
abline(v = c(20, 60), lty = 'longdash', lwd = 1.4)

legend("top", c( "wt_rep1", "wt_rep2", "mutantA_rep1", "mutantA_rep2" ), title.col = "black", col = c("red4", "red1", "steelblue4", "steelblue1"), inset=0.01, cex = 0.6, text.col = c("red4", "red1", "steelblue4", "steelblue1"), lty = c(1), lwd=c(4), pch = c(-1), bty="n")


## Plots a boxplot with means in as asterisks

boxplot(rowMeans(TRACK_matrix_te1[, 21:61], na.rm = TRUE), rowMeans(TRACK_matrix_te2[, 21:61], na.rm = TRUE), rowMeans(TRACK_matrix_te3[, 21:61], na.rm = TRUE), rowMeans(TRACK_matrix_te4[, 21:61], na.rm = TRUE), main = main_title_te, ylab = "DNA methylation ratio (CG)", names = c("wt_rep1", "wt_rep2", "mutantA_rep1", "mutantA_rep2"),  col = c("red4", "red1", "steelblue4", "steelblue1"), outline = FALSE, cex.lab = 1, notch = T, las = 2)

mean.value <- c(  mean(rowMeans(TRACK_matrix_te1[, 21:61], na.rm = TRUE), na.rm = TRUE), mean(rowMeans(TRACK_matrix_te2[, 21:61], na.rm = TRUE), na.rm = TRUE), mean(rowMeans(TRACK_matrix_te3[, 21:61], na.rm = TRUE), na.rm = TRUE), mean(rowMeans(TRACK_matrix_te4[, 21:61], na.rm = TRUE), na.rm = TRUE)  )

points(seq(1:4), mean.value, pch = "*", col = "gold", cex = 2) 
