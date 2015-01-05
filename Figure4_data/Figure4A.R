x = c(0.975, 0.71, 0.51, 0.49)
par(mar=c(7,20,6,1))

barplot(x, width = 1, names.arg = c("Genes with isoform-specific sgRNA sites", "Transcripts with isoform-specific sgRNA sites", "Sites in the sense strand", "Sites in the anti-sense strand"), horiz = TRUE, col = c(1,2,3,4), las=1, cex.names = 1.1, cex.axis = 1.2, xlim = c(0, 1))
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = rgb(1, 1, 1, alpha=0.2))
