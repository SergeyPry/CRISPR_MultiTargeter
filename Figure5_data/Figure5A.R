x = c(0.382, 0.439, 0.561, 0.571, 0.429)
par(mar=c(7,18,6,1))

barplot(x, width = 1, names.arg = c("Pairs with common targets", "Sites with a mismatch", "Sites without a mismatch","Sites in the sense strand", "Sites in the anti-sense strand"), horiz = TRUE, col = c(1,2,3,4,5), las=1, cex.names = 1.1, cex.axis = 1.2, xlim = c(0, 0.6))
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = rgb(1, 1, 1, alpha=0.2))