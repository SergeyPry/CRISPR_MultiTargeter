library(ggplot2)

data = read.csv("total_unique_sites.txt")

m <- ggplot(data, aes(x=Frequency))
m + geom_histogram(aes(fill = ..count..), binwidth = 0.15)+scale_x_sqrt(breaks=c(10,50,100,200,400,600,800,1000,1250, 1500,1750, 2000, 2500))+scale_y_sqrt(breaks=c(10,100,250,500,750,1000,1250,1500,2000,2500))+scale_fill_gradient("Count", low = "Orange", high = "Blue")+geom_hline(yintercept=48.74, linetype="dashed")+ylab("Frequency")+xlab("Number of transcript isoform-specific sgRNA target sites")+guides(fill=guide_legend(title="Frequency"))
ggsave(file="Figure4B.tiff", width=12, height=5)
