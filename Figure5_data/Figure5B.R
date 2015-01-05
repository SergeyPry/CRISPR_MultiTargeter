library(ggplot2)

data = read.csv("ohnologs_total_sites.txt")

m <- ggplot(data, aes(x=Frequency))
m + geom_histogram(aes(fill = ..count..), binwidth = 0.13)+scale_x_sqrt(breaks=c(10,50,100,150,200, 250, 300,350))+scale_y_sqrt(breaks=c(10,100,200,300,400,500,600))+scale_fill_gradient("Count", low = "Orange", high = "Blue")+geom_hline(yintercept=6.47, linetype="dashed")+ylab("Frequency")+xlab('Number of common sgRNA target sites')+guides(fill=guide_legend(title="Frequency"))

ggsave(file="Figure5B.tiff", width=12, height=5)
