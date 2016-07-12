mq <- read.csv("mq",header=T,sep=",")
mq_matrix <- data.matrix(mq)
mq_headtmap <- heatmap(mq_matrix, Rowv=NA, Colv=NA, col =brewer.pal(9,"Blues"), scale="column", margins=c(5,10))
