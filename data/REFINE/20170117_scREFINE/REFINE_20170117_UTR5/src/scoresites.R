# Daniel Riordan
# Copyright 2010
# driordan@stanford.edu
#
# scoresites.R

data <- read.table("tempfile.txt")
scores <- unique(data[data[,6]==1,2])
logp <- sapply(scores, function(y){phyper(sum(data[data[,2]>=y,6]), sum(data[,6]), length(data[,6])-sum(data[,6]), length(data[data[,2]>=y,6]), lower.tail=F, log.p=T)/log(10)})
outdata <- cbind(scores, logp)
write.table(outdata, file ="tempfile.ss", quote=F, sep="\t", row.names=F, col.names=F)

rm(data)
rm(scores)
rm(logp)
rm(outdata)