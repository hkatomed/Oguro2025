# Daniel Riordan
# Copyright 2010
# driordan@stanford.edu
#
# subsethyperpvals.R
# 
# This program takes in files with tab-delimited text format:
# name	a	b	c	d
# 
# and returns P(A >= a) = P(A > a-1)
# NOTE: ** ASSUMES a/b is a subset included in c/d **
#

data <- read.table("tempfile.txt")
logp <- phyper((data[,2]-1), data[,4], data[,5]-data[,4], data[,3], F, T)/log(10)
data <- cbind(data,logp)
write.table(data, file ="tempfile.p", quote=F, sep="\t", row.names=F, col.names=F)

rm(data)
rm(logp)