## original note: 	20250521_tif34_peaks_1.R.txt
## 			20250530_tif34_peaks_1.R.txt


library(openxlsx)

HOMEDIR <- ""

setwd(HOMEDIR)
setwd("data")
setwd("motifHits")
setwd("XLSX")
threeEndPeaks <- read.xlsx(xlsxFile = "20250519_iRPF_3endPeaks-ka.xlsx", colNames = TRUE)

names(threeEndPeaks)[1] <- "Gene"

EightBase_hits <- read.xlsx(xlsxFile = "20250427_EightBase_4940.xlsx", colNames = TRUE)

SEQS <-  EightBase_hits[seq(1, nrow(EightBase_hits)-1, 2),]$Seq_Match
EightBase_hits <- EightBase_hits[seq(2, nrow(EightBase_hits), 2),]

threeEndPeaks$sequence <- as.character(NA)
threeEndPeaks$Seq_Match <- as.character(NA)
threeEndPeaks$hit <- as.integer(NA)
threeEndPeaks$motif5ends <- as.character(NA)
threeEndPeaks$motif4th <- as.character(NA)
threeEndPeaks$relPos <- as.character(NA)
threeEndPeaks$nearestPos <- as.integer(NA)

for(i in 1:46){					
	geneName <- threeEndPeaks$Name[i]
	targetIndex <- which(EightBase_hits$Gene == geneName)
	threeEndPeaks$hit[i] <- as.integer(EightBase_hits$hit[targetIndex])

	SEQ <- SEQS[targetIndex]
	threeEndPeaks$sequence[i] <- SEQ
	hits <- EightBase_hits$Seq_Match[targetIndex]
	threeEndPeaks$Seq_Match[i] <- hits

	motif5ends <- which(as.integer(strsplit(hits, split = "")[[1]]) == 1)
	threeEndPeaks$motif5ends[i] <- paste(motif5ends, collapse = ",", sep = "")
	motif4th <- motif5ends + 3
	threeEndPeaks$motif4th[i] <- paste(motif4th, collapse = ",", sep = "")

	peakPos <- threeEndPeaks$from_ATG[i]
	relPos <- motif4th - peakPos
	threeEndPeaks$relPos[i] <- paste(relPos, collapse = ",", sep = "")
	threeEndPeaks$nearestPos[i] <- relPos[which.min(abs(relPos))]
}

threeEndPeaks.noNA <- subset(threeEndPeaks, is.na(nearestPos) == FALSE)

dep <- subset(threeEndPeaks.noNA, tif34_dep == TRUE)
indep <- subset(threeEndPeaks.noNA, tif34_indep == TRUE)

FROM <- min(threeEndPeaks$nearestPos, na.rm = TRUE) - 1
TO <- max(threeEndPeaks$nearestPos, na.rm = TRUE) + 1
motifCounts <- data.frame(pos = seq(FROM, TO, 1), dep = 0, indep = 0)

for(i in 1:nrow(dep)){
	nearestPos <- dep$nearestPos[i]
	motifCounts$dep[motifCounts$pos == nearestPos] <- 
		motifCounts$dep[motifCounts$pos == nearestPos] + 1
}

for(i in 1:nrow(indep)){
	nearestPos <- indep$nearestPos[i]
	motifCounts$indep[motifCounts$pos == nearestPos] <- 
		motifCounts$indep[motifCounts$pos == nearestPos] + 1
}

outList <- list(threeEndPeaks = threeEndPeaks, motifCounts = motifCounts)

setwd(HOMEDIR)
setwd("data")
setwd("motifHits")
setwd("XLSX")
write.xlsx(outList, file = "20250521_iRPF_3endPeaks_relPos_add5genes.xlsx")


## 20250530_tif34_peaks_1.R.txt from here ####################

HOMEDIR <- ""

library(openxlsx)
setwd(HOMEDIR)
setwd("data")
setwd("motifHits")
setwd("XLSX")
motifCounts <- read.xlsx(xlsxFile = "20250521_iRPF_3endPeaks_relPos_add5genes.xlsx", sheet = "motifCounts")

library(ggplot2)

df.dep <- motifCounts[,c(1,2)]
df.dep$tif34 <- "dep"
colnames(df.dep)[2] <- "count"

df.indep <- motifCounts[,c(1,3)]
df.indep$tif34 <- "indep"
colnames(df.indep)[2] <- "count"

df2 <- rbind(df.dep, df.indep)

setwd(HOMEDIR)
setwd("data")
setwd("motifHits")
setwd("PDF")
pdf(file = "20250530_tif34_peaks_1.pdf", width = 4, height = 1.8)
ggplot(data = df2, aes(x = pos, y = count, fill = tif34)) + geom_bar(stat = "identity")
dev.off()



