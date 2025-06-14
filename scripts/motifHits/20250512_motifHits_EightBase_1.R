## original note: 20250512_motifHits_EightBase_1.R.txt


HOMEDIR <- ""

library(openxlsx)
setwd(HOMEDIR)
setwd("data")
setwd("motifHits")
setwd("RData")
load(file = "20250427_matchDFlist_EightBase_4940.RData")	# matchDFlist
load(file = "20250427_EightBase_4940.RData")			# compareDF

geneNamesComp <- compareDF$first51$Gene

geneNamesComp <- geneNamesComp[seq(1, 9880, 2)]

geneNamesMatch <- matchDFlist$first51$Gene

length(geneNamesComp)
length(geneNamesMatch)

for(i in 1:4940){
	if(i < 6){
		cat("geneNamesComp: ", geneNamesComp[i], ", ", "geneNamesMatch: ", geneNamesMatch[i], "\n", sep = "")
	}
	if(geneNamesComp[i] != geneNamesMatch[i]){
		cat("i: ", i, "geneNamesComp[i] != geneNamesMatch[i] == TRUE")
	}
}

sheetNames <- names(matchDFlist)

compareDFplusMatchDFlist <- list()
for(s in 1:length(sheetNames)){
	sheetName <- sheetNames[s]
	cat("sheetName: ", sheetName, "\n", sep = "")
	matchDFeach <- matchDFlist[[sheetName]]
	matchDFeach$seqLine <- TRUE
	temp <- rbind(matchDFeach, matchDFeach)
	colnames(temp)[ncol(temp)-1] <- "merge"
	empty <- matchDFeach[1,]
	empty[1, 1] <- ""
	empty[1, 2:ncol(empty)] <- NA
	for(i in 1:nrow(matchDFeach)){
		temp[i*2-1,] <- matchDFeach[i,]
		temp[i*2,] <- empty
	}

	compareDFeach <- compareDF[[sheetName]]
	compareDFeach$hit2 <- compareDFeach$hit
	compareDFeach$hit[1:(nrow(compareDFeach)-1)] <- compareDFeach$hit[2:nrow(compareDFeach)]

	compareDFplusMatchDF <- cbind(compareDFeach, temp)
	compareDFplusMatchDFlist[[s]] <- compareDFplusMatchDF
	names(compareDFplusMatchDFlist)[s] <- sheetName
}

setwd(HOMEDIR)
setwd("data")
setwd("motifHits")
setwd("RData")
save(compareDFplusMatchDFlist, file = "20250512_compareDFplusMatchDFlist_EightBase.RData")

setwd("../XLSX")
write.xlsx(compareDFplusMatchDFlist, file = "20250512_compareDFplusMatchDFlist_EightBase.xlsx")


setwd(HOMEDIR)
setwd("data")
setwd("motifHits")
setwd("RData")
load(file = "20250512_compareDFplusMatchDFlist_EightBase.RData")

setwd(HOMEDIR)
setwd("data/motifHits/XLSX")
oldS4 <- read.xlsx(xlsxFile = "Table\ S4.xlsx")

S4geneNames <- oldS4$Name

temp <- rbind(oldS4, oldS4)
colnames(temp)[ncol(temp)] <- "merge"
empty <- oldS4[1,]
empty[1, 1] <- ""
empty[1, 2:ncol(empty)] <- NA
for(i in 1:nrow(oldS4)){
	temp[i*2-1,] <- oldS4[i,]
	temp[i*2,] <- empty
}

for(i in 1:nrow(oldS4)){
	temp[i*2,1] <- oldS4[i,1]
	temp[i*2,2] <- oldS4[i,2]
}

sheetNames <- names(compareDFplusMatchDFlist)
S4_246_compareDFplusMatchDFlist <- compareDFplusMatchDFlist
for(s in 1:length(sheetNames)){
	tempDF <- compareDFplusMatchDFlist[[s]][0,]
	for(i in 1:length(S4geneNames)){
		S4geneName <- S4geneNames[i]
		selectedDF <- subset(compareDFplusMatchDFlist[[s]], Gene == S4geneName)
		tempDF <- rbind(tempDF, selectedDF)
	}
	S4_246_compareDFplusMatchDFlist[[s]] <- cbind(temp, tempDF)
}


setwd(HOMEDIR)
setwd("data")
setwd("motifHits")
setwd("RData")
save(S4_246_compareDFplusMatchDFlist, file = "20250512_compareDFplusMatchDFlist_EightBase_S4_246.RData")

setwd("../XLSX")
write.xlsx(S4_246_compareDFplusMatchDFlist, file = "20250512_compareDFplusMatchDFlist_EightBase_S4_246.xlsx")

