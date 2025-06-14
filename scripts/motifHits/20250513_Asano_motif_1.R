## original note: 20250513_Asano_motif_1.R.txt


library(openxlsx)
library(Biostrings)


HOMEDIR <- ""

setwd(HOMEDIR)
setwd("data")
setwd("motifHits")
setwd("XLSX")
AsanoDF <- read.xlsx(xlsxFile = "20230728_Asano_motif.xlsx", colNames = FALSE)
names(AsanoDF) <- c("Gene", "CDS_51bp", "GUCG", "Name")

which(nchar(AsanoDF$CDS_51bp) != 51)

AsanoDF[which(nchar(AsanoDF$CDS_51bp) != 51),]

AHP1new <- "ATGTCTGACTTAGTTAACAAGAAATTCCCAGCTGGCGACTACAAATTCCAA"
GND1new <- "ATGTCTGCTGATTTCGGTTTGATTGGTTTGGCCGTCATGGGTCAAAATTTG"
CDC19new <- "ATGTCTAGATTAGAAAGATTGACCTCATTAAACGTTGTTGCTGGTTCTGAC"
SSE1new <- "ATGAGTACTCCATTTGGTTTAGATTTAGGTAACAATAACTCTGTCCTTGCC"

AsanoDF$CDS_51bp[18] <- AHP1new
AsanoDF$CDS_51bp[19] <- GND1new
AsanoDF$CDS_51bp[22] <- CDC19new
AsanoDF$CDS_51bp[26] <- SSE1new

AsanoDF[c(18, 19, 22, 26),]

for(i in 1:nrow(AsanoDF)){
	SEQ <- AsanoDF$CDS_51bp[i]
	if(length(grep(pattern = "U", x = SEQ)) > 0){
		cat("# ", i, " U: ", AsanoDF$Gene[i], "\n", sep = "")
	}
}

AsanoDF$CDS_51bp[41] <- as.character(DNAString(RNAString(AsanoDF$CDS_51bp[41])))
AsanoDF$CDS_51bp[43] <- as.character(DNAString(RNAString(AsanoDF$CDS_51bp[43])))

dupDF <- data.frame(Gene = AsanoDF$Gene, count = 0)
for(i in 1:length(AsanoDF$Gene)){
	dupDF$count[dupDF$Gene == AsanoDF$Gene[i]] <- 
		dupDF$count[dupDF$Gene == AsanoDF$Gene[i]] + 1
}

dupDF[dupDF$count > 1,]

AsanoDF[c(34, 66, 80, 162, 170, 204),]

AsanoDF$Gene[34] <- AsanoDF$Name[34]
AsanoDF$Gene[66] <- AsanoDF$Name[66]
AsanoDF$Gene[80] <- AsanoDF$Name[80]
AsanoDF$Gene[162] <- AsanoDF$Name[162]
AsanoDF$Gene[170] <- AsanoDF$Name[170]
AsanoDF$Gene[204] <- AsanoDF$Name[204]

AsanoDF$Gene[123] <- AsanoDF$Name[123]

SEQ <- AsanoDF$CDS_51bp[1]

setwd(HOMEDIR)
setwd("data")
setwd("motifHits")
setwd("TXT")

PATTERNS0 <- readLines(con = "20240427_Final8base_m0.txt")
PATTERNS1 <- ""
PATTERNS2 <- readLines(con = "20240427_Final8base_m2.txt")


matchList <- list()
for(i in 1:nrow(AsanoDF)){
	GENE <- AsanoDF$Gene[i]
	matchList[[GENE]] <- integer(length = 51)
	SEQ <- AsanoDF$CDS_51bp[i]
	if(length(grep(pattern = "U", SEQ)) > 0) SEQ <- RNAString(SEQ)

if(PATTERNS0[1] != ""){
for(j in 1:length(PATTERNS0)){
	PATTERN <- PATTERNS0[j]
	temp <- matchPattern(pattern = PATTERN, subject = DNAString(SEQ), 
			fixed = FALSE, max.mismatch = 0)
	if(length(temp@ranges@start) > 0){
		matchList[[GENE]][temp@ranges@start] <- 1
	}
}
}

if(PATTERNS1[1] != ""){
for(j in 1:length(PATTERNS1)){
	PATTERN <- PATTERNS1[j]
	temp <- matchPattern(pattern = PATTERN, subject = DNAString(SEQ), 
			fixed = FALSE, max.mismatch = 1)
	if(length(temp@ranges@start) > 0){
		matchList[[GENE]][temp@ranges@start] <- 1
	}
}
}

if(PATTERNS2[1] != ""){
for(j in 1:length(PATTERNS2)){
	PATTERN <- PATTERNS2[j]
	temp <- matchPattern(pattern = PATTERN, subject = DNAString(SEQ), 
			fixed = FALSE, max.mismatch = 2)
	if(length(temp@ranges@start) > 0){
		matchList[[GENE]][temp@ranges@start] <- 1
	}
}
}
}

length(matchList)

matchLines <- character(length = length(matchList) *2)
geneNames <- character(length = length(matchList) *2)
hit <- character(length = length(matchList) *2)
for(i in 1:length(matchList)){
	temp <- paste(matchList[[i]], collapse = "")
	matchLines[(i-1)*2+1] <- AsanoDF$CDS_51bp[i]
	matchLines[(i-1)*2+2] <- gsub(pattern = "0", replacement = " ", temp)
	geneNames[(i-1)*2+1] <- AsanoDF$Gene[i]
	geneNames[(i-1)*2+2] <- AsanoDF$Gene[i]
	if(sum(matchList[[i]]) > 0){
		hit[(i-1)*2+2] <- as.character(sum(matchList[[i]]))
	}
}

compareDF <- data.frame(Gene = geneNames, Seq_Match = matchLines, hit = hit)

setwd(HOMEDIR)
setwd("data")
setwd("motifHits")
setwd("XLSX")
write.xlsx(compareDF, file = "20250513_EightBase.xlsx")


matchDF <- AsanoDF

if(PATTERNS0[1] != ""){
for(j in 1:length(PATTERNS0)){
	PATTERN <- PATTERNS0[j]
	matchDF[[PATTERN]] <- as.integer(0)
	for(i in 1:nrow(AsanoDF)){
		SEQ <- AsanoDF$CDS_51bp[i]
		temp <- matchPattern(pattern = PATTERN, subject = DNAString(SEQ), 
				fixed = FALSE, max.mismatch = 0)
		if(length(temp@ranges@start) > 0){
			matchDF[[PATTERN]][i] <- length(temp@ranges@start)
		}
	}
}
}

if(PATTERNS1[1] != ""){
for(j in 1:length(PATTERNS1)){
	PATTERN <- PATTERNS1[j]
	colName <- paste(PATTERN, "m1", sep = "")
	matchDF[[colName]] <- as.integer(0)
	for(i in 1:nrow(AsanoDF)){
		SEQ <- AsanoDF$CDS_51bp[i]
		temp <- matchPattern(pattern = PATTERN, subject = DNAString(SEQ), 
				fixed = FALSE, max.mismatch = 1)
		if(length(temp@ranges@start) > 0){
			matchDF[[colName]][i] <- length(temp@ranges@start)
		}
	}
}
}

if(PATTERNS2[1] != ""){
for(j in 1:length(PATTERNS2)){
	PATTERN <- PATTERNS2[j]
	colName <- paste(PATTERN, "m2", sep = "")
	matchDF[[colName]] <- as.integer(0)
	for(i in 1:nrow(AsanoDF)){
		SEQ <- AsanoDF$CDS_51bp[i]
		temp <- matchPattern(pattern = PATTERN, subject = DNAString(SEQ), 
				fixed = FALSE, max.mismatch = 2)
		if(length(temp@ranges@start) > 0){
			matchDF[[colName]][i] <- length(temp@ranges@start)
		}
	}
}
}

setwd(HOMEDIR)
setwd("data")
setwd("motifHits")
setwd("XLSX")
write.xlsx(matchDF, file = "20250513_matchDF_EightBase.xlsx")


geneNumber <- 20
# colNames <- c(PATTERNS0, paste(PATTERNS1, "m1", sep = ""))
colNames <- PATTERNS0
if(PATTERNS1[1] != "") colNames <- c(colNames, paste(PATTERNS1, "m1", sep = ""))
if(PATTERNS2[1] != "") colNames <- c(colNames, paste(PATTERNS2, "m2", sep = ""))

outList <- list()
for(n in c(seq(10, 240, 10), 246)){
	geneNumber <- n
	enrichDF <- data.frame(PATTERN = colNames)
	enrichDF$geneNumber <- as.integer(NA)
	enrichDF$totalHits <- as.integer(NA)
	enrichDF$topHits <- as.integer(NA)
	enrichDF$bottomHits <- as.integer(NA)
	enrichDF$expectHits <- as.numeric(NA)
	enrichDF$topEnrichment <- as.numeric(NA)
	enrichDF$bottomEnrichment <- as.numeric(NA)
	enrichDF$pvalueTop <- as.numeric(NA)
	enrichDF$pvalueBottom <- as.numeric(NA)
	
	for(i in 1:length(colNames)){
		colName <- colNames[i]
		totalHits <- sum(matchDF[colName])
		topHits <- sum(head(matchDF[[colName]], n = geneNumber))
		bottomHits <- sum(tail(matchDF[[colName]], n = geneNumber))
		meanHits <- totalHits / nrow(matchDF)
		expectHits <- meanHits * geneNumber
		topEnrichment <- topHits / expectHits
		bottomEnrichment <- bottomHits / expectHits
		enrichDF$geneNumber[i] <- geneNumber
		enrichDF$totalHits[i] <- totalHits
		enrichDF$topHits[i] <- topHits
		enrichDF$bottomHits[i] <- bottomHits
		enrichDF$expectHits[i] <- round(expectHits, digits = 6)
		enrichDF$topEnrichment[i] <- round(topEnrichment, digits = 6)
		enrichDF$bottomEnrichment[i] <- round(bottomEnrichment, digits = 6)
		enrichDF$pvalueTop[i] <- ppois(q = topHits, lambda = expectHits, lower.tail = FALSE)
		enrichDF$pvalueBottom[i] <- ppois(q = bottomHits, lambda = expectHits, lower.tail = FALSE)
	}
	sheetName <- paste(n, "_genes", sep = "")
	outList[[sheetName]] <- enrichDF
}


setwd(HOMEDIR)
setwd("data")
setwd("motifHits")
setwd("XLSX")
write.xlsx(outList, file = "20250513_motifEnrichment_EightBase.xlsx")



outList2 <- outList
for(n in c(seq(10, 240, 10), 246)){
	geneNumber <- n
	sheetName <- paste(n, "_genes", sep = "")
	enrichDF <- outList[[sheetName]]
	temp <- subset(enrichDF, topHits >= 2 & pvalueTop < 0.05)
	temp <- temp[order(temp$pvalueTop),]
	outList2[[sheetName]] <- temp
}

setwd(HOMEDIR)
setwd("data")
setwd("motifHits")
setwd("XLSX")
write.xlsx(outList2, file = "20250513_motifEnrichment_EightBase_p0.05_topHits2.xlsx")


library(openxlsx)
setwd(HOMEDIR)
setwd("data")
setwd("motifHits")
setwd("XLSX")
inDF <- read.xlsx(xlsxFile = "20250513_EightBase.xlsx")

Index <- seq(2, 492, 2)
outDF <- inDF[Index, ]

write.xlsx(outDF, file = "20250513_EightBase_hits.xlsx")

