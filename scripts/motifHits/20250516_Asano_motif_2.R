## original note: 20250516_Asano_motif_2.R.txt

library(openxlsx)
library(Biostrings)

HOMEDIR <- ""

setwd(HOMEDIR)
setwd("data")
setwd("motifHits")
setwd("XLSX")
AsanoDF <- read.xlsx(xlsxFile = "20230728_Asano_motif.xlsx", colNames = FALSE)
names(AsanoDF) <- c("Gene", "CDS_51bp", "GUCG", "Name")

setwd(HOMEDIR)
setwd("data")
setwd("motifDistribution")
setwd("SGD_20240615")
orf_coding <- readDNAStringSet(filepath = "orf_coding.fasta.gz")

Genes6039 <- character(length = length(orf_coding))
Names6039 <- character(length = length(orf_coding))
geneNames <- names(orf_coding)
for(i in 1:length(geneNames)){
	temp <- strsplit(geneNames[i], split = " ")[[1]][1:2]
	Names6039[i] <- temp[1]
	Genes6039[i] <- temp[2]
}

AsanoDF$existNames6039 <- FALSE
AsanoDF$existGenes6039 <- FALSE
for(i in 1:nrow(AsanoDF)){
	index <- which(Names6039 == AsanoDF$Name[i])
	if(length(index) == 1)	AsanoDF$existNames6039[i] <- TRUE
	if(length(index) == 0)	cat("# i=", i, ", not found in Names6039: ", AsanoDF$Name[i], ", ", AsanoDF$Gene[i], "\n", sep = "")
	if(length(index) > 1)	cat("# i=", i, ", found many in Names6039: ", AsanoDF$Name[i], ", ", AsanoDF$Gene[i], "\n", sep = "")
	index <- which(Genes6039 == AsanoDF$Gene[i])
	if(length(index) == 1)	AsanoDF$existGenes6039[i] <- TRUE
	if(length(index) == 0)	cat("# i=", i, ", not found in Genes6039: ", AsanoDF$Gene[i], ", ", AsanoDF$Name[i], "\n", sep = "")
	if(length(index) > 1)	cat("# i=", i, ", found many in Genes6039: ", AsanoDF$Gene[i], ", ", AsanoDF$Name[i], "\n", sep = "")
}

AsanoDF$orf_coding <- as.character(NA)
AsanoDF$CDSlen <- as.integer(NA)
AsanoDF$CDSover200 <- FALSE

for(i in 1:nrow(AsanoDF)){
	index <- which(Names6039 == AsanoDF$Name[i])
	if(length(index) == 1){
		CDS <- orf_coding[[index]]
		AsanoDF$orf_coding[i] <- as.character(CDS)
		AsanoDF$CDSlen[i] <- length(CDS)
		if(length(CDS) > 200) AsanoDF$CDSover200[i] <- TRUE	## 20250516 100 を 200 に修正
	}
}

setwd(HOMEDIR)
setwd("data")
setwd("motifDistribution")
setwd("SGD_20250514")
YDR278C_CDS <- readDNAStringSet(filepath = "S288C_YDR278C_YDR278C_coding.fsa")[[1]]

i <- 34
AsanoDF$orf_coding[i] <- as.character(YDR278C_CDS)
AsanoDF$CDSlen[i] <- length(YDR278C_CDS)
if(length(YDR278C_CDS) > 200) AsanoDF$CDSover200[i] <- TRUE	## 20250516 100 を 200 に修正

AsanoDF$Gene[34] <- AsanoDF$Name[34]
AsanoDF$Gene[66] <- AsanoDF$Name[66]
AsanoDF$Gene[80] <- AsanoDF$Name[80]
AsanoDF$Gene[162] <- AsanoDF$Name[162]
AsanoDF$Gene[170] <- AsanoDF$Name[170]
AsanoDF$Gene[204] <- AsanoDF$Name[204]

AsanoDF[AsanoDF$CDSover200 == FALSE, c(1:6, 8, 9)]

AsanoDF2 <- AsanoDF[AsanoDF$CDSover200 == TRUE,]

AsanoDF2$CDS_200bp <- as.character(NA)
for(i in 1:nrow(AsanoDF2)){
	CDS <- AsanoDF2$orf_coding[i]
	AsanoDF2$CDS_200bp[i] <- substr(x = CDS, start = nchar(CDS) - 199, stop = nchar(CDS))
}

setwd(HOMEDIR)
setwd("data")
setwd("motifHits")
setwd("TXT")

PATTERNS0 <- readLines(con = "20240427_Final8base_m0.txt")
PATTERNS1 <- ""
PATTERNS2 <- readLines(con = "20240427_Final8base_m2.txt")


matchList <- list()
for(i in 1:nrow(AsanoDF2)){
	GENE <- AsanoDF2$Gene[i]
	matchList[[GENE]] <- integer(length = 200)
	SEQ <- AsanoDF2$CDS_200bp[i]

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


matchLines <- character(length = length(matchList) *2)
geneNames <- character(length = length(matchList) *2)
hit <- character(length = length(matchList) *2)
for(i in 1:length(matchList)){
	temp <- paste(matchList[[i]], collapse = "")
	matchLines[(i-1)*2+1] <- AsanoDF2$CDS_200bp[i]
	matchLines[(i-1)*2+2] <- gsub(pattern = "0", replacement = " ", temp)
	geneNames[(i-1)*2+1] <- AsanoDF2$Gene[i]
	geneNames[(i-1)*2+2] <- AsanoDF2$Gene[i]
	if(sum(matchList[[i]]) > 0){
		hit[(i-1)*2+2] <- as.character(sum(matchList[[i]]))
	}
}

compareDF <- data.frame(Gene = geneNames, Seq_Match = matchLines, hit = hit)

setwd(HOMEDIR)
setwd("data")
setwd("motifHits")
setwd("XLSX")
write.xlsx(compareDF, file = "20250516_EightBase_200base_STOP.xlsx")


matchDF <- AsanoDF2

if(PATTERNS0[1] != ""){
for(j in 1:length(PATTERNS0)){
	PATTERN <- PATTERNS0[j]
	matchDF[[PATTERN]] <- as.integer(0)
	for(i in 1:nrow(AsanoDF2)){
		SEQ <- AsanoDF2$CDS_200bp[i]
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
	for(i in 1:nrow(AsanoDF2)){
		SEQ <- AsanoDF2$CDS_200bp[i]
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
	for(i in 1:nrow(AsanoDF2)){
		SEQ <- AsanoDF2$CDS_200bp[i]
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
write.xlsx(matchDF, file = "20250516_matchDF_EightBase_200base_STOP.xlsx")

geneNumber <- 20
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
		# LAMBDA <- round(expectHits, digits = 0)
		# enrichDF$pvalueTop[i] <- ppois(q = topHits, lambda = LAMBDA, lower.tail = FALSE)
		# enrichDF$pvalueBottom[i] <- ppois(q = bottomHits, lambda = LAMBDA, lower.tail = FALSE)
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
write.xlsx(outList, file = "20250516_motifEnrichment_EightBase_200base_STOP.xlsx")


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
write.xlsx(outList2, file = "20250516_motifEnrichment_EightBase_p0.05_topHits2_200base_STOP.xlsx")


library(openxlsx)
setwd(HOMEDIR)
setwd("data")
setwd("motifHits")
setwd("XLSX")
inDF <- read.xlsx(xlsxFile = "20250516_EightBase_200base_STOP.xlsx")

Index <- seq(2, 492, 2)
outDF <- inDF[Index, ]

write.xlsx(outDF, file = "20250516_EightBase_hits_200base_STOP.xlsx")





