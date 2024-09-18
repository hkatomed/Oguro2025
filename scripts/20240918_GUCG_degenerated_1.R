## This file is based on 20230805_Asano_motif_1.txt (6th trial)

HOMEDIR <- "xxx/yyy/GUCG"	# Define home directory
setwd(HOMEDIR)
setwd("XLSX")
library(openxlsx)
AsanoDF <- read.xlsx(xlsxFile = "20230728_Asano_motif.xlsx", colNames = FALSE)
names(AsanoDF) <- c("Gene", "CDS_51bp", "GUCG", "Name")


AsanoDF[which(nchar(AsanoDF$CDS_51bp) != 51),]

# > AsanoDF[which(nchar(AsanoDF$CDS_51bp) != 51),]
#     Gene                                             CDS_51bp GUCG    Name
# 18  AHP1   ATGTCTGACTTAGTTAACAAGAAATTCCCAGCTGGCGACTACAAATTCCA    3 YLR109W
# 19  GND1   ATGTCTGCTGATTTCGGTTTGATTGGTTTGGCCGTCATGGGTCAAAATTT    5 YHR183W
# 22 CDC19   ATGTCTAGATTAGAAAGATTGACCTCATTAAACGTTGTTGCTGGTTCTGA    2 YAL038W
# 26  SSE1 ATGAGTACTCCATTTGGTTTAGATTTAGGTAACAATAACTCTGTCCTTGCCG    2 YPL106C

# Add A to the 3' end
AHP1new <- "ATGTCTGACTTAGTTAACAAGAAATTCCCAGCTGGCGACTACAAATTCCAA"

# Add G to the 3' end
GND1new <- "ATGTCTGCTGATTTCGGTTTGATTGGTTTGGCCGTCATGGGTCAAAATTTG"

# Add C to the 3' end
CDC19new <- "ATGTCTAGATTAGAAAGATTGACCTCATTAAACGTTGTTGCTGGTTCTGAC"

# Remove G from the 3' end
SSE1new <- "ATGAGTACTCCATTTGGTTTAGATTTAGGTAACAATAACTCTGTCCTTGCC"

AsanoDF$CDS_51bp[18] <- AHP1new
AsanoDF$CDS_51bp[19] <- GND1new
AsanoDF$CDS_51bp[22] <- CDC19new
AsanoDF$CDS_51bp[26] <- SSE1new

AsanoDF[c(18, 19, 22, 26),]

# > AsanoDF[c(18, 19, 22, 26),]
#     Gene                                            CDS_51bp GUCG    Name
# 18  AHP1 ATGTCTGACTTAGTTAACAAGAAATTCCCAGCTGGCGACTACAAATTCCAA    3 YLR109W
# 19  GND1 ATGTCTGCTGATTTCGGTTTGATTGGTTTGGCCGTCATGGGTCAAAATTTG    5 YHR183W
# 22 CDC19 ATGTCTAGATTAGAAAGATTGACCTCATTAAACGTTGTTGCTGGTTCTGAC    2 YAL038W
# 26  SSE1 ATGAGTACTCCATTTGGTTTAGATTTAGGTAACAATAACTCTGTCCTTGCC    2 YPL106C

# Substitute U with T
for(i in 1:nrow(AsanoDF)){
	SEQ <- AsanoDF$CDS_51bp[i]
	if(length(grep(pattern = "U", x = SEQ)) > 0){
		cat("# ", i, " U: ", AsanoDF$Gene[i], "\n", sep = "")
	}
}

# 41 U: TSA1
# 43 U: SOD1

library(Biostrings)
AsanoDF$CDS_51bp[41] <- as.character(DNAString(RNAString(AsanoDF$CDS_51bp[41])))
AsanoDF$CDS_51bp[43] <- as.character(DNAString(RNAString(AsanoDF$CDS_51bp[43])))

# > AsanoDF$CDS_51bp[41]
# [1] "ATGGTCGCTCAAGTTCAAAAGCAAGCTCCAACTTTTAAGAAAACTGCCGTC"
# > AsanoDF$CDS_51bp[43]
# [1] "ATGGTTCAAGCAGTCGCAGTGTTAAAGGGTGATGCCGGTGTCTCTGGTGTT"

# Search duplicated names
dupDF <- data.frame(Gene = AsanoDF$Gene, count = 0)
for(i in 1:length(AsanoDF$Gene)){
	dupDF$count[dupDF$Gene == AsanoDF$Gene[i]] <- 
		dupDF$count[dupDF$Gene == AsanoDF$Gene[i]] + 1
}

dupDF[dupDF$count > 1,]

# > dupDF[dupDF$count > 1,]
#                Gene count
# 34  Uncharacterized     4
# 66  uncharacterized     2
# 80  Uncharacterized     4
# 162 Uncharacterized     4
# 170 uncharacterized     2
# 204 Uncharacterized     4

AsanoDF[c(34, 66, 80, 162, 170, 204),]

# > AsanoDF[c(34, 66, 80, 162, 170, 204),]
#                Gene                                            CDS_51bp GUCG      Name
# 34  Uncharacterized ATGAACGTTTCTTTTAGGTCAAACTATCACCAAGGTACAATCCGATGTAGT    1   YDR278C
# 66  uncharacterized ATGCAATTCAAGACCATCGTCGCTGCCTTCGCTACTGTTGCCGCTGTCCAA    0 YDR524C-B
# 80  Uncharacterized ATGCGTGAGCAATTGAAGCTTTTTACGAGGGAAATAGTCGATTTTACATTT    2   YNL143C
# 162 Uncharacterized ATGTCATCTGCTCTATACAAACAAAGCACAAATTTTACTCATTCTACCGGT    1 YBR085C-A
# 170 uncharacterized ATGCATGAGGTTACCCGCACTTATTATTTTTTTCTCTTTTTTTTTCTTAGC    0   YGL088W
# 204 Uncharacterized ATGAAGATAAAAATTTCCATCGAAATTAGTCTTTCGCTCCTATCTGAACAT    0 YJL047C-A

AsanoDF$Gene[34] <- AsanoDF$Name[34]
AsanoDF$Gene[66] <- AsanoDF$Name[66]
AsanoDF$Gene[80] <- AsanoDF$Name[80]
AsanoDF$Gene[162] <- AsanoDF$Name[162]
AsanoDF$Gene[170] <- AsanoDF$Name[170]
AsanoDF$Gene[204] <- AsanoDF$Name[204]

# > AsanoDF[c(34, 66, 80, 162, 170, 204),]
#          Gene                                            CDS_51bp GUCG      Name
# 34    YDR278C ATGAACGTTTCTTTTAGGTCAAACTATCACCAAGGTACAATCCGATGTAGT    1   YDR278C
# 66  YDR524C-B ATGCAATTCAAGACCATCGTCGCTGCCTTCGCTACTGTTGCCGCTGTCCAA    0 YDR524C-B
# 80    YNL143C ATGCGTGAGCAATTGAAGCTTTTTACGAGGGAAATAGTCGATTTTACATTT    2   YNL143C
# 162 YBR085C-A ATGTCATCTGCTCTATACAAACAAAGCACAAATTTTACTCATTCTACCGGT    1 YBR085C-A
# 170   YGL088W ATGCATGAGGTTACCCGCACTTATTATTTTTTTCTCTTTTTTTTTCTTAGC    0   YGL088W
# 204 YJL047C-A ATGAAGATAAAAATTTCCATCGAAATTAGTCTTTCGCTCCTATCTGAACAT    0 YJL047C-A

# > AsanoDF[123,]
#               Gene                                            CDS_51bp GUCG    Name
# 123 (mitochodrial) ATGTCTGCAAACGAATTCTACTCAAGTGGCCAACAAGGTCAATATAACCAG    1 YNL208W

AsanoDF$Gene[123] <- AsanoDF$Name[123]

# > AsanoDF[123,]
#        Gene                                            CDS_51bp GUCG    Name
# 123 YNL208W ATGTCTGCAAACGAATTCTACTCAAGTGGCCAACAAGGTCAATATAACCAG    1 YNL208W

## No duplication confirmed
# > length(AsanoDF$Gene)
# [1] 246
# > length(unique(AsanoDF$Gene))
# [1] 246

# > nrow(AsanoDF)
# [1] 246

##############################


#################

## Obtain 51 bp sequence of SSA1
SEQ <- AsanoDF$CDS_51bp[1]

# > SEQ
# [1] "ATGTCAAAAGCTGTCGGTATTGATTTAGGTACAACATACTCGTGTGTTGCT"


library(Biostrings)
etwd(HOMEDIR)
setwd("TXT")

SELECTED <- readLines(con = "20230805_SELECTED.txt")		# 9 motifs removed
SELECTEDm1 <- readLines(con = "20230805_SELECTEDm1.txt")	# no change

## 20230801_Asano_motif_1.txt の結果（９モチーフ削除前）
# > SELECTED
#  [1] "AGYYGNN" "NGYYGRN" "NNYYGRY" "NGYYGRY" "NAYYGRY" "NGYYARY" "NAYYARY" "AGYYGRN" "GAYYGRN"
# [10] "AGYYGRY" "AGYYARY" "GAYYGRY" "GGYYARY" "GGYYGGN" "AAYYGGN" "AGYYAGN" "GGYYAGN" "AAYYAGN"
# [19] "GGYYGGY" "AAYYGGY" "AAYYAGY" "GNYYGAY" "GGNYGAY" "GGYNGAY" "GNYYGGY" "GGNYGGY" "GGYNGGY"
# [28] "GGYYNGY" "GGYYGNY" "AANYGGY" "AAYNGGY" "AAYYNGY" "AAYYGNY" "AGNYAAY" "AGYNAAY" "ANYYAGY"
# [37] "AGNYAGY" "AGYNAGY" "GANYGAY" "GAYNGAY" "GAYYNAY" "GANYGGY" "GAYNGGY" "GAYYNGY" "GAYYGNY"
# [46] "GGYYANY" "AANYAGY" "AAYNAGY" "AAYYANY"

## 今回（９モチーフ削除後）
# > SELECTED
#  [1] "AGYYGNN" "NGYYGRN" "NNYYGRY" "NGYYGRY" "NAYYGRY" "NGYYARY" "NAYYARY" "AGYYGRN" "GAYYGRN"
# [10] "AGYYGRY" "AGYYARY" "GAYYGRY" "GGYYARY" "AAYYGGN" "AGYYAGN" "GGYYAGN" "AAYYAGN" "GGYYGGY"
# [19] "AAYYGGY" "AAYYAGY" "GNYYGAY" "GGNYGAY" "GGYNGAY" "GNYYGGY" "GGYYGNY" "AANYGGY" "AAYNGGY"
# [28] "AAYYNGY" "AAYYGNY" "AGNYAAY" "AGYNAAY" "ANYYAGY" "AGNYAGY" "AGYNAGY" "GANYGAY" "GANYGGY"
# [37] "GAYNGGY" "AANYAGY" "AAYNAGY" "AAYYANY"

# > unique(SELECTED)
#  [1] "AGYYGNN" "NGYYGRN" "NNYYGRY" "NGYYGRY" "NAYYGRY" "NGYYARY" "NAYYARY" "AGYYGRN" "GAYYGRN"
# [10] "AGYYGRY" "AGYYARY" "GAYYGRY" "GGYYARY" "AAYYGGN" "AGYYAGN" "GGYYAGN" "AAYYAGN" "GGYYGGY"
# [19] "AAYYGGY" "AAYYAGY" "GNYYGAY" "GGNYGAY" "GGYNGAY" "GNYYGGY" "GGYYGNY" "AANYGGY" "AAYNGGY"
# [28] "AAYYNGY" "AAYYGNY" "AGNYAAY" "AGYNAAY" "ANYYAGY" "AGNYAGY" "AGYNAGY" "GANYGAY" "GANYGGY"
# [37] "GAYNGGY" "AANYAGY" "AAYNAGY" "AAYYANY"


# > SELECTEDm1
# [1] "AGYYGRY"


## 1: no mismatch
PATTERNS0 <- SELECTED

## 2: allow one mismatch
PATTERNS1 <- SELECTEDm1

######################################

matchList <- list()
for(i in 1:nrow(AsanoDF)){
	GENE <- AsanoDF$Gene[i]
	matchList[[GENE]] <- integer(length = 51)
	SEQ <- AsanoDF$CDS_51bp[i]
	
	# U に対応
	if(length(grep(pattern = "U", SEQ)) > 0) SEQ <- RNAString(SEQ)

# matchIndex <- integer(length = 51)
for(j in 1:length(PATTERNS0)){
	PATTERN <- PATTERNS0[j]
	temp <- matchPattern(pattern = PATTERN, subject = DNAString(SEQ), 
			fixed = FALSE, max.mismatch = 0)
	if(length(temp@ranges@start) > 0){
		matchList[[GENE]][temp@ranges@start] <- 1
	}
}

for(j in 1:length(PATTERNS1)){
	PATTERN <- PATTERNS1[j]
	temp <- matchPattern(pattern = PATTERN, subject = DNAString(SEQ), 
			fixed = FALSE, max.mismatch = 1)
	if(length(temp@ranges@start) > 0){
		matchList[[GENE]][temp@ranges@start] <- 1
	}
}
}

length(matchList)

# > length(matchList)
# [1] 246


## For checking hit positions in sequences
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

## save
setwd(HOMEDIR)
setwd("XLSX")
write.xlsx(compareDF, file = "20230805_AGYYGRY_6thTrial.xlsx")


######################################

matchDF <- AsanoDF
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


## save
setwd(HOMEDIR)
setwd("XLSX")
write.xlsx(matchDF, file = "20230805_matchDF_6thTrial.xlsx")


#################################################################
## Motif enrichment
geneNumber <- 20
colNames <- c(PATTERNS0, paste(PATTERNS1, "m1", sep = ""))

# > colNames
#  [1] "AGYYGNN"   "NGYYGRN"   "NNYYGRY"   "NGYYGRY"   "NAYYGRY"   "NGYYARY"   "NAYYARY"   "AGYYGRN"  
#  [9] "GAYYGRN"   "AGYYGRY"   "AGYYARY"   "GAYYGRY"   "GGYYARY"   "AAYYGGN"   "AGYYAGN"   "GGYYAGN"  
# [17] "AAYYAGN"   "GGYYGGY"   "AAYYGGY"   "AAYYAGY"   "GNYYGAY"   "GGNYGAY"   "GGYNGAY"   "GNYYGGY"  
# [25] "GGYYGNY"   "AANYGGY"   "AAYNGGY"   "AAYYNGY"   "AAYYGNY"   "AGNYAAY"   "AGYNAAY"   "ANYYAGY"  
# [33] "AGNYAGY"   "AGYNAGY"   "GANYGAY"   "GANYGGY"   "GAYNGGY"   "AANYAGY"   "AAYNAGY"   "AAYYANY"  
# [41] "AGYYGRYm1"

# > matchDF[1:3, 1:4]
#    Gene                                            CDS_51bp GUCG    Name
# 1  SSA1 ATGTCAAAAGCTGTCGGTATTGATTTAGGTACAACATACTCGTGTGTTGCT    3 YAL005C
# 2 HSP82 ATGGCTAGTGAAACTTTTGAATTTCAAGCTGAAATTACTCAGTTGATGAGT    3 YPL240C
# 3  TDH2 ATGGTTAGAGTTGCTATTAACGGTTTCGGTAGAATCGGTAGATTGGTTATG    6 YJR009C

# > matchDF[1:3, 5:12]
#   AGYYGNN NGYYGRN NNYYGRY NGYYGRY NAYYGRY NGYYARY NAYYARY AGYYGRN
# 1       1       1       2       1       1       0       0       0
# 2       2       2       1       1       0       1       0       2
# 3       1       0       3       0       2       0       1       0

# > ncol(matchDF)
# [1] 45
# > nrow(matchDF)
# [1] 246



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

## save
setwd(HOMEDIR)
setwd("XLSX")
write.xlsx(outList, file = "20230805_motifEnrichment_6thTrial.xlsx")




## topHits >= 2 & pvalueTop < 0.05
outList2 <- outList
for(n in c(seq(10, 240, 10), 246)){
	geneNumber <- n
	sheetName <- paste(n, "_genes", sep = "")
	enrichDF <- outList[[sheetName]]
	temp <- subset(enrichDF, topHits >= 2 & pvalueTop < 0.05)
	temp <- temp[order(temp$pvalueTop),]
	outList2[[sheetName]] <- temp
}
## save
setwd(HOMEDIR)
setwd("XLSX")
write.xlsx(outList2, file = "20230805_motifEnrichment_6thTrial_p0.05_topHits2.xlsx")
