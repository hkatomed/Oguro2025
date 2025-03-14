
#### 6thTrial only ####################################
## finally generates 20240717_motifHits_4940.xlsx
## 20240717_6thTrial_notes.zip
## original note: 20240619_Asano_motif_1_CDSeach2.txt
## original note: 20240619_Asano_motif_2.txt
## original note: 20240717_Asano_motif_1.txt


## 20240619_Asano_motif_1_CDSeach2.txt

HOMEDIR <- ""

library(openxlsx)
library(Biostrings)

setwd(HOMEDIR)
setwd("data")
setwd("motifHits")
setwd("XLSX")
RPKMtable <- read.xlsx(xlsxFile = "Yeast_CHX_KSU_RPKM_table_2_hk2.xlsx")

setwd(HOMEDIR)
setwd("data")
setwd("motifDistribution")
setwd("SGD_20240615")
ORFs <- readDNAStringSet(filepath = "orf_coding.fasta")

orf.names <- names(ORFs)
for(i in 1:length(orf.names)){
	temp <- orf.names[i]
	orf.names[i] <- strsplit(temp, split = " ")[[1]][1]
}

names(ORFs) <- orf.names

RPKMtable$CDS_first51bp <- as.character(NA)
for(i in 1:nrow(RPKMtable)){
	geneName <- RPKMtable$Gene[i]
	index <- which(names(ORFs) == geneName)
	if(length(index) > 1) cat("i:", i, "length(index) > 1", "\n", sep = "")
	if(length(index) == 0) cat("i:", i, "length(index) == 0", "\n", sep = "")
	if(length(index) == 1){
		RPKMtable$CDS_first51bp[i] <- as.character(ORFs[[index]][1:51])
	}
}

RPKMtable$CDS_middle51bp <- as.character(NA)
for(i in 1:nrow(RPKMtable)){
	geneName <- RPKMtable$Gene[i]
	index <- which(names(ORFs) == geneName)
	if(length(index) > 1) cat("i:", i, "length(index) > 1", "\n", sep = "")
	if(length(index) == 0) cat("i:", i, "length(index) == 0", "\n", sep = "")
	if(length(index) == 1){
		CDSend <- width(ORFs)[index]
		middleStart <- round(CDSend/2) - 25
		middleEnd <- round(CDSend/2) + 25
		RPKMtable$CDS_middle51bp[i] <- as.character(ORFs[[index]][middleStart:middleEnd])
	}
}

RPKMtable$CDS_last51bp <- as.character(NA)
for(i in 1:nrow(RPKMtable)){
	geneName <- RPKMtable$Gene[i]
	index <- which(names(ORFs) == geneName)
	if(length(index) > 1) cat("i:", i, "length(index) > 1", "\n", sep = "")
	if(length(index) == 0) cat("i:", i, "length(index) == 0", "\n", sep = "")
	if(length(index) == 1){
		CDSend <- width(ORFs)[index]
		RPKMtable$CDS_last51bp[i] <- as.character(ORFs[[index]][(CDSend-50):CDSend])
	}
}

setwd(HOMEDIR)
setwd("data")
setwd("motifDistribution")
setwd("SGD_20240615")
UTR5 <- readDNAStringSet(filepath = "UTR5.fasta.gz")
UTR3 <- readDNAStringSet(filepath = "UTR3.fasta.gz")
CDS_total <- readDNAStringSet(filepath = "CDS_total.fasta.gz")
CDS_1toLast52 <- readDNAStringSet(filepath = "CDS_1toLast52.fasta.gz")
CDS_52toLast <- readDNAStringSet(filepath = "CDS_52toLast.fasta.gz")

RPKMtable$UTR5 <- as.character(NA)
for(i in 1:nrow(RPKMtable)){
	geneName <- RPKMtable$Gene[i]
	index <- which(names(UTR5) == geneName)
	if(length(index) > 1) cat("i:", i, "length(index) > 1", "\n", sep = "")
	if(length(index) == 0) cat("i:", i, "length(index) == 0", "\n", sep = "")
	if(length(index) == 1){
		CDSend <- width(UTR5)[index]
		RPKMtable$UTR5[i] <- as.character(UTR5[[index]])
	}
}

RPKMtable$UTR3 <- as.character(NA)
for(i in 1:nrow(RPKMtable)){
	geneName <- RPKMtable$Gene[i]
	index <- which(names(UTR3) == geneName)
	if(length(index) > 1) cat("i:", i, "length(index) > 1", "\n", sep = "")
	if(length(index) == 0) cat("i:", i, "length(index) == 0", "\n", sep = "")
	if(length(index) == 1){
		CDSend <- width(UTR3)[index]
		RPKMtable$UTR3[i] <- as.character(UTR3[[index]])
	}
}

RPKMtable$CDS_total <- as.character(NA)
for(i in 1:nrow(RPKMtable)){
	geneName <- RPKMtable$Gene[i]
	index <- which(names(CDS_total) == geneName)
	if(length(index) > 1) cat("i:", i, "length(index) > 1", "\n", sep = "")
	if(length(index) == 0) cat("i:", i, "length(index) == 0", "\n", sep = "")
	if(length(index) == 1){
		CDSend <- width(CDS_total)[index]
		RPKMtable$CDS_total[i] <- as.character(CDS_total[[index]])
	}
}

RPKMtable$CDS_1toLast52 <- as.character(NA)
for(i in 1:nrow(RPKMtable)){
	geneName <- RPKMtable$Gene[i]
	index <- which(names(CDS_1toLast52) == geneName)
	if(length(index) > 1) cat("i:", i, "length(index) > 1", "\n", sep = "")
	if(length(index) == 0) cat("i:", i, "length(index) == 0", "\n", sep = "")
	if(length(index) == 1){
		CDSend <- width(CDS_1toLast52)[index]
		RPKMtable$CDS_1toLast52[i] <- as.character(CDS_1toLast52[[index]])
	}
}

RPKMtable$CDS_52toLast <- as.character(NA)
for(i in 1:nrow(RPKMtable)){
	geneName <- RPKMtable$Gene[i]
	index <- which(names(CDS_52toLast) == geneName)
	if(length(index) > 1) cat("i:", i, "length(index) > 1", "\n", sep = "")
	if(length(index) == 0) cat("i:", i, "length(index) == 0", "\n", sep = "")
	if(length(index) == 1){
		CDSend <- width(CDS_52toLast)[index]
		RPKMtable$CDS_52toLast[i] <- as.character(CDS_52toLast[[index]])
	}
}

AsanoDF <- RPKMtable[,c(1, 94:101)]

setwd(HOMEDIR)
setwd("data")
setwd("motifDistribution")
setwd("TXT")

SELECTED <- readLines(con = "20230805_SELECTED.txt")		# 9 motifs removed
SELECTEDm1 <- readLines(con = "20230805_SELECTEDm1.txt")	# no change

PATTERNS0 <- SELECTED
PATTERNS1 <- SELECTEDm1


matchList <- list()
matchList[["first51"]] <- list()
matchList[["middle51"]] <- list()
matchList[["last51"]] <- list()

matchList[["UTR5"]] <- list()
matchList[["UTR3"]] <- list()
matchList[["CDS_total"]] <- list()
matchList[["CDS_1toLast52"]] <- list()
matchList[["CDS_52toLast"]] <- list()

for(i in 1:nrow(AsanoDF)){
	if(i%%100 == 0) cat(i, ", ", sep = "")
	if(i == nrow(AsanoDF)) cat("\n")
	GENE <- AsanoDF$Gene[i]

	SEQ <- AsanoDF$CDS_first51bp[i]
	if(is.na(SEQ) == TRUE){
		matchList[["first51"]][[GENE]] <- integer(length = 0)
	}
	
	# 配列がある時のみモチーフを探す
	if(is.na(SEQ) == FALSE){
		matchList[["first51"]][[GENE]] <- integer(length = nchar(SEQ))
	
		# U に対応
		if(length(grep(pattern = "U", SEQ)) > 0) SEQ <- RNAString(SEQ)
	
		# matchIndex <- integer(length = 51)
		for(j in 1:length(PATTERNS0)){
			PATTERN <- PATTERNS0[j]
			temp <- matchPattern(pattern = PATTERN, subject = DNAString(SEQ), 
					fixed = FALSE, max.mismatch = 0)
			if(length(temp@ranges@start) > 0){
				matchList[["first51"]][[GENE]][temp@ranges@start] <- 1
			}
		}
		for(j in 1:length(PATTERNS1)){
			PATTERN <- PATTERNS1[j]
			temp <- matchPattern(pattern = PATTERN, subject = DNAString(SEQ), 
					fixed = FALSE, max.mismatch = 1)
			if(length(temp@ranges@start) > 0){
				matchList[["first51"]][[GENE]][temp@ranges@start] <- 1
			}
		}
	}
	
	SEQ <- AsanoDF$CDS_middle51bp[i]
	if(is.na(SEQ) == TRUE){
		matchList[["middle51"]][[GENE]] <- integer(length = 0)
	}
		
	# 配列がある時のみモチーフを探す
	if(is.na(SEQ) == FALSE){
		matchList[["middle51"]][[GENE]] <- integer(length = nchar(SEQ))
	
		# U に対応
		if(length(grep(pattern = "U", SEQ)) > 0) SEQ <- RNAString(SEQ)

		# matchIndex <- integer(length = 51)
		for(j in 1:length(PATTERNS0)){
			PATTERN <- PATTERNS0[j]
			temp <- matchPattern(pattern = PATTERN, subject = DNAString(SEQ), 
					fixed = FALSE, max.mismatch = 0)
			if(length(temp@ranges@start) > 0){
				matchList[["middle51"]][[GENE]][temp@ranges@start] <- 1
			}
		}	
		for(j in 1:length(PATTERNS1)){
			PATTERN <- PATTERNS1[j]
			temp <- matchPattern(pattern = PATTERN, subject = DNAString(SEQ), 
					fixed = FALSE, max.mismatch = 1)
			if(length(temp@ranges@start) > 0){
				matchList[["middle51"]][[GENE]][temp@ranges@start] <- 1
			}
		}
	}

	SEQ <- AsanoDF$CDS_last51bp[i]
	if(is.na(SEQ) == TRUE){
		matchList[["last51"]][[GENE]] <- integer(length = 0)
	}
	
	# 配列がある時のみモチーフを探す
	if(is.na(SEQ) == FALSE){
		matchList[["last51"]][[GENE]] <- integer(length = nchar(SEQ))
	
		# U に対応
		if(length(grep(pattern = "U", SEQ)) > 0) SEQ <- RNAString(SEQ)

		# matchIndex <- integer(length = 51)
		for(j in 1:length(PATTERNS0)){
			PATTERN <- PATTERNS0[j]
			temp <- matchPattern(pattern = PATTERN, subject = DNAString(SEQ), 
					fixed = FALSE, max.mismatch = 0)
			if(length(temp@ranges@start) > 0){
				matchList[["last51"]][[GENE]][temp@ranges@start] <- 1
			}
		}
		for(j in 1:length(PATTERNS1)){
			PATTERN <- PATTERNS1[j]
			temp <- matchPattern(pattern = PATTERN, subject = DNAString(SEQ), 
					fixed = FALSE, max.mismatch = 1)
			if(length(temp@ranges@start) > 0){
				matchList[["last51"]][[GENE]][temp@ranges@start] <- 1
			}
		}
	}

	SEQ <- AsanoDF$UTR5[i]
	if(is.na(SEQ) == TRUE){
		matchList[["UTR5"]][[GENE]] <- integer(length = 0)
	}
	
	# 配列がある時のみモチーフを探す
	if(is.na(SEQ) == FALSE){
		matchList[["UTR5"]][[GENE]] <- integer(length = nchar(SEQ))
	
		# U に対応
		if(length(grep(pattern = "U", SEQ)) > 0) SEQ <- RNAString(SEQ)

		# matchIndex <- integer(length = 51)
		for(j in 1:length(PATTERNS0)){
			PATTERN <- PATTERNS0[j]
			temp <- matchPattern(pattern = PATTERN, subject = DNAString(SEQ), 
					fixed = FALSE, max.mismatch = 0)
			if(length(temp@ranges@start) > 0){
				matchList[["UTR5"]][[GENE]][temp@ranges@start] <- 1
			}
		}
		for(j in 1:length(PATTERNS1)){
			PATTERN <- PATTERNS1[j]
			temp <- matchPattern(pattern = PATTERN, subject = DNAString(SEQ), 
					fixed = FALSE, max.mismatch = 1)
			if(length(temp@ranges@start) > 0){
				matchList[["UTR5"]][[GENE]][temp@ranges@start] <- 1
			}
		}
	}

	SEQ <- AsanoDF$UTR3[i]
	if(is.na(SEQ) == TRUE){
		matchList[["UTR3"]][[GENE]] <- integer(length = 0)
	}
	
	# 配列がある時のみモチーフを探す
	if(is.na(SEQ) == FALSE){
		matchList[["UTR3"]][[GENE]] <- integer(length = nchar(SEQ))
	
		# U に対応
		if(length(grep(pattern = "U", SEQ)) > 0) SEQ <- RNAString(SEQ)

		# matchIndex <- integer(length = 51)
		for(j in 1:length(PATTERNS0)){
			PATTERN <- PATTERNS0[j]
			temp <- matchPattern(pattern = PATTERN, subject = DNAString(SEQ), 
					fixed = FALSE, max.mismatch = 0)
			if(length(temp@ranges@start) > 0){
				matchList[["UTR3"]][[GENE]][temp@ranges@start] <- 1
			}
		}
		for(j in 1:length(PATTERNS1)){
			PATTERN <- PATTERNS1[j]
			temp <- matchPattern(pattern = PATTERN, subject = DNAString(SEQ), 
					fixed = FALSE, max.mismatch = 1)
			if(length(temp@ranges@start) > 0){
				matchList[["UTR3"]][[GENE]][temp@ranges@start] <- 1
			}
		}
	}

	SEQ <- AsanoDF$CDS_total[i]
	if(is.na(SEQ) == TRUE){
		matchList[["CDS_total"]][[GENE]] <- integer(length = 0)
	}
	
	# 配列がある時のみモチーフを探す
	if(is.na(SEQ) == FALSE){
		matchList[["CDS_total"]][[GENE]] <- integer(length = nchar(SEQ))
	
		# U に対応
		if(length(grep(pattern = "U", SEQ)) > 0) SEQ <- RNAString(SEQ)

		# matchIndex <- integer(length = 51)
		for(j in 1:length(PATTERNS0)){
			PATTERN <- PATTERNS0[j]
			temp <- matchPattern(pattern = PATTERN, subject = DNAString(SEQ), 
					fixed = FALSE, max.mismatch = 0)
			if(length(temp@ranges@start) > 0){
				matchList[["CDS_total"]][[GENE]][temp@ranges@start] <- 1
			}
		}
		for(j in 1:length(PATTERNS1)){
			PATTERN <- PATTERNS1[j]
			temp <- matchPattern(pattern = PATTERN, subject = DNAString(SEQ), 
					fixed = FALSE, max.mismatch = 1)
			if(length(temp@ranges@start) > 0){
				matchList[["CDS_total"]][[GENE]][temp@ranges@start] <- 1
			}
		}
	}

	SEQ <- AsanoDF$CDS_1toLast52[i]
	if(is.na(SEQ) == TRUE){
		matchList[["CDS_1toLast52"]][[GENE]] <- integer(length = 0)
	}
	
	# 配列がある時のみモチーフを探す
	if(is.na(SEQ) == FALSE){
		matchList[["CDS_1toLast52"]][[GENE]] <- integer(length = nchar(SEQ))
	
		# U に対応
		if(length(grep(pattern = "U", SEQ)) > 0) SEQ <- RNAString(SEQ)

		# matchIndex <- integer(length = 51)
		for(j in 1:length(PATTERNS0)){
			PATTERN <- PATTERNS0[j]
			temp <- matchPattern(pattern = PATTERN, subject = DNAString(SEQ), 
					fixed = FALSE, max.mismatch = 0)
			if(length(temp@ranges@start) > 0){
				matchList[["CDS_1toLast52"]][[GENE]][temp@ranges@start] <- 1
			}
		}
		for(j in 1:length(PATTERNS1)){
			PATTERN <- PATTERNS1[j]
			temp <- matchPattern(pattern = PATTERN, subject = DNAString(SEQ), 
					fixed = FALSE, max.mismatch = 1)
			if(length(temp@ranges@start) > 0){
				matchList[["CDS_1toLast52"]][[GENE]][temp@ranges@start] <- 1
			}
		}
	}

	SEQ <- AsanoDF$CDS_52toLast[i]
	if(is.na(SEQ) == TRUE){
		matchList[["CDS_52toLast"]][[GENE]] <- integer(length = 0)
	}
	
	# 配列がある時のみモチーフを探す
	if(is.na(SEQ) == FALSE){
		matchList[["CDS_52toLast"]][[GENE]] <- integer(length = nchar(SEQ))
	
		# U に対応
		if(length(grep(pattern = "U", SEQ)) > 0) SEQ <- RNAString(SEQ)

		# matchIndex <- integer(length = 51)
		for(j in 1:length(PATTERNS0)){
			PATTERN <- PATTERNS0[j]
			temp <- matchPattern(pattern = PATTERN, subject = DNAString(SEQ), 
					fixed = FALSE, max.mismatch = 0)
			if(length(temp@ranges@start) > 0){
				matchList[["CDS_52toLast"]][[GENE]][temp@ranges@start] <- 1
			}
		}
		for(j in 1:length(PATTERNS1)){
			PATTERN <- PATTERNS1[j]
			temp <- matchPattern(pattern = PATTERN, subject = DNAString(SEQ), 
					fixed = FALSE, max.mismatch = 1)
			if(length(temp@ranges@start) > 0){
				matchList[["CDS_52toLast"]][[GENE]][temp@ranges@start] <- 1
			}
		}
	}
}

compareDF <- list()

matchLines <- character(length = length(matchList[["first51"]]) *2)
geneNames <- character(length = length(matchList[["first51"]]) *2)
hit <- character(length = length(matchList[["first51"]]) *2)
for(i in 1:length(matchList[["first51"]])){
	if(i%%100 == 0) cat(i, ", ", sep = "")
	if(i == length(matchList[["first51"]])) cat("\n")
	temp <- paste(matchList[["first51"]][[i]], collapse = "")
	matchLines[(i-1)*2+1] <- AsanoDF$CDS_first51bp[i]
	matchLines[(i-1)*2+2] <- gsub(pattern = "0", replacement = " ", temp)
	geneNames[(i-1)*2+1] <- AsanoDF$Gene[i]
	geneNames[(i-1)*2+2] <- AsanoDF$Gene[i]
	if(sum(matchList[["first51"]][[i]]) > 0){
		hit[(i-1)*2+2] <- as.character(sum(matchList[["first51"]][[i]]))
	}
}
compareDF[["first51"]] <- data.frame(Gene = geneNames, Seq_Match = matchLines, hit = hit)

matchLines <- character(length = length(matchList[["middle51"]]) *2)
geneNames <- character(length = length(matchList[["middle51"]]) *2)
hit <- character(length = length(matchList[["middle51"]]) *2)
for(i in 1:length(matchList[["middle51"]])){
	if(i%%100 == 0) cat(i, ", ", sep = "")
	if(i == length(matchList[["middle51"]])) cat("\n")
	temp <- paste(matchList[["middle51"]][[i]], collapse = "")
	matchLines[(i-1)*2+1] <- AsanoDF$CDS_middle51bp[i]
	matchLines[(i-1)*2+2] <- gsub(pattern = "0", replacement = " ", temp)
	geneNames[(i-1)*2+1] <- AsanoDF$Gene[i]
	geneNames[(i-1)*2+2] <- AsanoDF$Gene[i]
	if(sum(matchList[["middle51"]][[i]]) > 0){
		hit[(i-1)*2+2] <- as.character(sum(matchList[["middle51"]][[i]]))
	}
}
compareDF[["middle51"]] <- data.frame(Gene = geneNames, Seq_Match = matchLines, hit = hit)

matchLines <- character(length = length(matchList[["last51"]]) *2)
geneNames <- character(length = length(matchList[["last51"]]) *2)
hit <- character(length = length(matchList[["last51"]]) *2)
for(i in 1:length(matchList[["last51"]])){
	if(i%%100 == 0) cat(i, ", ", sep = "")
	if(i == length(matchList[["last51"]])) cat("\n")
	temp <- paste(matchList[["last51"]][[i]], collapse = "")
	matchLines[(i-1)*2+1] <- AsanoDF$CDS_last51bp[i]
	matchLines[(i-1)*2+2] <- gsub(pattern = "0", replacement = " ", temp)
	geneNames[(i-1)*2+1] <- AsanoDF$Gene[i]
	geneNames[(i-1)*2+2] <- AsanoDF$Gene[i]
	if(sum(matchList[["last51"]][[i]]) > 0){
		hit[(i-1)*2+2] <- as.character(sum(matchList[["last51"]][[i]]))
	}
}
compareDF[["last51"]] <- data.frame(Gene = geneNames, Seq_Match = matchLines, hit = hit)

matchLines <- character(length = length(matchList[["UTR5"]]) *2)
geneNames <- character(length = length(matchList[["UTR5"]]) *2)
hit <- character(length = length(matchList[["UTR5"]]) *2)
for(i in 1:length(matchList[["UTR5"]])){
	if(i%%100 == 0) cat(i, ", ", sep = "")
	if(i == length(matchList[["UTR5"]])) cat("\n")
	temp <- paste(matchList[["UTR5"]][[i]], collapse = "")
	matchLines[(i-1)*2+1] <- AsanoDF$UTR5[i]
	matchLines[(i-1)*2+2] <- gsub(pattern = "0", replacement = " ", temp)
	geneNames[(i-1)*2+1] <- AsanoDF$Gene[i]
	geneNames[(i-1)*2+2] <- AsanoDF$Gene[i]
	if(sum(matchList[["UTR5"]][[i]]) > 0){
		hit[(i-1)*2+2] <- as.character(sum(matchList[["UTR5"]][[i]]))
	}
}
compareDF[["UTR5"]] <- data.frame(Gene = geneNames, Seq_Match = matchLines, hit = hit)

matchLines <- character(length = length(matchList[["UTR3"]]) *2)
geneNames <- character(length = length(matchList[["UTR3"]]) *2)
hit <- character(length = length(matchList[["UTR3"]]) *2)
for(i in 1:length(matchList[["UTR3"]])){
	if(i%%100 == 0) cat(i, ", ", sep = "")
	if(i == length(matchList[["UTR3"]])) cat("\n")
	temp <- paste(matchList[["UTR3"]][[i]], collapse = "")
	matchLines[(i-1)*2+1] <- AsanoDF$UTR3[i]
	matchLines[(i-1)*2+2] <- gsub(pattern = "0", replacement = " ", temp)
	geneNames[(i-1)*2+1] <- AsanoDF$Gene[i]
	geneNames[(i-1)*2+2] <- AsanoDF$Gene[i]
	if(sum(matchList[["UTR3"]][[i]]) > 0){
		hit[(i-1)*2+2] <- as.character(sum(matchList[["UTR3"]][[i]]))
	}
}
compareDF[["UTR3"]] <- data.frame(Gene = geneNames, Seq_Match = matchLines, hit = hit)

matchLines <- character(length = length(matchList[["CDS_total"]]) *2)
geneNames <- character(length = length(matchList[["CDS_total"]]) *2)
hit <- character(length = length(matchList[["CDS_total"]]) *2)
for(i in 1:length(matchList[["CDS_total"]])){
	if(i%%100 == 0) cat(i, ", ", sep = "")
	if(i == length(matchList[["CDS_total"]])) cat("\n")
	temp <- paste(matchList[["CDS_total"]][[i]], collapse = "")
	matchLines[(i-1)*2+1] <- AsanoDF$CDS_total[i]
	matchLines[(i-1)*2+2] <- gsub(pattern = "0", replacement = " ", temp)
	geneNames[(i-1)*2+1] <- AsanoDF$Gene[i]
	geneNames[(i-1)*2+2] <- AsanoDF$Gene[i]
	if(sum(matchList[["CDS_total"]][[i]]) > 0){
		hit[(i-1)*2+2] <- as.character(sum(matchList[["CDS_total"]][[i]]))
	}
}
compareDF[["CDS_total"]] <- data.frame(Gene = geneNames, Seq_Match = matchLines, hit = hit)

matchLines <- character(length = length(matchList[["CDS_1toLast52"]]) *2)
geneNames <- character(length = length(matchList[["CDS_1toLast52"]]) *2)
hit <- character(length = length(matchList[["CDS_1toLast52"]]) *2)
for(i in 1:length(matchList[["CDS_1toLast52"]])){
	if(i%%100 == 0) cat(i, ", ", sep = "")
	if(i == length(matchList[["CDS_1toLast52"]])) cat("\n")
	temp <- paste(matchList[["CDS_1toLast52"]][[i]], collapse = "")
	matchLines[(i-1)*2+1] <- AsanoDF$CDS_1toLast52[i]
	matchLines[(i-1)*2+2] <- gsub(pattern = "0", replacement = " ", temp)
	geneNames[(i-1)*2+1] <- AsanoDF$Gene[i]
	geneNames[(i-1)*2+2] <- AsanoDF$Gene[i]
	if(sum(matchList[["CDS_1toLast52"]][[i]]) > 0){
		hit[(i-1)*2+2] <- as.character(sum(matchList[["CDS_1toLast52"]][[i]]))
	}
}
compareDF[["CDS_1toLast52"]] <- data.frame(Gene = geneNames, Seq_Match = matchLines, hit = hit)

matchLines <- character(length = length(matchList[["CDS_52toLast"]]) *2)
geneNames <- character(length = length(matchList[["CDS_52toLast"]]) *2)
hit <- character(length = length(matchList[["CDS_52toLast"]]) *2)
for(i in 1:length(matchList[["CDS_52toLast"]])){
	if(i%%100 == 0) cat(i, ", ", sep = "")
	if(i == length(matchList[["CDS_52toLast"]])) cat("\n")
	temp <- paste(matchList[["CDS_52toLast"]][[i]], collapse = "")
	matchLines[(i-1)*2+1] <- AsanoDF$CDS_52toLast[i]
	matchLines[(i-1)*2+2] <- gsub(pattern = "0", replacement = " ", temp)
	geneNames[(i-1)*2+1] <- AsanoDF$Gene[i]
	geneNames[(i-1)*2+2] <- AsanoDF$Gene[i]
	if(sum(matchList[["CDS_52toLast"]][[i]]) > 0){
		hit[(i-1)*2+2] <- as.character(sum(matchList[["CDS_52toLast"]][[i]]))
	}
}
compareDF[["CDS_52toLast"]] <- data.frame(Gene = geneNames, Seq_Match = matchLines, hit = hit)

setwd(HOMEDIR)
setwd("data")
setwd("motifHits")
setwd("XLSX")
write.xlsx(compareDF, file = "20240619_AGYYGRY_6thTrial_4940.xlsx")

setwd("../RData")
save(compareDF, file = "20240619_AGYYGRY_6thTrial_4940.RData")



matchDFlist <- list()

matchDF <- AsanoDF
for(j in 1:length(PATTERNS0)){
	PATTERN <- PATTERNS0[j]
	cat(j, ", ", PATTERN, "\n", sep = "")
	matchDF[[PATTERN]] <- as.integer(NA)
	for(i in 1:nrow(AsanoDF)){
		if(i%%100 == 0) cat(i, ", ", sep = "")
		if(i == nrow(AsanoDF)) cat("\n")
		SEQ <- AsanoDF$CDS_first51bp[i]
		# 配列がある時のみモチーフを探す
		if(is.na(SEQ) == FALSE){
			matchDF[[PATTERN]][i] <- 0
			temp <- matchPattern(pattern = PATTERN, subject = DNAString(SEQ), 
					fixed = FALSE, max.mismatch = 0)
			if(length(temp@ranges@start) > 0){
				matchDF[[PATTERN]][i] <- length(temp@ranges@start) / nchar(SEQ) * 51
			}
		}
	}
}

for(j in 1:length(PATTERNS1)){
	PATTERN <- PATTERNS1[j]
	cat(j, ", ", PATTERN, "\n", sep = "")
	colName <- paste(PATTERN, "m1", sep = "")
	matchDF[[colName]] <- as.integer(NA)
	for(i in 1:nrow(AsanoDF)){
		if(i%%100 == 0) cat(i, ", ", sep = "")
		if(i == nrow(AsanoDF)) cat("\n")
		SEQ <- AsanoDF$CDS_first51bp[i]
		# 配列がある時のみモチーフを探す
		if(is.na(SEQ) == FALSE){
			matchDF[[colName]][i] <- 0
			temp <- matchPattern(pattern = PATTERN, subject = DNAString(SEQ), 
					fixed = FALSE, max.mismatch = 1)
			if(length(temp@ranges@start) > 0){
				matchDF[[colName]][i] <- length(temp@ranges@start) / nchar(SEQ) * 51
			}
		}
	}
}


matchDF$hit <- as.integer(0)

index <- seq(2, nrow(compareDF[["first51"]]), 2)

matchDF$hit <- as.integer(compareDF[["first51"]]$hit[index]) /  nchar(compareDF[["first51"]]$Seq_Match[index]) * 51
matchDF$hit[is.na(matchDF$hit) == TRUE] <- 0
matchDF$hit[is.na(matchDF$CDS_first51bp) == TRUE] <- NA

colnames(matchDF)[ncol(matchDF)] <- "hit_first51"
matchDF[, colnames(matchDF) == "CDS_first51bp"] <- NULL
matchDF[, colnames(matchDF) == "CDS_middle51bp"] <- NULL
matchDF[, colnames(matchDF) == "CDS_last51bp"] <- NULL

matchDF[, colnames(matchDF) == "UTR5"] <- NULL
matchDF[, colnames(matchDF) == "UTR3"] <- NULL
matchDF[, colnames(matchDF) == "CDS_total"] <- NULL
matchDF[, colnames(matchDF) == "CDS_1toLast52"] <- NULL
matchDF[, colnames(matchDF) == "CDS_52toLast"] <- NULL

matchDFlist[["first51"]] <- matchDF


matchDF <- AsanoDF
for(j in 1:length(PATTERNS0)){
	PATTERN <- PATTERNS0[j]
	cat(j, ", ", PATTERN, "\n", sep = "")
	matchDF[[PATTERN]] <- as.integer(NA)
	for(i in 1:nrow(AsanoDF)){
		if(i%%100 == 0) cat(i, ", ", sep = "")
		if(i == nrow(AsanoDF)) cat("\n")
		SEQ <- AsanoDF$CDS_middle51bp[i]
		# 配列がある時のみモチーフを探す
		if(is.na(SEQ) == FALSE){
			matchDF[[PATTERN]][i] <- 0
			temp <- matchPattern(pattern = PATTERN, subject = DNAString(SEQ), 
					fixed = FALSE, max.mismatch = 0)
			if(length(temp@ranges@start) > 0){
				matchDF[[PATTERN]][i] <- length(temp@ranges@start) / nchar(SEQ) * 51
			}
		}
	}
}

for(j in 1:length(PATTERNS1)){
	PATTERN <- PATTERNS1[j]
	cat(j, ", ", PATTERN, "\n", sep = "")
	colName <- paste(PATTERN, "m1", sep = "")
	matchDF[[colName]] <- as.integer(NA)
	for(i in 1:nrow(AsanoDF)){
		if(i%%100 == 0) cat(i, ", ", sep = "")
		if(i == nrow(AsanoDF)) cat("\n")
		SEQ <- AsanoDF$CDS_middle51bp[i]
		# 配列がある時のみモチーフを探す
		if(is.na(SEQ) == FALSE){
			matchDF[[colName]][i] <- 0
			temp <- matchPattern(pattern = PATTERN, subject = DNAString(SEQ), 
					fixed = FALSE, max.mismatch = 1)
			if(length(temp@ranges@start) > 0){
				matchDF[[colName]][i] <- length(temp@ranges@start) / nchar(SEQ) * 51
			}
		}
	}
}

matchDF$hit <- as.integer(0)

index <- seq(2, nrow(compareDF[["middle51"]]), 2)

matchDF$hit <- as.integer(compareDF[["middle51"]]$hit[index]) /  nchar(compareDF[["middle51"]]$Seq_Match[index]) * 51
matchDF$hit[is.na(matchDF$hit) == TRUE] <- 0
matchDF$hit[is.na(matchDF$CDS_middle51bp) == TRUE] <- NA

colnames(matchDF)[ncol(matchDF)] <- "hit_middle51"
matchDF[, colnames(matchDF) == "CDS_first51bp"] <- NULL
matchDF[, colnames(matchDF) == "CDS_middle51bp"] <- NULL
matchDF[, colnames(matchDF) == "CDS_last51bp"] <- NULL

matchDF[, colnames(matchDF) == "UTR5"] <- NULL
matchDF[, colnames(matchDF) == "UTR3"] <- NULL
matchDF[, colnames(matchDF) == "CDS_total"] <- NULL
matchDF[, colnames(matchDF) == "CDS_1toLast52"] <- NULL
matchDF[, colnames(matchDF) == "CDS_52toLast"] <- NULL

matchDFlist[["middle51"]] <- matchDF


matchDF <- AsanoDF
for(j in 1:length(PATTERNS0)){
	PATTERN <- PATTERNS0[j]
	cat(j, ", ", PATTERN, "\n", sep = "")
	matchDF[[PATTERN]] <- as.integer(NA)
	for(i in 1:nrow(AsanoDF)){
		if(i%%100 == 0) cat(i, ", ", sep = "")
		if(i == nrow(AsanoDF)) cat("\n")
		SEQ <- AsanoDF$CDS_last51bp[i]
		# 配列がある時のみモチーフを探す
		if(is.na(SEQ) == FALSE){
			matchDF[[PATTERN]][i] <- 0
			temp <- matchPattern(pattern = PATTERN, subject = DNAString(SEQ), 
					fixed = FALSE, max.mismatch = 0)
			if(length(temp@ranges@start) > 0){
				matchDF[[PATTERN]][i] <- length(temp@ranges@start) / nchar(SEQ) * 51
			}
		}
	}
}

for(j in 1:length(PATTERNS1)){
	PATTERN <- PATTERNS1[j]
	cat(j, ", ", PATTERN, "\n", sep = "")
	colName <- paste(PATTERN, "m1", sep = "")
	matchDF[[colName]] <- as.integer(NA)
	for(i in 1:nrow(AsanoDF)){
		if(i%%100 == 0) cat(i, ", ", sep = "")
		if(i == nrow(AsanoDF)) cat("\n")
		SEQ <- AsanoDF$CDS_last51bp[i]
		# 配列がある時のみモチーフを探す
		if(is.na(SEQ) == FALSE){
			matchDF[[colName]][i] <- 0
			temp <- matchPattern(pattern = PATTERN, subject = DNAString(SEQ), 
					fixed = FALSE, max.mismatch = 1)
			if(length(temp@ranges@start) > 0){
				matchDF[[colName]][i] <- length(temp@ranges@start) / nchar(SEQ) * 51
			}
		}
	}
}

matchDF$hit <- as.integer(0)

index <- seq(2, nrow(compareDF[["last51"]]), 2)

matchDF$hit <- as.integer(compareDF[["last51"]]$hit[index]) /  nchar(compareDF[["last51"]]$Seq_Match[index]) * 51
matchDF$hit[is.na(matchDF$hit) == TRUE] <- 0
matchDF$hit[is.na(matchDF$CDS_last51bp) == TRUE] <- NA

colnames(matchDF)[ncol(matchDF)] <- "hit_last51"
matchDF[, colnames(matchDF) == "CDS_first51bp"] <- NULL
matchDF[, colnames(matchDF) == "CDS_middle51bp"] <- NULL
matchDF[, colnames(matchDF) == "CDS_last51bp"] <- NULL

matchDF[, colnames(matchDF) == "UTR5"] <- NULL
matchDF[, colnames(matchDF) == "UTR3"] <- NULL
matchDF[, colnames(matchDF) == "CDS_total"] <- NULL
matchDF[, colnames(matchDF) == "CDS_1toLast52"] <- NULL
matchDF[, colnames(matchDF) == "CDS_52toLast"] <- NULL

matchDFlist[["last51"]] <- matchDF



matchDF <- AsanoDF
for(j in 1:length(PATTERNS0)){
	PATTERN <- PATTERNS0[j]
	cat(j, ", ", PATTERN, "\n", sep = "")
	matchDF[[PATTERN]] <- as.integer(NA)
	for(i in 1:nrow(AsanoDF)){
		if(i%%100 == 0) cat(i, ", ", sep = "")
		if(i == nrow(AsanoDF)) cat("\n")
		SEQ <- AsanoDF$UTR5[i]
		# 配列がある時のみモチーフを探す
		if(is.na(SEQ) == FALSE){
			matchDF[[PATTERN]][i] <- 0
			temp <- matchPattern(pattern = PATTERN, subject = DNAString(SEQ), 
					fixed = FALSE, max.mismatch = 0)
			if(length(temp@ranges@start) > 0){
				matchDF[[PATTERN]][i] <- length(temp@ranges@start) / nchar(SEQ) * 51
			}
		}
	}
}

for(j in 1:length(PATTERNS1)){
	PATTERN <- PATTERNS1[j]
	cat(j, ", ", PATTERN, "\n", sep = "")
	colName <- paste(PATTERN, "m1", sep = "")
	matchDF[[colName]] <- as.integer(NA)
	for(i in 1:nrow(AsanoDF)){
		if(i%%100 == 0) cat(i, ", ", sep = "")
		if(i == nrow(AsanoDF)) cat("\n")
		SEQ <- AsanoDF$UTR5[i]
		# 配列がある時のみモチーフを探す
		if(is.na(SEQ) == FALSE){
			matchDF[[colName]][i] <- 0
			temp <- matchPattern(pattern = PATTERN, subject = DNAString(SEQ), 
					fixed = FALSE, max.mismatch = 1)
			if(length(temp@ranges@start) > 0){
				matchDF[[colName]][i] <- length(temp@ranges@start) / nchar(SEQ) * 51
			}
		}
	}
}

matchDF$hit <- as.integer(0)

index <- seq(2, nrow(compareDF[["UTR5"]]), 2)

matchDF$hit <- as.integer(compareDF[["UTR5"]]$hit[index]) /  nchar(compareDF[["UTR5"]]$Seq_Match[index]) * 51
matchDF$hit[is.na(matchDF$hit) == TRUE] <- 0
matchDF$hit[is.na(matchDF$UTR5) == TRUE] <- NA

colnames(matchDF)[ncol(matchDF)] <- "hit_UTR5"
matchDF[, colnames(matchDF) == "CDS_first51bp"] <- NULL
matchDF[, colnames(matchDF) == "CDS_middle51bp"] <- NULL
matchDF[, colnames(matchDF) == "CDS_last51bp"] <- NULL

matchDF[, colnames(matchDF) == "UTR5"] <- NULL
matchDF[, colnames(matchDF) == "UTR3"] <- NULL
matchDF[, colnames(matchDF) == "CDS_total"] <- NULL
matchDF[, colnames(matchDF) == "CDS_1toLast52"] <- NULL
matchDF[, colnames(matchDF) == "CDS_52toLast"] <- NULL

matchDFlist[["UTR5"]] <- matchDF



matchDF <- AsanoDF
for(j in 1:length(PATTERNS0)){
	PATTERN <- PATTERNS0[j]
	cat(j, ", ", PATTERN, "\n", sep = "")
	matchDF[[PATTERN]] <- as.integer(NA)
	for(i in 1:nrow(AsanoDF)){
		if(i%%100 == 0) cat(i, ", ", sep = "")
		if(i == nrow(AsanoDF)) cat("\n")
		SEQ <- AsanoDF$UTR3[i]
		# 配列がある時のみモチーフを探す
		if(is.na(SEQ) == FALSE){
			matchDF[[PATTERN]][i] <- 0
			temp <- matchPattern(pattern = PATTERN, subject = DNAString(SEQ), 
					fixed = FALSE, max.mismatch = 0)
			if(length(temp@ranges@start) > 0){
				matchDF[[PATTERN]][i] <- length(temp@ranges@start) / nchar(SEQ) * 51
			}
		}
	}
}

for(j in 1:length(PATTERNS1)){
	PATTERN <- PATTERNS1[j]
	cat(j, ", ", PATTERN, "\n", sep = "")
	colName <- paste(PATTERN, "m1", sep = "")
	matchDF[[colName]] <- as.integer(NA)
	for(i in 1:nrow(AsanoDF)){
		if(i%%100 == 0) cat(i, ", ", sep = "")
		if(i == nrow(AsanoDF)) cat("\n")
		SEQ <- AsanoDF$UTR3[i]
		# 配列がある時のみモチーフを探す
		if(is.na(SEQ) == FALSE){
			matchDF[[colName]][i] <- 0
			temp <- matchPattern(pattern = PATTERN, subject = DNAString(SEQ), 
					fixed = FALSE, max.mismatch = 1)
			if(length(temp@ranges@start) > 0){
				matchDF[[colName]][i] <- length(temp@ranges@start) / nchar(SEQ) * 51
			}
		}
	}
}

matchDF$hit <- as.integer(0)

index <- seq(2, nrow(compareDF[["UTR3"]]), 2)

matchDF$hit <- as.integer(compareDF[["UTR3"]]$hit[index]) /  nchar(compareDF[["UTR3"]]$Seq_Match[index]) * 51
matchDF$hit[is.na(matchDF$hit) == TRUE] <- 0
matchDF$hit[is.na(matchDF$UTR3) == TRUE] <- NA

colnames(matchDF)[ncol(matchDF)] <- "hit_UTR3"
matchDF[, colnames(matchDF) == "CDS_first51bp"] <- NULL
matchDF[, colnames(matchDF) == "CDS_middle51bp"] <- NULL
matchDF[, colnames(matchDF) == "CDS_last51bp"] <- NULL

matchDF[, colnames(matchDF) == "UTR5"] <- NULL
matchDF[, colnames(matchDF) == "UTR3"] <- NULL
matchDF[, colnames(matchDF) == "CDS_total"] <- NULL
matchDF[, colnames(matchDF) == "CDS_1toLast52"] <- NULL
matchDF[, colnames(matchDF) == "CDS_52toLast"] <- NULL

matchDFlist[["UTR3"]] <- matchDF



matchDF <- AsanoDF
for(j in 1:length(PATTERNS0)){
	PATTERN <- PATTERNS0[j]
	cat(j, ", ", PATTERN, "\n", sep = "")
	matchDF[[PATTERN]] <- as.integer(NA)
	for(i in 1:nrow(AsanoDF)){
		if(i%%100 == 0) cat(i, ", ", sep = "")
		if(i == nrow(AsanoDF)) cat("\n")
		SEQ <- AsanoDF$CDS_total[i]
		# 配列がある時のみモチーフを探す
		if(is.na(SEQ) == FALSE){
			matchDF[[PATTERN]][i] <- 0
			temp <- matchPattern(pattern = PATTERN, subject = DNAString(SEQ), 
					fixed = FALSE, max.mismatch = 0)
			if(length(temp@ranges@start) > 0){
				matchDF[[PATTERN]][i] <- length(temp@ranges@start) / nchar(SEQ) * 51
			}
		}
	}
}

for(j in 1:length(PATTERNS1)){
	PATTERN <- PATTERNS1[j]
	cat(j, ", ", PATTERN, "\n", sep = "")
	colName <- paste(PATTERN, "m1", sep = "")
	matchDF[[colName]] <- as.integer(NA)
	for(i in 1:nrow(AsanoDF)){
		if(i%%100 == 0) cat(i, ", ", sep = "")
		if(i == nrow(AsanoDF)) cat("\n")
		SEQ <- AsanoDF$CDS_total[i]
		# 配列がある時のみモチーフを探す
		if(is.na(SEQ) == FALSE){
			matchDF[[colName]][i] <- 0
			temp <- matchPattern(pattern = PATTERN, subject = DNAString(SEQ), 
					fixed = FALSE, max.mismatch = 1)
			if(length(temp@ranges@start) > 0){
				matchDF[[colName]][i] <- length(temp@ranges@start) / nchar(SEQ) * 51
			}
		}
	}
}

matchDF$hit <- as.integer(0)

index <- seq(2, nrow(compareDF[["CDS_total"]]), 2)

matchDF$hit <- as.integer(compareDF[["CDS_total"]]$hit[index]) /  nchar(compareDF[["CDS_total"]]$Seq_Match[index]) * 51
matchDF$hit[is.na(matchDF$hit) == TRUE] <- 0
matchDF$hit[is.na(matchDF$CDS_total) == TRUE] <- NA

colnames(matchDF)[ncol(matchDF)] <- "hit_CDS_total"
matchDF[, colnames(matchDF) == "CDS_first51bp"] <- NULL
matchDF[, colnames(matchDF) == "CDS_middle51bp"] <- NULL
matchDF[, colnames(matchDF) == "CDS_last51bp"] <- NULL

matchDF[, colnames(matchDF) == "UTR5"] <- NULL
matchDF[, colnames(matchDF) == "UTR3"] <- NULL
matchDF[, colnames(matchDF) == "CDS_total"] <- NULL
matchDF[, colnames(matchDF) == "CDS_1toLast52"] <- NULL
matchDF[, colnames(matchDF) == "CDS_52toLast"] <- NULL

matchDFlist[["CDS_total"]] <- matchDF



matchDF <- AsanoDF
for(j in 1:length(PATTERNS0)){
	PATTERN <- PATTERNS0[j]
	cat(j, ", ", PATTERN, "\n", sep = "")
	matchDF[[PATTERN]] <- as.integer(NA)
	for(i in 1:nrow(AsanoDF)){
		if(i%%100 == 0) cat(i, ", ", sep = "")
		if(i == nrow(AsanoDF)) cat("\n")
		SEQ <- AsanoDF$CDS_1toLast52[i]
		# 配列がある時のみモチーフを探す
		if(is.na(SEQ) == FALSE){
			matchDF[[PATTERN]][i] <- 0
			temp <- matchPattern(pattern = PATTERN, subject = DNAString(SEQ), 
					fixed = FALSE, max.mismatch = 0)
			if(length(temp@ranges@start) > 0){
				matchDF[[PATTERN]][i] <- length(temp@ranges@start) / nchar(SEQ) * 51
			}
		}
	}
}

for(j in 1:length(PATTERNS1)){
	PATTERN <- PATTERNS1[j]
	cat(j, ", ", PATTERN, "\n", sep = "")
	colName <- paste(PATTERN, "m1", sep = "")
	matchDF[[colName]] <- as.integer(NA)
	for(i in 1:nrow(AsanoDF)){
		if(i%%100 == 0) cat(i, ", ", sep = "")
		if(i == nrow(AsanoDF)) cat("\n")
		SEQ <- AsanoDF$CDS_1toLast52[i]
		# 配列がある時のみモチーフを探す
		if(is.na(SEQ) == FALSE){
			matchDF[[colName]][i] <- 0
			temp <- matchPattern(pattern = PATTERN, subject = DNAString(SEQ), 
					fixed = FALSE, max.mismatch = 1)
			if(length(temp@ranges@start) > 0){
				matchDF[[colName]][i] <- length(temp@ranges@start) / nchar(SEQ) * 51
			}
		}
	}
}

matchDF$hit <- as.integer(0)

index <- seq(2, nrow(compareDF[["CDS_1toLast52"]]), 2)

matchDF$hit <- as.integer(compareDF[["CDS_1toLast52"]]$hit[index]) /  nchar(compareDF[["CDS_1toLast52"]]$Seq_Match[index]) * 51
matchDF$hit[is.na(matchDF$hit) == TRUE] <- 0
matchDF$hit[is.na(matchDF$CDS_1toLast52) == TRUE] <- NA

colnames(matchDF)[ncol(matchDF)] <- "hit_CDS_1toLast52"
matchDF[, colnames(matchDF) == "CDS_first51bp"] <- NULL
matchDF[, colnames(matchDF) == "CDS_middle51bp"] <- NULL
matchDF[, colnames(matchDF) == "CDS_last51bp"] <- NULL

matchDF[, colnames(matchDF) == "UTR5"] <- NULL
matchDF[, colnames(matchDF) == "UTR3"] <- NULL
matchDF[, colnames(matchDF) == "CDS_total"] <- NULL
matchDF[, colnames(matchDF) == "CDS_1toLast52"] <- NULL
matchDF[, colnames(matchDF) == "CDS_52toLast"] <- NULL

matchDFlist[["CDS_1toLast52"]] <- matchDF



matchDF <- AsanoDF
for(j in 1:length(PATTERNS0)){
	PATTERN <- PATTERNS0[j]
	cat(j, ", ", PATTERN, "\n", sep = "")
	matchDF[[PATTERN]] <- as.integer(NA)
	for(i in 1:nrow(AsanoDF)){
		if(i%%100 == 0) cat(i, ", ", sep = "")
		if(i == nrow(AsanoDF)) cat("\n")
		SEQ <- AsanoDF$CDS_52toLast[i]
		# 配列がある時のみモチーフを探す
		if(is.na(SEQ) == FALSE){
			matchDF[[PATTERN]][i] <- 0
			temp <- matchPattern(pattern = PATTERN, subject = DNAString(SEQ), 
					fixed = FALSE, max.mismatch = 0)
			if(length(temp@ranges@start) > 0){
				matchDF[[PATTERN]][i] <- length(temp@ranges@start) / nchar(SEQ) * 51
			}
		}
	}
}

for(j in 1:length(PATTERNS1)){
	PATTERN <- PATTERNS1[j]
	cat(j, ", ", PATTERN, "\n", sep = "")
	colName <- paste(PATTERN, "m1", sep = "")
	matchDF[[colName]] <- as.integer(NA)
	for(i in 1:nrow(AsanoDF)){
		if(i%%100 == 0) cat(i, ", ", sep = "")
		if(i == nrow(AsanoDF)) cat("\n")
		SEQ <- AsanoDF$CDS_52toLast[i]
		# 配列がある時のみモチーフを探す
		if(is.na(SEQ) == FALSE){
			matchDF[[colName]][i] <- 0
			temp <- matchPattern(pattern = PATTERN, subject = DNAString(SEQ), 
					fixed = FALSE, max.mismatch = 1)
			if(length(temp@ranges@start) > 0){
				matchDF[[colName]][i] <- length(temp@ranges@start) / nchar(SEQ) * 51
			}
		}
	}
}

matchDF$hit <- as.integer(0)

index <- seq(2, nrow(compareDF[["CDS_52toLast"]]), 2)

matchDF$hit <- as.integer(compareDF[["CDS_52toLast"]]$hit[index]) /  nchar(compareDF[["CDS_52toLast"]]$Seq_Match[index]) * 51
matchDF$hit[is.na(matchDF$hit) == TRUE] <- 0
matchDF$hit[is.na(matchDF$CDS_52toLast) == TRUE] <- NA

colnames(matchDF)[ncol(matchDF)] <- "hit_CDS_52toLast"
matchDF[, colnames(matchDF) == "CDS_first51bp"] <- NULL
matchDF[, colnames(matchDF) == "CDS_middle51bp"] <- NULL
matchDF[, colnames(matchDF) == "CDS_last51bp"] <- NULL

matchDF[, colnames(matchDF) == "UTR5"] <- NULL
matchDF[, colnames(matchDF) == "UTR3"] <- NULL
matchDF[, colnames(matchDF) == "CDS_total"] <- NULL
matchDF[, colnames(matchDF) == "CDS_1toLast52"] <- NULL
matchDF[, colnames(matchDF) == "CDS_52toLast"] <- NULL

matchDFlist[["CDS_52toLast"]] <- matchDF


setwd(HOMEDIR)
setwd("data")
setwd("motifHits")
setwd("XLSX")
write.xlsx(matchDFlist, file = "20240619_matchDFlist_6thTrial_4940.xlsx")

setwd("../RData")
save(matchDFlist, file = "20240619_matchDFlist_6thTrial_4940.RData")



## original note: 20240619_Asano_motif_2.txt

setwd(HOMEDIR)
setwd("data")
setwd("motifHits")
setwd("RData")
load(file = "20240619_NewDF.RData")
load(file = "20240619_matchDFlist_6thTrial_4940.RData")

NewDF2 <- NewDF[, c(1, 6, 15, 20, 33:42)]
NewDF2$hit_first51 <- matchDFlist[["first51"]][,ncol(matchDFlist[["first51"]])]
NewDF2$hit_middle51 <- matchDFlist[["middle51"]][,ncol(matchDFlist[["middle51"]])]
NewDF2$hit_last51 <- matchDFlist[["last51"]][,ncol(matchDFlist[["last51"]])]
NewDF2$hit_UTR5 <- matchDFlist[["UTR5"]][,ncol(matchDFlist[["UTR5"]])]
NewDF2$hit_UTR3 <- matchDFlist[["UTR3"]][,ncol(matchDFlist[["UTR3"]])]
NewDF2$hit_CDS_total <- matchDFlist[["CDS_total"]][,ncol(matchDFlist[["CDS_total"]])]
NewDF2$hit_CDS_1toLast52 <- matchDFlist[["CDS_1toLast52"]][,ncol(matchDFlist[["CDS_1toLast52"]])]
NewDF2$hit_CDS_52toLast <- matchDFlist[["CDS_52toLast"]][,ncol(matchDFlist[["CDS_52toLast"]])]

setwd(HOMEDIR)
setwd("data")
setwd("motifHits")
setwd("RData")
save(NewDF2, file = "20240619_NewDF2.RData")
setwd("../XLSX")
write.xlsx(NewDF2, file = "20240619_NewDF2.xlsx")



## original note: 20240717_Asano_motif_1.txt

setwd(HOMEDIR)
setwd("data")
setwd("motifHits")
setwd("RData")
load(file = "20240619_NewDF2.RData")	# NewDF2

setwd("../TXT")
Top64 <- readLines(con = "tif34 down top 64.txt")

NewDF2$Top64 <- FALSE
for(i in 1:64){
	NewDF2$Top64[NewDF2$Gene == Top64[i]] <- TRUE
}

DF_hit_first51 <- subset(NewDF2, is.na(hit_first51) == FALSE)
RPS_hit_first51 <- subset(DF_hit_first51, RPS == TRUE)$hit_first51
RPL_hit_first51 <- subset(DF_hit_first51, RPL == TRUE)$hit_first51
top50_hit_first51 <- subset(DF_hit_first51, tif34down_top50 == TRUE)$hit_first51
bottom50_hit_first51 <- subset(DF_hit_first51, tif34down_bottom50 == TRUE)$hit_first51
S4_246_hit_first51 <- subset(DF_hit_first51, S4_246 == TRUE)$hit_first51
All_hit_first51 <- DF_hit_first51$hit_first51

HSR64_hit_first51 <- subset(DF_hit_first51, Top64 == TRUE)$hit_first51

targets <- c("hit_first51", "hit_middle51", "hit_last51", "hit_UTR5", "hit_UTR3", "hit_CDS_total", "hit_CDS_1toLast52", "hit_CDS_52toLast")
pvaluesDF <- data.frame(target = targets, All_vs_RPS = as.numeric(NA), All_vs_RPL = as.numeric(NA), 
			All_vs_top50 = as.numeric(NA), All_vs_bottom50 = as.numeric(NA), All_vs_S4_246 = as.numeric(NA), 
			All_vs_HSR64 = as.numeric(NA))

meanvaluesDF <- data.frame(target = targets, All = as.numeric(NA), RPS = as.numeric(NA), RPL = as.numeric(NA), 
			top50 = as.numeric(NA), bottom50 = as.numeric(NA), S4_246 = as.numeric(NA), HSR64 = as.numeric(NA))

sdvaluesDF <- data.frame(target = targets, All = as.numeric(NA), RPS = as.numeric(NA), RPL = as.numeric(NA), 
			top50 = as.numeric(NA), bottom50 = as.numeric(NA), S4_246 = as.numeric(NA), HSR64 = as.numeric(NA))


for(t in 1:length(targets)){
	target <- targets[t]

	pvalues <- numeric(length = 6)

	# target <- "hit_first51"
	DF_hit <- subset(NewDF2, is.na(NewDF2[[target]]) == FALSE)
	RPS_hit <- subset(DF_hit, RPS == TRUE)[[target]]
	RPL_hit <- subset(DF_hit, RPL == TRUE)[[target]]
	top50_hit <- subset(DF_hit, tif34down_top50 == TRUE)[[target]]
	bottom50_hit <- subset(DF_hit, tif34down_bottom50 == TRUE)[[target]]
	S4_246_hit <- subset(DF_hit, S4_246 == TRUE)[[target]]
	All_hit <- DF_hit[[target]]
	HSR64_hit <- subset(DF_hit, Top64 == TRUE)[[target]]

	meanvaluesDF$All[t] <- mean(All_hit)
	meanvaluesDF$RPS[t] <- mean(RPS_hit)
	meanvaluesDF$RPL[t] <- mean(RPL_hit)
	meanvaluesDF$top50[t] <- mean(top50_hit)
	meanvaluesDF$bottom50[t] <- mean(bottom50_hit)
	meanvaluesDF$S4_246[t] <- mean(S4_246_hit)
	meanvaluesDF$HSR64[t] <- mean(HSR64_hit)

	sdvaluesDF$All[t] <- sd(All_hit)
	sdvaluesDF$RPS[t] <- sd(RPS_hit)
	sdvaluesDF$RPL[t] <- sd(RPL_hit)
	sdvaluesDF$top50[t] <- sd(top50_hit)
	sdvaluesDF$bottom50[t] <- sd(bottom50_hit)
	sdvaluesDF$S4_246[t] <- sd(S4_246_hit)
	sdvaluesDF$HSR64[t] <- sd(HSR64_hit)
	
	cat("## ", target, "\n", sep = "")
	cat("#  target==TRUE:         ", nrow(DF_hit), "\n", sep = "")
	cat("#  length(RPS_hit):      ", length(RPS_hit), "\n", sep = "")
	cat("#  length(RPL_hit):      ", length(RPL_hit), "\n", sep = "")
	cat("#  length(top50_hit):    ", length(top50_hit), "\n", sep = "")
	cat("#  length(bottom50_hit): ", length(bottom50_hit), "\n", sep = "")
	cat("#  length(S4_246_hit):   ", length(S4_246_hit), "\n", sep = "")
	cat("#  length(All_hit):      ", length(All_hit), "\n", sep = "")
	cat("#  length(HSR64_hit):      ", length(HSR64_hit), "\n", sep = "")


	cat("\n## t-test All_hit vs RPS_hit\n")
	temp <- t.test(x = All_hit, y = RPS_hit)
	pvalues[1] <- temp$p.value
	x <- paste("# ", temp)
	for(i in 1:length(x)){
		cat(x[i], "\n")
	}

	cat("\n## t-test All_hit vs RPL_hit\n")
	temp <- t.test(x = All_hit, y = RPL_hit)
	pvalues[2] <- temp$p.value
	x <- paste("# ", temp)
	for(i in 1:length(x)){
		cat(x[i], "\n")
	}

	cat("\n## t-test All_hit vs top50_hit\n")
	temp <- t.test(x = All_hit, y = top50_hit)
	pvalues[3] <- temp$p.value
	x <- paste("# ", temp)
	for(i in 1:length(x)){
		cat(x[i], "\n")
	}

	cat("\n## t-test All_hit vs bottom50_hit\n")
	temp <- t.test(x = All_hit, y = bottom50_hit)
	pvalues[4] <- temp$p.value
	x <- paste("# ", temp)
	for(i in 1:length(x)){
		cat(x[i], "\n")
	}

	cat("\n## t-test All_hit vs S4_246_hit\n")
	temp <- t.test(x = All_hit, y = S4_246_hit)
	pvalues[5] <- temp$p.value
	x <- paste("# ", temp)
	for(i in 1:length(x)){
		cat(x[i], "\n")
	}

	cat("\n## t-test All_hit vs HSR64_hit\n")
	temp <- t.test(x = All_hit, y = HSR64_hit)
	pvalues[6] <- temp$p.value
	x <- paste("# ", temp)
	for(i in 1:length(x)){
		cat(x[i], "\n")
	}
	
	pvaluesDF[t, 2:7] <- pvalues
	cat("\n\n")
}

pvaluesDFsig <- pvaluesDF
pvaluesDFsig[, 2:7][pvaluesDFsig[, 2:7] > 0.05] <- NA

outList <- list()
outList$pvalues <- pvaluesDF
outList$significant <- pvaluesDFsig
outList$mean <- meanvaluesDF
outList$sd <- sdvaluesDF


setwd(HOMEDIR)
setwd("data")
setwd("motifHits")
setwd("XLSX")
write.xlsx(outList, file = "20240717_motifHits_4940.xlsx")


