## original note: 20250420_motifHits_given_motif_mismatch.R.txt


args <- commandArgs(trailingOnly = TRUE)
MOTIF <- args[1]
MAX.MISMATCH <- as.integer(args[2])
DAY <- gsub(pattern = "-", replacement = "", x = as.character(Sys.Date()))

cat("args", args, "\n")
cat("MOTIF", MOTIF, "\n")
cat("MAX.MISMATCH", MAX.MISMATCH, "\n")
cat("DAY", DAY, "\n")

outDir <- paste(DAY, "_", MOTIF, "_m", MAX.MISMATCH, sep = "")

HOMEDIR <- ""

library(openxlsx)
library(Biostrings)


setwd(HOMEDIR)
setwd("data")
setwd("motifHits")
setwd("RData")
load(file = "20250420_AsanoDF.RData")

PATTERNS0 <- MOTIF

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
	
	if(is.na(SEQ) == FALSE){
		matchList[["first51"]][[GENE]] <- integer(length = nchar(SEQ))
	
		if(length(grep(pattern = "U", SEQ)) > 0) SEQ <- RNAString(SEQ)
	
		for(j in 1:length(PATTERNS0)){
			PATTERN <- PATTERNS0[j]
			temp <- matchPattern(pattern = PATTERN, subject = DNAString(SEQ), 
					fixed = FALSE, max.mismatch = MAX.MISMATCH)
			if(length(temp@ranges@start) > 0){
				starts <- temp@ranges@start[temp@ranges@start > 0]
				matchList[["first51"]][[GENE]][starts] <- 1
			}
		}
	}
	
	SEQ <- AsanoDF$CDS_middle51bp[i]
	if(is.na(SEQ) == TRUE){
		matchList[["middle51"]][[GENE]] <- integer(length = 0)
	}
		
	if(is.na(SEQ) == FALSE){
		matchList[["middle51"]][[GENE]] <- integer(length = nchar(SEQ))
	
		if(length(grep(pattern = "U", SEQ)) > 0) SEQ <- RNAString(SEQ)

		for(j in 1:length(PATTERNS0)){
			PATTERN <- PATTERNS0[j]
			temp <- matchPattern(pattern = PATTERN, subject = DNAString(SEQ), 
					fixed = FALSE, max.mismatch = MAX.MISMATCH)
			if(length(temp@ranges@start) > 0){
				starts <- temp@ranges@start[temp@ranges@start > 0]
				matchList[["middle51"]][[GENE]][starts] <- 1
			}
		}
	}

	SEQ <- AsanoDF$CDS_last51bp[i]
	if(is.na(SEQ) == TRUE){
		matchList[["last51"]][[GENE]] <- integer(length = 0)
	}
	
	if(is.na(SEQ) == FALSE){
		matchList[["last51"]][[GENE]] <- integer(length = nchar(SEQ))
	
		if(length(grep(pattern = "U", SEQ)) > 0) SEQ <- RNAString(SEQ)

		for(j in 1:length(PATTERNS0)){
			PATTERN <- PATTERNS0[j]
			temp <- matchPattern(pattern = PATTERN, subject = DNAString(SEQ), 
					fixed = FALSE, max.mismatch = MAX.MISMATCH)
			if(length(temp@ranges@start) > 0){
				starts <- temp@ranges@start[temp@ranges@start > 0]
				matchList[["last51"]][[GENE]][starts] <- 1
			}
		}
	}

	SEQ <- AsanoDF$UTR5[i]
	if(is.na(SEQ) == TRUE){
		matchList[["UTR5"]][[GENE]] <- integer(length = 0)
	}
	
	if(is.na(SEQ) == FALSE){
		matchList[["UTR5"]][[GENE]] <- integer(length = nchar(SEQ))
	
		if(length(grep(pattern = "U", SEQ)) > 0) SEQ <- RNAString(SEQ)

		for(j in 1:length(PATTERNS0)){
			PATTERN <- PATTERNS0[j]
			temp <- matchPattern(pattern = PATTERN, subject = DNAString(SEQ), 
					fixed = FALSE, max.mismatch = MAX.MISMATCH)
			if(length(temp@ranges@start) > 0){
				starts <- temp@ranges@start[temp@ranges@start > 0]
				matchList[["UTR5"]][[GENE]][starts] <- 1
			}
		}
	}

	SEQ <- AsanoDF$UTR3[i]
	if(is.na(SEQ) == TRUE){
		matchList[["UTR3"]][[GENE]] <- integer(length = 0)
	}
	
	if(is.na(SEQ) == FALSE){
		matchList[["UTR3"]][[GENE]] <- integer(length = nchar(SEQ))
	
		if(length(grep(pattern = "U", SEQ)) > 0) SEQ <- RNAString(SEQ)

		for(j in 1:length(PATTERNS0)){
			PATTERN <- PATTERNS0[j]
			temp <- matchPattern(pattern = PATTERN, subject = DNAString(SEQ), 
					fixed = FALSE, max.mismatch = MAX.MISMATCH)
			if(length(temp@ranges@start) > 0){
				starts <- temp@ranges@start[temp@ranges@start > 0]
				matchList[["UTR3"]][[GENE]][starts] <- 1
			}
		}
	}

	SEQ <- AsanoDF$CDS_total[i]
	if(is.na(SEQ) == TRUE){
		matchList[["CDS_total"]][[GENE]] <- integer(length = 0)
	}
	
	if(is.na(SEQ) == FALSE){
		matchList[["CDS_total"]][[GENE]] <- integer(length = nchar(SEQ))
	
		if(length(grep(pattern = "U", SEQ)) > 0) SEQ <- RNAString(SEQ)

		for(j in 1:length(PATTERNS0)){
			PATTERN <- PATTERNS0[j]
			temp <- matchPattern(pattern = PATTERN, subject = DNAString(SEQ), 
					fixed = FALSE, max.mismatch = MAX.MISMATCH)
			if(length(temp@ranges@start) > 0){
				starts <- temp@ranges@start[temp@ranges@start > 0]
				matchList[["CDS_total"]][[GENE]][starts] <- 1
			}
		}
	}

	SEQ <- AsanoDF$CDS_1toLast52[i]
	if(is.na(SEQ) == TRUE){
		matchList[["CDS_1toLast52"]][[GENE]] <- integer(length = 0)
	}
	
	if(is.na(SEQ) == FALSE){
		matchList[["CDS_1toLast52"]][[GENE]] <- integer(length = nchar(SEQ))
	
		if(length(grep(pattern = "U", SEQ)) > 0) SEQ <- RNAString(SEQ)

		for(j in 1:length(PATTERNS0)){
			PATTERN <- PATTERNS0[j]
			temp <- matchPattern(pattern = PATTERN, subject = DNAString(SEQ), 
					fixed = FALSE, max.mismatch = MAX.MISMATCH)
			if(length(temp@ranges@start) > 0){
				starts <- temp@ranges@start[temp@ranges@start > 0]
				matchList[["CDS_1toLast52"]][[GENE]][starts] <- 1
			}
		}
	}

	SEQ <- AsanoDF$CDS_52toLast[i]
	if(is.na(SEQ) == TRUE){
		matchList[["CDS_52toLast"]][[GENE]] <- integer(length = 0)
	}
	
	if(is.na(SEQ) == FALSE){
		matchList[["CDS_52toLast"]][[GENE]] <- integer(length = nchar(SEQ))
	
		if(length(grep(pattern = "U", SEQ)) > 0) SEQ <- RNAString(SEQ)

		for(j in 1:length(PATTERNS0)){
			PATTERN <- PATTERNS0[j]
			temp <- matchPattern(pattern = PATTERN, subject = DNAString(SEQ), 
					fixed = FALSE, max.mismatch = MAX.MISMATCH)
			if(length(temp@ranges@start) > 0){
				starts <- temp@ranges@start[temp@ranges@start > 0]
				matchList[["CDS_52toLast"]][[GENE]][starts] <- 1
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
dir.create(outDir)
setwd(outDir)
fileName <- paste(DAY, "_", MOTIF, "_m", MAX.MISMATCH, "_4940.xlsx", sep = "")
write.xlsx(compareDF, file = fileName)

setwd("../../RData")
dir.create(outDir)
setwd(outDir)
fileName <- paste(DAY, "_", MOTIF, "_m", MAX.MISMATCH, "_4940.RData", sep = "")
save(compareDF, file = fileName)


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
		if(is.na(SEQ) == FALSE){
			matchDF[[PATTERN]][i] <- 0
			temp <- matchPattern(pattern = PATTERN, subject = DNAString(SEQ), 
					fixed = FALSE, max.mismatch = MAX.MISMATCH)
			if(length(temp@ranges@start) > 0){
				starts <- temp@ranges@start[temp@ranges@start > 0]
				matchDF[[PATTERN]][i] <- length(starts) / nchar(SEQ) * 51
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
		if(is.na(SEQ) == FALSE){
			matchDF[[PATTERN]][i] <- 0
			temp <- matchPattern(pattern = PATTERN, subject = DNAString(SEQ), 
					fixed = FALSE, max.mismatch = MAX.MISMATCH)
			if(length(temp@ranges@start) > 0){
				starts <- temp@ranges@start[temp@ranges@start > 0]
				matchDF[[PATTERN]][i] <- length(starts) / nchar(SEQ) * 51
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
		if(is.na(SEQ) == FALSE){
			matchDF[[PATTERN]][i] <- 0
			temp <- matchPattern(pattern = PATTERN, subject = DNAString(SEQ), 
					fixed = FALSE, max.mismatch = MAX.MISMATCH)
			if(length(temp@ranges@start) > 0){
				starts <- temp@ranges@start[temp@ranges@start > 0]
				matchDF[[PATTERN]][i] <- length(starts) / nchar(SEQ) * 51
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
		if(is.na(SEQ) == FALSE){
			matchDF[[PATTERN]][i] <- 0
			temp <- matchPattern(pattern = PATTERN, subject = DNAString(SEQ), 
					fixed = FALSE, max.mismatch = MAX.MISMATCH)
			if(length(temp@ranges@start) > 0){
				starts <- temp@ranges@start[temp@ranges@start > 0]
				matchDF[[PATTERN]][i] <- length(starts) / nchar(SEQ) * 51
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
		if(is.na(SEQ) == FALSE){
			matchDF[[PATTERN]][i] <- 0
			temp <- matchPattern(pattern = PATTERN, subject = DNAString(SEQ), 
					fixed = FALSE, max.mismatch = MAX.MISMATCH)
			if(length(temp@ranges@start) > 0){
				starts <- temp@ranges@start[temp@ranges@start > 0]
				matchDF[[PATTERN]][i] <- length(starts) / nchar(SEQ) * 51
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
		if(is.na(SEQ) == FALSE){
			matchDF[[PATTERN]][i] <- 0
			temp <- matchPattern(pattern = PATTERN, subject = DNAString(SEQ), 
					fixed = FALSE, max.mismatch = MAX.MISMATCH)
			if(length(temp@ranges@start) > 0){
				starts <- temp@ranges@start[temp@ranges@start > 0]
				matchDF[[PATTERN]][i] <- length(starts) / nchar(SEQ) * 51
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
		if(is.na(SEQ) == FALSE){
			matchDF[[PATTERN]][i] <- 0
			temp <- matchPattern(pattern = PATTERN, subject = DNAString(SEQ), 
					fixed = FALSE, max.mismatch = MAX.MISMATCH)
			if(length(temp@ranges@start) > 0){
				starts <- temp@ranges@start[temp@ranges@start > 0]
				matchDF[[PATTERN]][i] <- length(starts) / nchar(SEQ) * 51
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
		if(is.na(SEQ) == FALSE){
			matchDF[[PATTERN]][i] <- 0
			temp <- matchPattern(pattern = PATTERN, subject = DNAString(SEQ), 
					fixed = FALSE, max.mismatch = MAX.MISMATCH)
			if(length(temp@ranges@start) > 0){
				starts <- temp@ranges@start[temp@ranges@start > 0]
				matchDF[[PATTERN]][i] <- length(starts) / nchar(SEQ) * 51
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
setwd(outDir)
fileName <- paste(DAY, "_matchDFlist_", MOTIF, "_m", MAX.MISMATCH, "_4940.xlsx", sep = "")
write.xlsx(matchDFlist, file = fileName)

setwd("../../RData")
setwd(outDir)
fileName <- paste(DAY, "_matchDFlist_", MOTIF, "_m", MAX.MISMATCH, "_4940.RData", sep = "")
save(matchDFlist, file = fileName)

setwd(HOMEDIR)
setwd("data")
setwd("motifHits")
setwd("RData")
load(file = "20250418_NewDF.RData")

setwd(outDir)
fileName <- paste(DAY, "_matchDFlist_", MOTIF, "_m", MAX.MISMATCH, "_4940.RData", sep = "")
load(file = fileName)


NewDF$tif34down_ALLtop50 <- NewDF$tif34_down_rank <= 50
NewDF$tif34down_ALLtop100 <- NewDF$tif34_down_rank <= 100
NewDF$tif34down_ALLtop500 <- NewDF$tif34_down_rank <= 500
NewDF$tif34down_ALLtop1000 <- NewDF$tif34_down_rank <= 1000

NewDF$tif34down_ALLbottom50 <- NewDF$tif34_down_rank >= nrow(NewDF) - 49
NewDF$tif34down_ALLbottom100 <- NewDF$tif34_down_rank >= nrow(NewDF) - 99
NewDF$tif34down_ALLbottom500 <- NewDF$tif34_down_rank >= nrow(NewDF) - 499
NewDF$tif34down_ALLbottom1000 <- NewDF$tif34_down_rank >= nrow(NewDF) - 999

NewDF$tif34down_ALLmiddle50 <- NewDF$tif34_down_rank >= round(nrow(NewDF)/2, digits = 0) - 25 & 
				NewDF$tif34_down_rank <= round(nrow(NewDF)/2, digits = 0) + 24
NewDF$tif34down_ALLmiddle100 <- NewDF$tif34_down_rank >= round(nrow(NewDF)/2, digits = 0) - 50 & 
				NewDF$tif34_down_rank <= round(nrow(NewDF)/2, digits = 0) + 49
NewDF$tif34down_ALLmiddle500 <- NewDF$tif34_down_rank >= round(nrow(NewDF)/2, digits = 0) - 250 & 
				NewDF$tif34_down_rank <= round(nrow(NewDF)/2, digits = 0) + 249
NewDF$tif34down_ALLmiddle1000 <- NewDF$tif34_down_rank >= round(nrow(NewDF)/2, digits = 0) - 500 & 
				NewDF$tif34_down_rank <= round(nrow(NewDF)/2, digits = 0) + 499

NewDF2 <- NewDF[, c(1, 6, 15, 20, 43:54, 33:42)]	# 20241216 新しい範囲を追加する

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

setwd("XLSX")
setwd(outDir)
fileName <- paste(DAY, "_NewDF2_", MOTIF, "_m", MAX.MISMATCH, ".xlsx", sep = "")
write.xlsx(NewDF2, file = fileName)

setwd("../../RData")
setwd(outDir)
fileName <- paste(DAY, "_NewDF2_", MOTIF, "_m", MAX.MISMATCH, ".RData", sep = "")
save(NewDF2, file = fileName)


setwd(HOMEDIR)
setwd("data")
setwd("motifHits")
setwd("RData")
setwd(outDir)
fileName <- paste(DAY, "_NewDF2_", MOTIF, "_m", MAX.MISMATCH, ".RData", sep = "")
load(file = fileName)	# NewDF2

setwd("../../TXT")
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

ALLtop50_hit_first51 <- subset(DF_hit_first51, tif34down_ALLtop50 == TRUE)$hit_first51
ALLtop100_hit_first51 <- subset(DF_hit_first51, tif34down_ALLtop100 == TRUE)$hit_first51
ALLtop500_hit_first51 <- subset(DF_hit_first51, tif34down_ALLtop500 == TRUE)$hit_first51
ALLtop1000_hit_first51 <- subset(DF_hit_first51, tif34down_ALLtop1000 == TRUE)$hit_first51

ALLbottom50_hit_first51 <- subset(DF_hit_first51, tif34down_ALLbottom50 == TRUE)$hit_first51
ALLbottom100_hit_first51 <- subset(DF_hit_first51, tif34down_ALLbottom100 == TRUE)$hit_first51
ALLbottom500_hit_first51 <- subset(DF_hit_first51, tif34down_ALLbottom500 == TRUE)$hit_first51
ALLbottom1000_hit_first51 <- subset(DF_hit_first51, tif34down_ALLbottom1000 == TRUE)$hit_first51

ALLmiddle50_hit_first51 <- subset(DF_hit_first51, tif34down_ALLmiddle50 == TRUE)$hit_first51
ALLmiddle100_hit_first51 <- subset(DF_hit_first51, tif34down_ALLmiddle100 == TRUE)$hit_first51
ALLmiddle500_hit_first51 <- subset(DF_hit_first51, tif34down_ALLmiddle500 == TRUE)$hit_first51
ALLmiddle1000_hit_first51 <- subset(DF_hit_first51, tif34down_ALLmiddle1000 == TRUE)$hit_first51

targets <- c("hit_first51", "hit_middle51", "hit_last51", "hit_UTR5", "hit_UTR3", "hit_CDS_total", "hit_CDS_1toLast52", "hit_CDS_52toLast")
pvaluesDF <- data.frame(target = targets, All_vs_RPS = as.numeric(NA), All_vs_RPL = as.numeric(NA), 
			All_vs_top50 = as.numeric(NA), All_vs_bottom50 = as.numeric(NA), All_vs_S4_246 = as.numeric(NA), 
			All_vs_HSR64 = as.numeric(NA), 
			All_vs_ALLtop50 = as.numeric(NA), All_vs_ALLtop100 = as.numeric(NA), All_vs_ALLtop500 = as.numeric(NA), All_vs_ALLtop1000 = as.numeric(NA),
			All_vs_ALLbottom50 = as.numeric(NA), All_vs_ALLbottom100 = as.numeric(NA), All_vs_ALLbottom500 = as.numeric(NA), All_vs_ALLbottom1000 = as.numeric(NA),
			All_vs_ALLmiddle50 = as.numeric(NA), All_vs_ALLmiddle100 = as.numeric(NA), All_vs_ALLmiddle500 = as.numeric(NA), All_vs_ALLmiddle1000 = as.numeric(NA))

meanvaluesDF <- data.frame(target = targets, All = as.numeric(NA), RPS = as.numeric(NA), RPL = as.numeric(NA), 
			top50 = as.numeric(NA), bottom50 = as.numeric(NA), S4_246 = as.numeric(NA), HSR64 = as.numeric(NA), 
			ALLtop50 = as.numeric(NA), ALLtop100 = as.numeric(NA), ALLtop500 = as.numeric(NA), ALLtop1000 = as.numeric(NA), 
			ALLbottom50 = as.numeric(NA), ALLbottom100 = as.numeric(NA), ALLbottom500 = as.numeric(NA), ALLbottom1000 = as.numeric(NA), 
			ALLmiddle50 = as.numeric(NA), ALLmiddle100 = as.numeric(NA), ALLmiddle500 = as.numeric(NA), ALLmiddle1000 = as.numeric(NA))

sdvaluesDF <- data.frame(target = targets, All = as.numeric(NA), RPS = as.numeric(NA), RPL = as.numeric(NA), 
			top50 = as.numeric(NA), bottom50 = as.numeric(NA), S4_246 = as.numeric(NA), HSR64 = as.numeric(NA), 
			ALLtop50 = as.numeric(NA), ALLtop100 = as.numeric(NA), ALLtop500 = as.numeric(NA), ALLtop1000 = as.numeric(NA), 
			ALLbottom50 = as.numeric(NA), ALLbottom100 = as.numeric(NA), ALLbottom500 = as.numeric(NA), ALLbottom1000 = as.numeric(NA), 
			ALLmiddle50 = as.numeric(NA), ALLmiddle100 = as.numeric(NA), ALLmiddle500 = as.numeric(NA), ALLmiddle1000 = as.numeric(NA))

for(t in 1:length(targets)){
	target <- targets[t]

	pvalues <- numeric(length = 18)

	# target <- "hit_first51"
	DF_hit <- subset(NewDF2, is.na(NewDF2[[target]]) == FALSE)
	RPS_hit <- subset(DF_hit, RPS == TRUE)[[target]]
	RPL_hit <- subset(DF_hit, RPL == TRUE)[[target]]
	top50_hit <- subset(DF_hit, tif34down_top50 == TRUE)[[target]]
	bottom50_hit <- subset(DF_hit, tif34down_bottom50 == TRUE)[[target]]
	S4_246_hit <- subset(DF_hit, S4_246 == TRUE)[[target]]
	All_hit <- DF_hit[[target]]
	HSR64_hit <- subset(DF_hit, Top64 == TRUE)[[target]]

	ALLtop50_hit <- subset(DF_hit, tif34down_ALLtop50 == TRUE)[[target]]
	ALLtop100_hit <- subset(DF_hit, tif34down_ALLtop100 == TRUE)[[target]]
	ALLtop500_hit <- subset(DF_hit, tif34down_ALLtop500 == TRUE)[[target]]
	ALLtop1000_hit <- subset(DF_hit, tif34down_ALLtop1000 == TRUE)[[target]]
	ALLbottom50_hit <- subset(DF_hit, tif34down_ALLbottom50 == TRUE)[[target]]
	ALLbottom100_hit <- subset(DF_hit, tif34down_ALLbottom100 == TRUE)[[target]]
	ALLbottom500_hit <- subset(DF_hit, tif34down_ALLbottom500 == TRUE)[[target]]
	ALLbottom1000_hit <- subset(DF_hit, tif34down_ALLbottom1000 == TRUE)[[target]]
	ALLmiddle50_hit <- subset(DF_hit, tif34down_ALLmiddle50 == TRUE)[[target]]
	ALLmiddle100_hit <- subset(DF_hit, tif34down_ALLmiddle100 == TRUE)[[target]]
	ALLmiddle500_hit <- subset(DF_hit, tif34down_ALLmiddle500 == TRUE)[[target]]
	ALLmiddle1000_hit <- subset(DF_hit, tif34down_ALLmiddle1000 == TRUE)[[target]]

	meanvaluesDF$All[t] <- mean(All_hit)
	meanvaluesDF$RPS[t] <- mean(RPS_hit)
	meanvaluesDF$RPL[t] <- mean(RPL_hit)
	meanvaluesDF$top50[t] <- mean(top50_hit)
	meanvaluesDF$bottom50[t] <- mean(bottom50_hit)
	meanvaluesDF$S4_246[t] <- mean(S4_246_hit)
	meanvaluesDF$HSR64[t] <- mean(HSR64_hit)

	meanvaluesDF$ALLtop50[t] <- mean(ALLtop50_hit)
	meanvaluesDF$ALLtop100[t] <- mean(ALLtop100_hit)
	meanvaluesDF$ALLtop500[t] <- mean(ALLtop500_hit)
	meanvaluesDF$ALLtop1000[t] <- mean(ALLtop1000_hit)
	meanvaluesDF$ALLbottom50[t] <- mean(ALLbottom50_hit)
	meanvaluesDF$ALLbottom100[t] <- mean(ALLbottom100_hit)
	meanvaluesDF$ALLbottom500[t] <- mean(ALLbottom500_hit)
	meanvaluesDF$ALLbottom1000[t] <- mean(ALLbottom1000_hit)
	meanvaluesDF$ALLmiddle50[t] <- mean(ALLmiddle50_hit)
	meanvaluesDF$ALLmiddle100[t] <- mean(ALLmiddle100_hit)
	meanvaluesDF$ALLmiddle500[t] <- mean(ALLmiddle500_hit)
	meanvaluesDF$ALLmiddle1000[t] <- mean(ALLmiddle1000_hit)

	sdvaluesDF$All[t] <- sd(All_hit)
	sdvaluesDF$RPS[t] <- sd(RPS_hit)
	sdvaluesDF$RPL[t] <- sd(RPL_hit)
	sdvaluesDF$top50[t] <- sd(top50_hit)
	sdvaluesDF$bottom50[t] <- sd(bottom50_hit)
	sdvaluesDF$S4_246[t] <- sd(S4_246_hit)
	sdvaluesDF$HSR64[t] <- sd(HSR64_hit)

	sdvaluesDF$ALLtop50[t] <- sd(ALLtop50_hit)
	sdvaluesDF$ALLtop100[t] <- sd(ALLtop100_hit)
	sdvaluesDF$ALLtop500[t] <- sd(ALLtop500_hit)
	sdvaluesDF$ALLtop1000[t] <- sd(ALLtop1000_hit)
	sdvaluesDF$ALLbottom50[t] <- sd(ALLbottom50_hit)
	sdvaluesDF$ALLbottom100[t] <- sd(ALLbottom100_hit)
	sdvaluesDF$ALLbottom500[t] <- sd(ALLbottom500_hit)
	sdvaluesDF$ALLbottom1000[t] <- sd(ALLbottom1000_hit)
	sdvaluesDF$ALLmiddle50[t] <- sd(ALLmiddle50_hit)
	sdvaluesDF$ALLmiddle100[t] <- sd(ALLmiddle100_hit)
	sdvaluesDF$ALLmiddle500[t] <- sd(ALLmiddle500_hit)
	sdvaluesDF$ALLmiddle1000[t] <- sd(ALLmiddle1000_hit)
	
	cat("## ", target, "\n", sep = "")
	cat("#  target==TRUE:         ", nrow(DF_hit), "\n", sep = "")
	cat("#  length(RPS_hit):      ", length(RPS_hit), "\n", sep = "")
	cat("#  length(RPL_hit):      ", length(RPL_hit), "\n", sep = "")
	cat("#  length(top50_hit):    ", length(top50_hit), "\n", sep = "")
	cat("#  length(bottom50_hit): ", length(bottom50_hit), "\n", sep = "")
	cat("#  length(S4_246_hit):   ", length(S4_246_hit), "\n", sep = "")
	cat("#  length(All_hit):      ", length(All_hit), "\n", sep = "")
	cat("#  length(HSR64_hit):      ", length(HSR64_hit), "\n", sep = "")

	cat("#  length(ALLtop50_hit):      ", length(ALLtop50_hit), "\n", sep = "")
	cat("#  length(ALLtop100_hit):      ", length(ALLtop100_hit), "\n", sep = "")
	cat("#  length(ALLtop500_hit):      ", length(ALLtop500_hit), "\n", sep = "")
	cat("#  length(ALLtop1000_hit):      ", length(ALLtop1000_hit), "\n", sep = "")
	cat("#  length(ALLbottom50_hit):      ", length(ALLbottom50_hit), "\n", sep = "")
	cat("#  length(ALLbottom100_hit):      ", length(ALLbottom100_hit), "\n", sep = "")
	cat("#  length(ALLbottom500_hit):      ", length(ALLbottom500_hit), "\n", sep = "")
	cat("#  length(ALLbottom1000_hit):      ", length(ALLbottom1000_hit), "\n", sep = "")
	cat("#  length(ALLmiddle50_hit):      ", length(ALLmiddle50_hit), "\n", sep = "")
	cat("#  length(ALLmiddle100_hit):      ", length(ALLmiddle100_hit), "\n", sep = "")
	cat("#  length(ALLmiddle500_hit):      ", length(ALLmiddle500_hit), "\n", sep = "")
	cat("#  length(ALLmiddle1000_hit):      ", length(ALLmiddle1000_hit), "\n", sep = "")

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
	
	cat("\n## t-test All_hit vs ALLtop50_hit\n")
	temp <- t.test(x = All_hit, y = ALLtop50_hit)
	pvalues[7] <- temp$p.value
	x <- paste("# ", temp)
	for(i in 1:length(x)){
		cat(x[i], "\n")
	}
	
	cat("\n## t-test All_hit vs ALLtop100_hit\n")
	temp <- t.test(x = All_hit, y = ALLtop100_hit)
	pvalues[8] <- temp$p.value
	x <- paste("# ", temp)
	for(i in 1:length(x)){
		cat(x[i], "\n")
	}
	
	cat("\n## t-test All_hit vs ALLtop500_hit\n")
	temp <- t.test(x = All_hit, y = ALLtop500_hit)
	pvalues[9] <- temp$p.value
	x <- paste("# ", temp)
	for(i in 1:length(x)){
		cat(x[i], "\n")
	}
	
	cat("\n## t-test All_hit vs ALLtop1000_hit\n")
	temp <- t.test(x = All_hit, y = ALLtop1000_hit)
	pvalues[10] <- temp$p.value
	x <- paste("# ", temp)
	for(i in 1:length(x)){
		cat(x[i], "\n")
	}
	
	cat("\n## t-test All_hit vs ALLbottom50_hit\n")
	temp <- t.test(x = All_hit, y = ALLbottom50_hit)
	pvalues[11] <- temp$p.value
	x <- paste("# ", temp)
	for(i in 1:length(x)){
		cat(x[i], "\n")
	}
	
	cat("\n## t-test All_hit vs ALLbottom100_hit\n")
	temp <- t.test(x = All_hit, y = ALLbottom100_hit)
	pvalues[12] <- temp$p.value
	x <- paste("# ", temp)
	for(i in 1:length(x)){
		cat(x[i], "\n")
	}
	
	cat("\n## t-test All_hit vs ALLbottom500_hit\n")
	temp <- t.test(x = All_hit, y = ALLbottom500_hit)
	pvalues[13] <- temp$p.value
	x <- paste("# ", temp)
	for(i in 1:length(x)){
		cat(x[i], "\n")
	}
	
	cat("\n## t-test All_hit vs ALLbottom1000_hit\n")
	temp <- t.test(x = All_hit, y = ALLbottom1000_hit)
	pvalues[14] <- temp$p.value
	x <- paste("# ", temp)
	for(i in 1:length(x)){
		cat(x[i], "\n")
	}
	
	cat("\n## t-test All_hit vs ALLmiddle50_hit\n")
	temp <- t.test(x = All_hit, y = ALLmiddle50_hit)
	pvalues[15] <- temp$p.value
	x <- paste("# ", temp)
	for(i in 1:length(x)){
		cat(x[i], "\n")
	}
	
	cat("\n## t-test All_hit vs ALLmiddle100_hit\n")
	temp <- t.test(x = All_hit, y = ALLmiddle100_hit)
	pvalues[16] <- temp$p.value
	x <- paste("# ", temp)
	for(i in 1:length(x)){
		cat(x[i], "\n")
	}
	
	cat("\n## t-test All_hit vs ALLmiddle500_hit\n")
	temp <- t.test(x = All_hit, y = ALLmiddle500_hit)
	pvalues[17] <- temp$p.value
	x <- paste("# ", temp)
	for(i in 1:length(x)){
		cat(x[i], "\n")
	}
	
	cat("\n## t-test All_hit vs ALLmiddle1000_hit\n")
	temp <- t.test(x = All_hit, y = ALLmiddle1000_hit)
	pvalues[18] <- temp$p.value
	x <- paste("# ", temp)
	for(i in 1:length(x)){
		cat(x[i], "\n")
	}

	pvaluesDF[t, 2:19] <- pvalues
	cat("\n\n")
}

pvaluesDFsig <- pvaluesDF
pvaluesDFsig[, 2:19][pvaluesDFsig[, 2:19] > 0.05] <- NA

outList <- list()
outList$pvalues <- pvaluesDF
outList$significant <- pvaluesDFsig
outList$mean <- meanvaluesDF
outList$sd <- sdvaluesDF

setwd(HOMEDIR)
setwd("data")
setwd("motifHits")
setwd("XLSX")
setwd(outDir)
fileName <- paste(DAY, "_motifHits_4940_", MOTIF, "_m", MAX.MISMATCH, ".xlsx", sep = "")
write.xlsx(outList, file = fileName)







