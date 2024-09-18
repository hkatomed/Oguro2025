## This file is based on 20240615_Asano_motif_4_CDS_v3.txt


## Files used in the analysis stored in SGD_20240615
## (generated with the codes in 20240918_genomic_distribution_of_GUCG_motifs_1.txt)

# sc_genome.fasta.gz
# Intergenic.fasta.gz

# ncRNA.fasta.gz
# rRNA.fasta.gz
# snoRNA.fasta.gz
# snRNA.fasta.gz
# teloRNA.fasta.gz
# tRNA.fasta.gz

# CDS_total.fasta.gz
# CDS_over101.fasta.gz
# CDS_first51.fasta.gz
# CDS_last51.fasta.gz
# CDS_52toLast.fasta.gz
# CDS_1toLast52.fasta.gz

# UTR5_selected.fasta.gz
# UTR3_selected.fasta.gz

# transcripts.fasta.gz


HOMEDIR <- "xxx/yyy/motifDistribution"	# define home directory

##################################################
setwd(HOMEDIR)
setwd("TXT")

SELECTED <- readLines(con = "20230805_SELECTED.txt")		# 9 motifs removed
SELECTEDm1 <- readLines(con = "20230805_SELECTEDm1.txt")	# no change

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

PATTERNS0 <- SELECTED

PATTERNS1 <- SELECTEDm1



################################################

# targetNames <- c("ncRNA", "tRNA", "snoRNA", "snRNA", "teloRNA", "tRNA")
# targetNames <- c("CDS_total", "transcripts", "sc_genome", "Intergenic")
# targetNames <- c("CDS_over101", "CDS_first51", "CDS_last51", "CDS_52toLast", "CDS_1toLast52", "UTR5_selected", "UTR3_selected")
# targetNames <- c("UTR5", "UTR3")
targetNames <- "CDS_52toLast"
fileNames <- paste(targetNames, ".fasta.gz", sep = "")


library(Biostrings)
setwd("/Users/hrk_kato/Library/Mobile Documents/com~apple~CloudDocs/20230728_Asano_sc_motif")
setwd("SGD_20240615")


matchDF <- data.frame(target = targetNames)
for(j in 1:length(PATTERNS0)){
	PATTERN <- PATTERNS0[j]
	matchDF[[PATTERN]] <- as.numeric(NA)
}
for(j in 1:length(PATTERNS1)){
	PATTERN <- PATTERNS1[j]
	matchDF[[paste(PATTERN, "m1", sep = "")]] <- as.numeric(NA)
}
matchDF[["ALL"]] <- as.numeric(NA)

for(t in 1:length(targetNames)){
	cat("target: ", targetNames[t], "\n", sep = "")
	SEQs <- readDNAStringSet(filepath = fileNames[t])

	coordListTemp <- list()
	for(i in 1:length(SEQs)){
		coordListTemp[[names(SEQs)[i]]] <- integer(length = width(SEQs)[i])
	}

	# > str(coordList[1:3])
	# List of 3
	#  $ YNCA0001W: int [1:564] 0 0 0 0 0 0 0 0 0 0 ...
	#  $ YNCB0008W: int [1:3841] 0 0 0 0 0 0 0 0 0 0 ...
	#  $ YNCB0014W: int [1:2024] 0 0 0 0 0 0 0 0 0 0 ...
	
	coordListAll <- coordListTemp

	for(j in 1:length(PATTERNS0)){
		PATTERN <- PATTERNS0[j]
		cat(PATTERN, ", ", sep = "")
		coordList <- coordListTemp
		for(i in 1:length(SEQs)){
			if(i%%100 == 0) cat(i, ", ", sep = "")		
			SEQ <- SEQs[[i]]
			temp <- matchPattern(pattern = PATTERN, subject = DNAString(SEQ), 
					fixed = FALSE, max.mismatch = 0)
			if(length(temp@ranges@start) > 0){
				coordList[[i]][temp@ranges@start] <- 1
				coordListAll[[i]][temp@ranges@start] <- 1
			}
		}
		patternSum <- 0
		for(i in 1:length(SEQs)){
			patternSum <- patternSum + sum(coordList[[i]])
		}
		patternPer51 <- patternSum / sum(width(SEQs)) * 51
		matchDF[[PATTERN]][t] <- patternPer51
	}
	cat("\n")

	for(j in 1:length(PATTERNS1)){
		PATTERN <- PATTERNS1[j]
		cat(PATTERN, "m1, ", sep = "")
		coordList <- coordListTemp
		for(i in 1:length(SEQs)){
			if(i%%100 == 0) cat(i, ", ", sep = "")		
			SEQ <- SEQs[[i]]
			temp <- matchPattern(pattern = PATTERN, subject = DNAString(SEQ), 
					fixed = FALSE, max.mismatch = 1)
			if(length(temp@ranges@start) > 0){
				coordList[[i]][temp@ranges@start] <- 1
				coordListAll[[i]][temp@ranges@start] <- 1
			}
		}
		patternSum <- 0
		for(i in 1:length(SEQs)){
			patternSum <- patternSum + sum(coordList[[i]])
		}
		patternPer51 <- patternSum / sum(width(SEQs)) * 51
		matchDF[[paste(PATTERN, "m1", sep = "")]][t] <- patternPer51
	}
	cat("\n")

	patternSum <- 0
	for(i in 1:length(SEQs)){
		patternSum <- patternSum + sum(coordListAll[[i]])
	}
	patternPer51 <- patternSum / sum(width(SEQs)) * 51
	matchDF[["ALL"]][t] <- patternPer51
	cat("\n")
}

library(openxlsx)
setwd(HOMEDIR)
dir.create("XLSX")
setwd("XLSX")
write.xlsx(format(matchDF, scientific = FALSE), file = "20240615_CDS_6thTrial_v3.xlsx")



######################################



