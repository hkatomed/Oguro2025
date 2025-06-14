## original note: 20250420_make_AsanoDF_before_motifHits.R.txt


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
ORFs <- readDNAStringSet(filepath = "orf_coding.fasta.gz")

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
setwd("motifHits")
setwd("RData")
save(AsanoDF, file = "20250420_AsanoDF.RData")




