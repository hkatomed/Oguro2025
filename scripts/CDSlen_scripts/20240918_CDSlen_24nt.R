## This file is based on 20240614_Asano_motif_1_24nt.txt


HOMEDIR <- "xxx/yyy/CDS_first12to51"	# Define home directory
library(openxlsx)
setwd(HOMEDIR)
setwd("20240522_AsanoTE")

####################
RPKMtable <- read.xlsx(xlsxFile = "Yeast_CHX_KSU_RPKM_table_2_hk2.xlsx")

# > str(RPKMtable[,1:12])
# 'data.frame':	4940 obs. of  12 variables:
#  $ Gene              : chr  "YJL153C" "YOL154W" "YDR206W" "YER091C" ...
#  $ Original.Count.C  : num  132 227 125 5951 4177 ...
#  $ Original.Count.TV1: num  10 31 18 654 421 ...
#  $ Original.Count.TV2: num  12 27 14 968 766 ...
#  $ Original.Count.TT : num  41 84 45 2462 2028 ...
#  $ RPKM.C            : num  22.9 84.2 13.1 718.7 1006.3 ...
#  $ RPKM.TV1          : num  4.29 28.43 4.66 195.22 250.69 ...
#  $ RPKM.TV2          : num  4.88 23.43 3.43 273.48 431.69 ...
#  $ RPKM.TT           : num  9.1 39.81 6.02 379.82 624.11 ...
#  $ tif34.down.fold   : num  5 3.25 3.24 3.07 2.95 ...
#  $ tif35.up.fold     : num  1.98 1.54 1.49 1.62 1.83 ...
#  $ down/up           : num  0.246 0.238 0.218 0.3 0.425 ...


# > length(unique(RPKMtable$Gene))
# [1] 4940

######################
library(Biostrings)
setwd(HOMEDIR)
setwd("SGD")
ORFs <- readDNAStringSet(filepath = "orf_coding.fasta")

# > ORFs
# DNAStringSet object of length 6034:
#        width seq                                           names               
#    [1]   363 ATGGTCAAATTAACTTCAATC...ATCTACACTATCGCAAACTAG YAL068C PAU8 SGDI...
#    [2]   228 ATGCCAATTATAGGGGTGCCG...TTGGGAGTCGTATACTGTTAG YAL067W-A YAL067W...
#    [3]  1782 ATGTATTCAATTGTTAAAGAG...TCAGTATCTGATGAAAAATAA YAL067C SEO1 SGDI...
#    [4]   387 ATGAACAGTGCTACCAGTGAG...TTGCTGGCAATCGTATGGTAA YAL065C YAL065C S...
#    [5]   381 ATGGCAGGTGAAGCAGTTTCG...ATGGATTCAGTGCACACATGA YAL064W-B YAL064W...
#    ...   ... ...
# [6030]  1197 ATGAAATTAAAATTATTAAAT...GTTAAATTAAACTTTATTTAA Q0140 VAR1 SGDID:...
# [6031]   708 ATGAAAAATATTAAAAAAAAT...TCCGAAACTTTTTTAAAATAA Q0160 SCEI SGDID:...
# [6032]   756 ATGTTAGATTTATTAAGATTA...GAATGATTAAATGAACAATAA Q0250 COX2 SGDID:...
# [6033]  1419 ATGATTAAATGAACAATAATT...AAAAATATTTTATTAGATTAA Q0255 Q0255 SGDID...
# [6034]   810 ATGACACATTTAGAAAGAAGT...TTCTACTGATGAGGAGTCTAA Q0275 COX3 SGDID:...

orf.names <- names(ORFs)
for(i in 1:length(orf.names)){
	temp <- orf.names[i]
	orf.names[i] <- strsplit(temp, split = " ")[[1]][1]
}

# > head(orf.names)
# [1] "YAL068C"   "YAL067W-A" "YAL067C"   "YAL065C"   "YAL064W-B" "YAL064C-A"

names(ORFs) <- orf.names

# > ORFs
# DNAStringSet object of length 6034:
#        width seq                                           names               
#    [1]   363 ATGGTCAAATTAACTTCAATC...ATCTACACTATCGCAAACTAG YAL068C
#    [2]   228 ATGCCAATTATAGGGGTGCCG...TTGGGAGTCGTATACTGTTAG YAL067W-A
#    [3]  1782 ATGTATTCAATTGTTAAAGAG...TCAGTATCTGATGAAAAATAA YAL067C
#    [4]   387 ATGAACAGTGCTACCAGTGAG...TTGCTGGCAATCGTATGGTAA YAL065C
#    [5]   381 ATGGCAGGTGAAGCAGTTTCG...ATGGATTCAGTGCACACATGA YAL064W-B
#    ...   ... ...
# [6030]  1197 ATGAAATTAAAATTATTAAAT...GTTAAATTAAACTTTATTTAA Q0140
# [6031]   708 ATGAAAAATATTAAAAAAAAT...TCCGAAACTTTTTTAAAATAA Q0160
# [6032]   756 ATGTTAGATTTATTAAGATTA...GAATGATTAAATGAACAATAA Q0250
# [6033]  1419 ATGATTAAATGAACAATAATT...AAAAATATTTTATTAGATTAA Q0255
# [6034]   810 ATGACACATTTAGAAAGAAGT...TTCTACTGATGAGGAGTCTAA Q0275

# > summary(width(ORFs))
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#      51     681    1194    1464    1871   14733 

# RPKMtable$CDS_51bp <- as.character(NA)
# RPKMtable$CDS_36bp <- as.character(NA)
# RPKMtable$CDS_33bp <- as.character(NA)
# RPKMtable$CDS_30bp <- as.character(NA)
RPKMtable$CDS_24bp <- as.character(NA)	
for(i in 1:nrow(RPKMtable)){
	geneName <- RPKMtable$Gene[i]
	index <- which(names(ORFs) == geneName)
	if(length(index) > 1) cat("i:", i, "length(index) > 1", "\n", sep = "")
	if(length(index) == 0) cat("i:", i, "length(index) == 0", "\n", sep = "")
	if(length(index) == 1){
		# RPKMtable$CDS_51bp[i] <- as.character(ORFs[[index]][1:51])
		# RPKMtable$CDS_36bp[i] <- as.character(ORFs[[index]][1:36])
		# RPKMtable$CDS_33bp[i] <- as.character(ORFs[[index]][1:33])
		# RPKMtable$CDS_30bp[i] <- as.character(ORFs[[index]][1:30])
		RPKMtable$CDS_24bp[i] <- as.character(ORFs[[index]][1:24])
	}
}

length(which(is.na(RPKMtable$CDS_24bp) == TRUE))

# > length(which(is.na(RPKMtable$CDS_24bp) == TRUE))
# [1] 49


##################
AsanoDF <- RPKMtable[,c(1, 94)]
head(AsanoDF, n = 3)

# > head(AsanoDF, n = 3)
#      Gene                 CDS_24bp
# 1 YJL153C ATGACAGAAGATAATATTGCTCCA
# 2 YOL154W ATGAAGTTCTCTTCCGGCAAATCT
# 3 YDR206W ATGGAACCATCGAATACCCAAAAA


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

######################################

matchList <- list()
for(i in 1:nrow(AsanoDF)){
	if(i%%100 == 0) cat(i, ", ", sep = "")
	if(i == nrow(AsanoDF)) cat("\n")
	GENE <- AsanoDF$Gene[i]
	matchList[[GENE]] <- integer(length = 51)
	SEQ <- AsanoDF$CDS_24bp[i]
	
	if(is.na(SEQ) == FALSE){
	
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
}

length(matchList)

# > length(matchList)
# [1] 4940


matchLines <- character(length = length(matchList) *2)
geneNames <- character(length = length(matchList) *2)
hit <- character(length = length(matchList) *2)
for(i in 1:length(matchList)){
	if(i%%100 == 0) cat(i, ", ", sep = "")
	if(i == length(matchList)) cat("\n")
	temp <- paste(matchList[[i]], collapse = "")
	matchLines[(i-1)*2+1] <- AsanoDF$CDS_24bp[i]
	matchLines[(i-1)*2+2] <- gsub(pattern = "0", replacement = " ", temp)
	geneNames[(i-1)*2+1] <- AsanoDF$Gene[i]
	geneNames[(i-1)*2+2] <- AsanoDF$Gene[i]
	if(sum(matchList[[i]]) > 0){
		hit[(i-1)*2+2] <- as.character(sum(matchList[[i]]))
	}
}

compareDF <- data.frame(Gene = geneNames, Seq_Match = matchLines, hit = hit)

setwd(HOMEDIR)
setwd("XLSX")
write.xlsx(compareDF, file = "20240614_AGYYGRY_6thTrial_4940_24bp.xlsx")




######################################

matchDF <- AsanoDF
for(j in 1:length(PATTERNS0)){
	PATTERN <- PATTERNS0[j]
	cat(j, ", ", PATTERN, "\n", sep = "")
	matchDF[[PATTERN]] <- as.integer(NA)
	for(i in 1:nrow(AsanoDF)){
		if(i%%100 == 0) cat(i, ", ", sep = "")
		if(i == nrow(AsanoDF)) cat("\n")
		SEQ <- AsanoDF$CDS_24bp[i]
		if(is.na(SEQ) == FALSE){
			matchDF[[PATTERN]][i] <- 0
			temp <- matchPattern(pattern = PATTERN, subject = DNAString(SEQ), 
					fixed = FALSE, max.mismatch = 0)
			if(length(temp@ranges@start) > 0){
				matchDF[[PATTERN]][i] <- length(temp@ranges@start)
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
		SEQ <- AsanoDF$CDS_24bp[i]
		if(is.na(SEQ) == FALSE){
			matchDF[[colName]][i] <- 0
			temp <- matchPattern(pattern = PATTERN, subject = DNAString(SEQ), 
					fixed = FALSE, max.mismatch = 1)
			if(length(temp@ranges@start) > 0){
				matchDF[[colName]][i] <- length(temp@ranges@start)
			}
		}
	}
}

# > colnames(matchDF)
#  [1] "Gene"      "CDS_24bp"  "AGYYGNN"   "NGYYGRN"   "NNYYGRY"   "NGYYGRY"  
#  [7] "NAYYGRY"   "NGYYARY"   "NAYYARY"   "AGYYGRN"   "GAYYGRN"   "AGYYGRY"  
# [13] "AGYYARY"   "GAYYGRY"   "GGYYARY"   "AAYYGGN"   "AGYYAGN"   "GGYYAGN"  
# [19] "AAYYAGN"   "GGYYGGY"   "AAYYGGY"   "AAYYAGY"   "GNYYGAY"   "GGNYGAY"  
# [25] "GGYNGAY"   "GNYYGGY"   "GGYYGNY"   "AANYGGY"   "AAYNGGY"   "AAYYNGY"  
# [31] "AAYYGNY"   "AGNYAAY"   "AGYNAAY"   "ANYYAGY"   "AGNYAGY"   "AGYNAGY"  
# [37] "GANYGAY"   "GANYGGY"   "GAYNGGY"   "AANYAGY"   "AAYNAGY"   "AAYYANY"  
# [43] "AGYYGRYm1"

matchDF$hit <- as.integer(0)
# matchDF$hit <- rowSums(matchDF[,3:43]) 

# > nrow(compareDF)
# [1] 9880
# > nrow(matchDF)
# [1] 4940
# > nrow(matchDF) *2
# [1] 9880

index <- seq(2, nrow(compareDF), 2)

# > head(index)
# [1]  2  4  6  8 10 12
# > tail(index)
# [1] 9870 9872 9874 9876 9878 9880

matchDF$hit <- as.integer(compareDF$hit[index])
matchDF$hit[is.na(matchDF$hit) == TRUE] <- 0
matchDF$hit[is.na(matchDF$CDS_24bp) == TRUE] <- NA


setwd(HOMEDIR)
setwd("XLSX")
write.xlsx(matchDF, file = "20240614_matchDF_6thTrial_4940_24bp.xlsx")


