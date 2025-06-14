## original note: 20250427_SELEX_ngs_motif_mismatch_1.R.txt

library(Biostrings)
HOMEDIR <- ""

setwd(HOMEDIR)
setwd("data")
setwd("SELEX")
setwd("onePercent")

seq01reconst <- readDNAStringSet(filepath = "seq01reconst.fasta.gz")
seq02reconst <- readDNAStringSet(filepath = "seq02reconst.fasta.gz")
seq01control <- readDNAStringSet(filepath = "seq01control.fasta.gz")
seq02control <- readDNAStringSet(filepath = "seq02control.fasta.gz")

getMotifHitsInt <- function(SEQ){
	SEQ <- DNAString(SEQ)
	motifHitsInt <- integer(length = length(SEQ))
	MAX.MISMATCH <- 0
	if(PATTERNS0[1] != ""){
	for(m in 1:length(PATTERNS0)){
		PATTERN <- PATTERNS0[m]
		temp <- matchPattern(pattern = PATTERN, SEQ, fixed = FALSE, max.mismatch = MAX.MISMATCH)
		index <- temp@ranges@start[temp@ranges@start > MAX.MISMATCH]
		motifHitsInt[index] <- 1
	}}

	if(PATTERNS1[1] != ""){
	MAX.MISMATCH <- 1
	for(m in 1:length(PATTERNS1)){
		PATTERN <- PATTERNS1[m]
		temp <- matchPattern(pattern = PATTERN, SEQ, fixed = FALSE, max.mismatch = MAX.MISMATCH)
		index <- temp@ranges@start[temp@ranges@start > MAX.MISMATCH]
		motifHitsInt[index] <- 1
	}}

	if(PATTERNS2[1] != ""){
	MAX.MISMATCH <- 2
	for(m in 1:length(PATTERNS2)){
		PATTERN <- PATTERNS2[m]
		temp <- matchPattern(pattern = PATTERN, SEQ, fixed = FALSE, max.mismatch = MAX.MISMATCH)
		index <- temp@ranges@start[temp@ranges@start > MAX.MISMATCH]
		motifHitsInt[index] <- 1
	}}
	return(sum(motifHitsInt))
}


getMotifHitsIntEach <- function(SEQ, PATTERN, MAX.MISMATCH){
	SEQ <- DNAString(SEQ)
	motifHitsInt <- integer(length = length(SEQ))
	temp <- matchPattern(pattern = PATTERN, SEQ, fixed = FALSE, max.mismatch = MAX.MISMATCH)
	index <- temp@ranges@start[temp@ranges@start > MAX.MISMATCH]
	motifHitsInt[index] <- 1
	return(sum(motifHitsInt))
}


setwd(HOMEDIR)
setwd("data")
setwd("motifHits")
setwd("TXT")

# PATTERNS0 <- readLines(con = "20240427_Final8base_m0_withNGTTARNN.txt")
PATTERNS0 <- readLines(con = "20240427_Final8base_m0.txt")
PATTERNS1 <- ""
PATTERNS2 <- readLines(con = "20240427_Final8base_m2.txt")


library(parallel)

SEQS <- as.character(seq01reconst)
motifHitsInts <- integer(length = length(SEQS))
motifHitsInts <- unlist(mcMap(f = getMotifHitsInt, SEQS, mc.cores = 6), use.names = FALSE)
hits.seq01reconst <- sum(motifHitsInts)

SEQS <- as.character(seq02reconst)
motifHitsInts <- integer(length(SEQS))
motifHitsInts <- unlist(mcMap(f = getMotifHitsInt, SEQS, mc.cores = 6), use.names = FALSE)
hits.seq02reconst <- sum(motifHitsInts)

SEQS <- as.character(seq01control)
motifHitsInts <- integer(length(SEQS))
motifHitsInts <- unlist(mcMap(f = getMotifHitsInt, SEQS, mc.cores = 6), use.names = FALSE)
hits.seq01control <- sum(motifHitsInts)

SEQS <- as.character(seq02control)
motifHitsInts <- integer(length(SEQS))
motifHitsInts <- unlist(mcMap(f = getMotifHitsInt, SEQS, mc.cores = 6), use.names = FALSE)
hits.seq02control <- sum(motifHitsInts)


outDF <- data.frame(motif = "merged", 
		seq01r.count = hits.seq01reconst, 
		seq02r.count = hits.seq02reconst, 
		seq01c.count = hits.seq01control, 
		seq02c.count = hits.seq02control)

MAX.MISMATCH <- 0
for(p in 1:length(PATTERNS0)){
	PATTERN <- PATTERNS0[p]
	
	SEQS <- as.character(seq01reconst)
	motifHitsInts <- integer(length(SEQS))
	motifHitsInts <- unlist(mcMap(f = getMotifHitsIntEach, 
		SEQS, PATTERN, MAX.MISMATCH = MAX.MISMATCH, mc.cores = 6), use.names = FALSE)
	hits.seq01reconst <- sum(motifHitsInts)

	SEQS <- as.character(seq02reconst)
	motifHitsInts <- integer(length(SEQS))
	motifHitsInts <- unlist(mcMap(f = getMotifHitsIntEach, 
		SEQS, PATTERN, MAX.MISMATCH = MAX.MISMATCH, mc.cores = 6), use.names = FALSE)
	hits.seq02reconst <- sum(motifHitsInts)
	
	SEQS <- as.character(seq01control)
	motifHitsInts <- integer(length(SEQS))
	motifHitsInts <- unlist(mcMap(f = getMotifHitsIntEach, 
		SEQS, PATTERN, MAX.MISMATCH = MAX.MISMATCH, mc.cores = 6), use.names = FALSE)
	hits.seq01control <- sum(motifHitsInts)

	SEQS <- as.character(seq02control)
	motifHitsInts <- integer(length(SEQS))
	motifHitsInts <- unlist(mcMap(f = getMotifHitsIntEach, 
		SEQS, PATTERN, MAX.MISMATCH = MAX.MISMATCH, mc.cores = 6), use.names = FALSE)
	hits.seq02control <- sum(motifHitsInts)

	tempDF <- data.frame(motif = PATTERN, 
			seq01r.count = hits.seq01reconst, 
			seq02r.count = hits.seq02reconst, 
			seq01c.count = hits.seq01control, 
			seq02c.count = hits.seq02control)

	outDF <- rbind(outDF, tempDF)
}


MAX.MISMATCH <- 1
if(PATTERNS1[1] != ""){
for(p in 1:length(PATTERNS1)){
	PATTERN <- PATTERNS1[p]
	
	SEQS <- as.character(seq01reconst)
	motifHitsInts <- integer(length(SEQS))
	motifHitsInts <- unlist(mcMap(f = getMotifHitsIntEach, 
		SEQS, PATTERN, MAX.MISMATCH = MAX.MISMATCH, mc.cores = 6), use.names = FALSE)
	hits.seq01reconst <- sum(motifHitsInts)

	SEQS <- as.character(seq02reconst)
	motifHitsInts <- integer(length(SEQS))
	motifHitsInts <- unlist(mcMap(f = getMotifHitsIntEach, 
		SEQS, PATTERN, MAX.MISMATCH = MAX.MISMATCH, mc.cores = 6), use.names = FALSE)
	hits.seq02reconst <- sum(motifHitsInts)
	
	SEQS <- as.character(seq01control)
	motifHitsInts <- integer(length(SEQS))
	motifHitsInts <- unlist(mcMap(f = getMotifHitsIntEach, 
		SEQS, PATTERN, MAX.MISMATCH = MAX.MISMATCH, mc.cores = 6), use.names = FALSE)
	hits.seq01control <- sum(motifHitsInts)

	SEQS <- as.character(seq02control)
	motifHitsInts <- integer(length(SEQS))
	motifHitsInts <- unlist(mcMap(f = getMotifHitsIntEach, 
		SEQS, PATTERN, MAX.MISMATCH = MAX.MISMATCH, mc.cores = 6), use.names = FALSE)
	hits.seq02control <- sum(motifHitsInts)

	tempDF <- data.frame(motif = paste(PATTERN, "_m1", sep = ""), 
			seq01r.count = hits.seq01reconst, 
			seq02r.count = hits.seq02reconst, 
			seq01c.count = hits.seq01control, 
			seq02c.count = hits.seq02control)

	outDF <- rbind(outDF, tempDF)
}
}


MAX.MISMATCH <- 2
if(PATTERNS2[1] != ""){
for(p in 1:length(PATTERNS2)){
	PATTERN <- PATTERNS2[p]
	
	SEQS <- as.character(seq01reconst)
	motifHitsInts <- integer(length(SEQS))
	motifHitsInts <- unlist(mcMap(f = getMotifHitsIntEach, 
		SEQS, PATTERN, MAX.MISMATCH = MAX.MISMATCH, mc.cores = 6), use.names = FALSE)
	hits.seq01reconst <- sum(motifHitsInts)

	SEQS <- as.character(seq02reconst)
	motifHitsInts <- integer(length(SEQS))
	motifHitsInts <- unlist(mcMap(f = getMotifHitsIntEach, 
		SEQS, PATTERN, MAX.MISMATCH = MAX.MISMATCH, mc.cores = 6), use.names = FALSE)
	hits.seq02reconst <- sum(motifHitsInts)
	
	SEQS <- as.character(seq01control)
	motifHitsInts <- integer(length(SEQS))
	motifHitsInts <- unlist(mcMap(f = getMotifHitsIntEach, 
		SEQS, PATTERN, MAX.MISMATCH = MAX.MISMATCH, mc.cores = 6), use.names = FALSE)
	hits.seq01control <- sum(motifHitsInts)

	SEQS <- as.character(seq02control)
	motifHitsInts <- integer(length(SEQS))
	motifHitsInts <- unlist(mcMap(f = getMotifHitsIntEach, 
		SEQS, PATTERN, MAX.MISMATCH = MAX.MISMATCH, mc.cores = 6), use.names = FALSE)
	hits.seq02control <- sum(motifHitsInts)

	tempDF <- data.frame(motif = paste(PATTERN, "_m2", sep = ""), 
			seq01r.count = hits.seq01reconst, 
			seq02r.count = hits.seq02reconst, 
			seq01c.count = hits.seq01control, 
			seq02c.count = hits.seq02control)

	outDF <- rbind(outDF, tempDF)
}
}

outDF$seq01.p.value <- unlist(Map(f = ppois, q = outDF$seq01r.count, lambda = outDF$seq01c.count, lower.tail = FALSE))
outDF$seq02.p.value <- unlist(Map(f = ppois, q = outDF$seq02r.count, lambda = outDF$seq02c.count, lower.tail = FALSE))


outDF.ordered01 <- outDF[order(outDF$seq01r.count, decreasing = TRUE), ]
outDF.ordered02 <- outDF[order(outDF$seq02r.count, decreasing = TRUE), ]

temp <- c(outDF$seq01.p.value, outDF$seq02.p.value)
temp <- temp[temp > 0]
summary(temp)
min(temp)


setwd(HOMEDIR)
setwd("data")
setwd("SELEX")
setwd("onePercent")
save(outDF, file = "20250427_SELEX_motifs_final8base_fixed.RData")
library(openxlsx)
write.xlsx(outDF, file = "20250427_SELEX_motifs_final8base_fixed.xlsx")



