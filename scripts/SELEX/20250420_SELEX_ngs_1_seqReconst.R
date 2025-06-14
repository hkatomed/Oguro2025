## original note: 20250420_SELEX_ngs_1_seqReconst.txt


library(Biostrings)
HOMEDIR <- ""

WIDTH <- 4

setwd(HOMEDIR)
setwd("data")
setwd("SELEX")
setwd("onePercent")

seq01uniq <- read.table(file = "01_count_data.txt", header = TRUE, sep = "\t")
seq02uniq <- read.table(file = "02_count_data.txt", header = TRUE, sep = "\t")

seq01reconst <- character(length = sum(seq01uniq$count))
l <- 0
for(i in 1:nrow(seq01uniq)){
	SEQ <- seq01uniq$sequence[i]
	COUNT <- seq01uniq$count[i]
	START <- l+1
	END <- l+COUNT
	seq01reconst[START:END] <- SEQ
	l <- END
}
seq01reconst <- DNAStringSet(seq01reconst)
names(seq01reconst) <- paste("seq01_", 1:length(seq01reconst), sep = "")

seq02reconst <- character(length = sum(seq02uniq$count))
l <- 0
for(i in 1:nrow(seq02uniq)){
	SEQ <- seq02uniq$sequence[i]
	COUNT <- seq02uniq$count[i]
	START <- l+1
	END <- l+COUNT
	seq02reconst[START:END] <- SEQ
	l <- END
}
seq02reconst <- DNAStringSet(seq02reconst)
names(seq02reconst) <- paste("seq02_", 1:length(seq02reconst), sep = "")

for(i in 1:10){
	cat(sample(x = c("A", "C", "G", "T"), size = 1, prob = c(0.25, 0.25, 0.25, 0.25)), ",")
}

set.seed(1)
for(i in 1:10){
	cat(sample(x = c("A", "C", "G", "T"), size = 1, prob = c(0.25, 0.25, 0.25, 0.25)), ",")
}

seq01control <- seq01reconst
names(seq01control) <- paste("seq01_c", 1:length(seq01reconst), sep = "")
set.seed(1)
for(i in 1:length(seq01control)){
	if(i %% 100 == 0) cat(i, ",")
	if(i %% 1000 == 0) cat(i, "\n")
	SEQ <- seq01control[[i]]
	LEN <- length(SEQ)
	newSEQ <- SEQ
	for(j in 1:LEN){
		newSEQ[j] <- sample(x = c("A", "C", "G", "T"), 
				size = 1, prob = c(0.25, 0.25, 0.25, 0.25))
	}
	seq01control[[i]] <- newSEQ
}

seq02control <- seq02reconst
names(seq02control) <- paste("seq02_c", 1:length(seq02reconst), sep = "")
set.seed(1)
for(i in 1:length(seq02control)){
	if(i %% 100 == 0) cat(i, ",")
	if(i %% 1000 == 0) cat(i, "\n")
	SEQ <- seq02control[[i]]
	LEN <- length(SEQ)
	newSEQ <- SEQ
	for(j in 1:LEN){
		newSEQ[j] <- sample(x = c("A", "C", "G", "T"), 
				size = 1, prob = c(0.25, 0.25, 0.25, 0.25))
	}
	seq02control[[i]] <- newSEQ
}


setwd(HOMEDIR)
setwd("data")
setwd("SELEX")
setwd("onePercent")
writeXStringSet(seq01reconst, filepath = "seq01reconst.fasta.gz", compress = TRUE)
writeXStringSet(seq02reconst, filepath = "seq02reconst.fasta.gz", compress = TRUE)
writeXStringSet(seq01control, filepath = "seq01control.fasta.gz", compress = TRUE)
writeXStringSet(seq02control, filepath = "seq02control.fasta.gz", compress = TRUE)




