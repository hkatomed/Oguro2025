## original note: 20250419_SELEX_ngs_1_freq4.txt


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

seq01reconst.freq4c <- oligonucleotideFrequency(x = seq01reconst, width = WIDTH, as.prob = FALSE)
seq02reconst.freq4c <- oligonucleotideFrequency(x = seq02reconst, width = WIDTH, as.prob = FALSE)
seq01control.freq4c <- oligonucleotideFrequency(x = seq01control, width = WIDTH, as.prob = FALSE)
seq02control.freq4c <- oligonucleotideFrequency(x = seq02control, width = WIDTH, as.prob = FALSE)

seq01reconst.freq4cT <- t(seq01reconst.freq4c)
seq02reconst.freq4cT <- t(seq02reconst.freq4c)
seq01control.freq4cT <- t(seq01control.freq4c)
seq02control.freq4cT <- t(seq02control.freq4c)

seq01reconst.freq4c.combined <- rowSums(seq01reconst.freq4cT)
seq02reconst.freq4c.combined <- rowSums(seq02reconst.freq4cT)
seq01control.freq4c.combined <- rowSums(seq01control.freq4cT)
seq02control.freq4c.combined <- rowSums(seq02control.freq4cT)

seq01reconst.freq4p <- oligonucleotideFrequency(x = seq01reconst, width = WIDTH, as.prob = TRUE)
seq02reconst.freq4p <- oligonucleotideFrequency(x = seq02reconst, width = WIDTH, as.prob = TRUE)
seq01control.freq4p <- oligonucleotideFrequency(x = seq01control, width = WIDTH, as.prob = TRUE)
seq02control.freq4p <- oligonucleotideFrequency(x = seq02control, width = WIDTH, as.prob = TRUE)

seq01reconst.freq4pT <- t(seq01reconst.freq4p)
seq02reconst.freq4pT <- t(seq02reconst.freq4p)
seq01control.freq4pT <- t(seq01control.freq4p)
seq02control.freq4pT <- t(seq02control.freq4p)

seq01reconst.freq4p.combined <- rowSums(seq01reconst.freq4pT) / sum(seq01reconst.freq4pT)
seq02reconst.freq4p.combined <- rowSums(seq02reconst.freq4pT) / sum(seq02reconst.freq4pT)
seq01control.freq4p.combined <- rowSums(seq01control.freq4pT) / sum(seq01control.freq4pT)
seq02control.freq4p.combined <- rowSums(seq02control.freq4pT) / sum(seq02control.freq4pT)



freq4DF <- data.frame(motif = names(seq01reconst.freq4p.combined), 
			seq01r.count = seq01reconst.freq4c.combined, 
			seq02r.count = seq02reconst.freq4c.combined, 
			seq01c.count = seq01control.freq4c.combined, 
			seq02c.count = seq02control.freq4c.combined, 
			seq01r.prob = seq01reconst.freq4p.combined, 
			seq02r.prob = seq02reconst.freq4p.combined, 
			seq01c.prob = seq01control.freq4p.combined, 
			seq02c.prob = seq02control.freq4p.combined)

ppois(q = freq4DF$seq01r.count[i], lambda = freq4DF$seq01c.count[i], lower.tail = FALSE)

freq4DF$seq01.p.value <- unlist(Map(f = ppois, q = freq4DF$seq01r.count, lambda = freq4DF$seq01c.count, lower.tail = FALSE))
freq4DF$seq02.p.value <- unlist(Map(f = ppois, q = freq4DF$seq02r.count, lambda = freq4DF$seq02c.count, lower.tail = FALSE))

freq4DF$seq01.p.value2 <- unlist(Map(f = ppois, q = freq4DF$seq01r.count, lambda = mean(freq4DF$seq01c.count), lower.tail = FALSE))
freq4DF$seq02.p.value2 <- unlist(Map(f = ppois, q = freq4DF$seq02r.count, lambda = mean(freq4DF$seq02c.count), lower.tail = FALSE))


freq4DF.ordered01 <- freq4DF[order(freq4DF$seq01r.prob, decreasing = TRUE), ]
freq4DF.ordered02 <- freq4DF[order(freq4DF$seq02r.prob, decreasing = TRUE), ]

temp <- freq4DF$seq01.p.value
temp <- temp[temp > 0]
summary(temp)
min(temp)


setwd(HOMEDIR)
setwd("data")
setwd("SELEX")
setwd("onePercent")
save(freq4DF, file = "20250419_freq4DF.RData")
library(openxlsx)
write.xlsx(freq4DF, file = "20250419_freq4DF.xlsx")
write.xlsx(freq4DF.ordered01, file = "20250419_freq4DF_ordered01.xlsx")
write.xlsx(freq4DF.ordered02, file = "20250419_freq4DF_ordered02.xlsx")




