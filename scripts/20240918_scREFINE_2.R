## This file is base on 20170117_scREFINE_templates_1.txt
## REFINE was obtained from Riordan et al. NAR 2011 (PMID: 20959291), and modifed to fit our mac.


########################
# summarize data in R and output as xls
library(WriteXLS)

output.xls <- function(filename, out.name, motif.number = 1, prefix = prefix){
	cat("output.xls\n")
	require(WriteXLS)
	motif.number <- as.character(motif.number)
	# sites から始まるファイルがあるかどうかをチェックする。
	if(file.exists(filename) == TRUE){
	motif <- read.table(file = filename, stringsAsFactors = FALSE)
	colnames(motif) <- c("gene_id", "score", "start", "end", "seq", get("prefix"))	## prefix
	filename2 <- paste(out.name, "_motif", motif.number, ".txt", sep = "")
	write(motif[motif[,6] == 1,]$seq, file = filename2)
	filename3 <- paste(out.name, "_motif", motif.number, "_BG.txt", sep = "")
	write(motif$seq, file = filename3)
	filename4 <- paste(out.name, "_motif", motif.number, ".xls", sep = "")
	SheetName <- paste("motif", motif.number, sep = "")
	WriteXLS(motif[motif[,6] == 1,], ExcelFileName = filename4, 
		SheetNames = SheetName, row.names = FALSE, col.names = TRUE)
	if(nrow(motif) < 50001){
		filename5 <- paste(out.name, "_motif", motif.number, "_BG.xls", sep = "")
		SheetName2 <- paste("motif", motif.number, sep = "")
		WriteXLS(motif, ExcelFileName = filename5, 
			SheetNames = SheetName2, row.names = FALSE, col.names = TRUE)
	}
	if(nrow(motif) >= 50001){
		ofn <- ceiling(nrow(motif)/50000)	# out file number
		for(n in 1:ofn){
		motif.out <- motif[((n-1)*50000+1):(n*50000),]
		filename5 <- paste(out.name, "_motif", motif.number, "_", n, "_BG.xls", sep = "")
		SheetName2 <- paste("motif", motif.number, "_", n, sep = "")
		WriteXLS(motif.out, ExcelFileName = filename5, 
			SheetNames = SheetName2, row.names = FALSE, col.names = TRUE)
		}
	}
	}
}

output.xls2 <- function(TYPE = "UTR5", genome = "sacCer2", out.xls.name = "test"){
	curr.dir <- getwd()
	files <- list.files(pattern = "list.")
	list.number <- length(files)
	prefixes <- files
	for(i in 1:list.number)	prefixes[i] <- strsplit(files[i], split = "list.")[[1]][2]
	Kmers = c("K5", "K6", "K7")
	Markovs = c("0th", "5th")
	
		n <- 0
		df.start.pos <- data.frame(meme_out_name = character(n), Motif = character(n), 
			width = numeric(n), sites = numeric(n), Evalue = numeric(n), sequence = numeric(n), 
			Min = numeric(n), firstQ = numeric(n), 
			Median = numeric(n), Mean = numeric(n), thirdQ = numeric(n), Max = numeric(n), 
			stringsAsFactors = FALSE)

	for(i in 1:list.number) for(j in 1:length(Kmers)) for(k in 1:length(Markovs)){
		prefix <- prefixes[i]
		dir.name <- paste("output", prefix, TYPE, genome, Kmers[j], Markovs[k], sep = "-")
		cat(dir.name, "\n")
		if(file.exists(dir.name) == TRUE){
		setwd(dir.name)
		sites <- list.files(pattern = "sites.")
		n <- length(sites)
		cat(n, "\n")
		if(n != 0){
		for(l in 1:n){
			cat(sites[l], "\n")
			out.name <- paste(prefix, TYPE, genome, Kmers[j], Markovs[k], sep = "-")
			cat(out.name, "\n")
			output.xls(filename = sites[l], out.name, motif.number = l, prefix = prefix)
			
		}
		setwd(curr.dir)
		out.name <- paste(prefix, TYPE, genome, Kmers[j], Markovs[k], sep = "-")
		meme.out.name <- paste("meme_out", out.name, sep = "-")
		cat("meme.out.name:", meme.out.name, "\n")
		df.start.pos <- rbind(df.start.pos, summary.meme(dir.name = meme.out.name))
		}}
		setwd(curr.dir)
	}
	setwd(curr.dir)
	require(WriteXLS)
	filename <- out.xls.name
	SheetName <- "meme_out_summary"
	WriteXLS(df.start.pos, ExcelFileName = filename, 
		SheetNames = SheetName, row.names = FALSE, col.names = TRUE)
	# return(df.start.pos)
}

summary.meme <- function(dir.name){
	cat("summary.meme\n")
	curr.dir2 <- getwd()
	cat("curr.dir2:", curr.dir2, "\n")
	cat("dir.name:", dir.name, "\n")
	setwd(dir.name)
	memetxt <- readLines("meme.txt")
	StartIndex <- grep(pattern = "Start", memetxt) + 2
	n <- length(StartIndex)
	lineIndex <- grep(pattern = "------------------------------------", memetxt)
	EndIndex <- StartIndex
	for(i in 1:length(StartIndex)){
		EndIndex[i] <- min(lineIndex[lineIndex > StartIndex[i]]) - 1
	}
	
	
	df.start.pos <- data.frame(meme_out_name = character(n), Motif = character(n), 
			width = numeric(n), sites = numeric(n), Evalue = numeric(n), sequence = numeric(n), 
			Min = numeric(n), firstQ = numeric(n), 
			Median = numeric(n), Mean = numeric(n), thirdQ = numeric(n), Max = numeric(n), 
			stringsAsFactors = FALSE)

	for(i in 1:length(StartIndex)){
		TABLE <- memetxt[StartIndex[i]:EndIndex[i]]
		start.pos <- integer(length(TABLE))
		for(j in 1:length(TABLE)){
			start.pos[j] <- as.integer(strsplit(TABLE[j], split = "[ ]+")[[1]][2])
		}
		df.start.pos[i, 7:12] <- summary(start.pos)
		df.start.pos$meme_out_name[i] <- dir.name
		df.start.pos$Motif[i] <- as.character(i)
		PATTERN <- paste("MOTIF  ", i, sep = "")
		motif.data <- memetxt[grep(pattern = PATTERN, memetxt)]
		df.start.pos$width[i] <- as.numeric(strsplit(motif.data, split = "[ ]+")[[1]][5])
		df.start.pos$sites[i] <- as.numeric(strsplit(motif.data, split = "[ ]+")[[1]][8])
		df.start.pos$Evalue[i] <- as.numeric(strsplit(motif.data, split = "[ ]+")[[1]][14])
		MULTILEVEL <- memetxt[grep(pattern = "Multilevel", memetxt)]
		df.start.pos$sequence[i] <- strsplit(MULTILEVEL[i], split = "[ ]+")[[1]][2]
	}
	setwd(curr.dir2)
	return(df.start.pos)
}


HOMEDIR <- "xxx/yyy/REFINE"	# define home directory

## 5'-UTR ##
setwd(HOMEDIR)
setwd("20170117_scREFINE/REFINE_20170117_UTR5")
output.xls2(TYPE = "UTR5", genome = "sacCer2", out.xls.name = "UTR5_sacCer2_summary.xls")

## 3'-UTR ##
setwd(HOMEDIR)
setwd("20170117_scREFINE/REFINE_20170117_UTR3")
output.xls2(TYPE = "UTR3", genome = "sacCer2", out.xls.name = "UTR3_sacCer2_summary.xls")

## CDS ##
setwd(HOMEDIR)
setwd("20170117_scREFINE/REFINE_20170117_CDS")
output.xls2(TYPE = "CDS", genome = "sacCer3", out.xls.name = "CDS_sacCer2_summary.xls")


