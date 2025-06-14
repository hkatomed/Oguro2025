## original note: 	20250516_EightBase_sliding9base_1base_in200base_STOP.R.txt
## 			20250516_EightBase_sliding12base_1base_in200base_STOP.R.txt
## 			20250516_EightBase_sliding15base_1base_in200base_STOP.R.txt
## 			20250516_EightBase_sliding18base_1base_in200base_STOP.R.txt
## 			20250516_EightBase_sliding21base_1base_in200base_STOP.R.txt
## 			20250516_EightBase_sliding24base_1base_in200base_STOP.R.txt
## 			20250516_EightBase_sliding27base_1base_in200base_STOP.R.txt

baseNum <- 1	
# slideLEN <- 9				
for(slideLEN in c(9, 12, 15, 18, 21, 24, 27)){		

library(openxlsx)
library(ggplot2)

HOMEDIR <- ""

setwd(HOMEDIR)
setwd("data")
setwd("motifHits")
setwd("XLSX")
EightBase_hits <- read.xlsx(xlsxFile = "20250516_EightBase_hits_200base_STOP.xlsx", colNames = TRUE)

EightBase_hits$hit[which(EightBase_hits$hit == "")] <- "0"
EightBase_hits$hit <- as.integer(EightBase_hits$hit)

setwd("../TXT")
sysNames <- readLines(con = "20250513_246_ranked_sysName.txt")
geneNames <- readLines(con = "20250513_246_ranked_geneName.txt")

RANK <- 1:246
RANK_243 <- RANK[-c(60, 62, 63)]


EightBase_hits <- cbind(data.frame(geneName = geneNames[RANK_243], sysName = sysNames[RANK_243], rank = RANK_243), 
			EightBase_hits[RANK_243,])

EightBase_hits <- EightBase_hits[, c(4, 2, 3, 5, 6)]

STARTsEND <- 200 - slideLEN + 1		
ENDsSTART <- slideLEN
STARTs <- seq(1, STARTsEND, baseNum)
ENDs <- seq(ENDsSTART, 200, baseNum)	

motifLen <- 8

for(i in 1:length(STARTs)){
	START <- STARTs[i]
	END <- ENDs[i]
	colName <- paste("r", START, "to", END, sep = "")
	EightBase_hits[[colName]] <- as.integer(NA)
	for(j in 1:nrow(EightBase_hits)){
		hits <- EightBase_hits$Seq_Match[j]
		temp <- as.integer(strsplit(hits, split = "")[[1]])
		selected_hits <- temp[START:(END-motifLen+1)]
		EightBase_hits[[colName]][j] <- sum(selected_hits, na.rm = TRUE)
	}
}

EightBase_hits_per51 <- EightBase_hits

for(i in 1:length(STARTs)){
	START <- STARTs[i]
	END <- ENDs[i]
	colName <- paste("r", START, "to", END, sep = "")
	for(j in 1:nrow(EightBase_hits_per51)){
		EightBase_hits_per51[[colName]][j] <- EightBase_hits_per51[[colName]][j] *
			51/(END - START + 1)
	}
}

out3 <- data.frame(start = STARTs -199 , end = ENDs -199, intercept = as.numeric(NA), slope = as.numeric(NA))

for(i in 1:length(STARTs)){
	START <- STARTs[i]
	END <- ENDs[i]
	colName <- paste("r", START, "to", END, sep = "")
	x = EightBase_hits_per51$rank
	y = EightBase_hits_per51[[colName]]
	temp <- lm(y ~ x)
	out3$intercept[i] <- temp$coefficients[1]
	out3$slope[i] <- temp$coefficients[2]
}

setwd(HOMEDIR)
setwd("data")
setwd("motifHits")
dir.create("PDF")
setwd("PDF")
DAY <- gsub(pattern = "-", replacement = "", x = as.character(Sys.Date()))
fileName <- paste(DAY, "_EightBase_slope_LEN", slideLEN, "_", baseNum, "base_in200base_STOP.pdf", sep = "")
pdf(file = fileName)

p <- ggplot(out3, aes(start, -slope)) + geom_smooth(method = "loess", span = 0.1)
p + geom_point()

dev.off()


DAY <- gsub(pattern = "-", replacement = "", x = as.character(Sys.Date()))
fileName2 <- paste(DAY, "_EightBase_hits_LEN", slideLEN, "_", baseNum, "base_in200base_STOP.pdf", sep = "")
# pdf(file = "20250514_EightBase_hits_LEN24_3base.pdf")
pdf(file = fileName2)

for(i in 1:length(STARTs)){
	START <- STARTs[i]
	END <- ENDs[i]
	colName <- paste("r", START, "to", END, sep = "")
	x = EightBase_hits_per51$rank
	y = EightBase_hits_per51[[colName]]
	a = out3$intercept[out3$end == END -199]
	b = out3$slope[out3$end == END -199]
	plot(x = x, y = y, main = paste("Search area, bases ", START, " to ", END, sep =""), 
		ylim = c(0, 12), cex.lab = 1.3, cex.axis = 1.2, 
		pch = 19, col = "dodgerblue3", cex = 1, 
		xlab = "Ranking of genes by footprint density ratio (tif34-1/TIF34)", 
		ylab = "Number of motifs identified within:", 
		sub = paste("intercept=", round(a, digits = 5), 
			", slope=", round(b, digits = 5), sep = ""))
	abline(a = a, b = b)
}

dev.off()


setwd(HOMEDIR)
setwd("data")
setwd("motifHits")
setwd("XLSX")
outList <- list(EightBase_hits = EightBase_hits, 
		EightBase_hits_per51 = EightBase_hits_per51, 
		bases_x_to_y = out3)
DAY <- gsub(pattern = "-", replacement = "", x = as.character(Sys.Date()))
fileName3 <- paste(DAY, "_EightBase_slope_LEN", slideLEN, "_", baseNum, "base_in200base_STOP.xlsx", sep = "")

write.xlsx(outList, file = fileName3)

}	# for-loop
