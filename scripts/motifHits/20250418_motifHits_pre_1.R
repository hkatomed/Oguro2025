## original note: 20250418_motifHits_pre_1.R.txt

HOMEDIR <- ""

setwd(HOMEDIR)
setwd("data")
setwd("motifHits")
setwd("XLSX")

pone <- read.xlsx(xlsxFile = "pone.0017272.s004_column_I_P_removed2.xlsx")
RPKMtable <- read.xlsx(xlsxFile = "Yeast_CHX_KSU_RPKM_table_2_hk2.xlsx")

pone$RPKMtable <- FALSE
geneNames <- unique(RPKMtable$Gene)
for(i in 1:length(geneNames)){
	pone$RPKMtable[pone$ORF.name == geneNames[i]] <- TRUE
}
poneSpecific <- pone$ORF.name[pone$RPKMtable == FALSE]
poneSpecificDF <- subset(pone,RPKMtable == FALSE)
write.xlsx(poneSpecificDF, file = "20250418_poneSpecific.xlsx")

common <- pone$ORF.name[pone$RPKMtable == TRUE]
commonDF <- subset(pone,RPKMtable == TRUE)

RPKMtable$pone <- FALSE
geneNames <- common
for(i in 1:length(geneNames)){
	RPKMtable$pone[RPKMtable$Gene == geneNames[i]] <- TRUE
}

RPKMtableSpecific <- RPKMtable$Gene[RPKMtable$pone == FALSE]
RPKMtableSpecificDF <- subset(RPKMtable, pone == FALSE)
write.xlsx(RPKMtableSpecificDF, file = "20250418_RPKMtableSpecific.xlsx")

NewDF <- RPKMtable[,c(1, ncol(RPKMtable))]
NewDF$order <- 1:nrow(NewDF)

colnames(commonDF)[1] <- "Gene"
NewDF <- merge(NewDF, commonDF, by = "Gene", all = TRUE)
NewDF <- NewDF[order(NewDF$order),]

NewDF$order <- NULL
NewDF$RPKMtable <- NULL



S4 <- read.xlsx(xlsxFile = "Table S4_246 genes.xlsx")
S4genes <- S4$Name
NewDF$S4_246 <- as.numeric(NA)
for(i in 1:length(S4genes)){
	NewDF$S4_246[NewDF$Gene == S4genes[i]] <- 1
}
write.xlsx(NewDF, file = "20250418_RPKMtable_merged_with_pone_match4496_with444.xlsx")


setwd(HOMEDIR)
setwd("data")
setwd("motifHits")
setwd("XLSX")
RPKMtable <- read.xlsx(xlsxFile = "Yeast_CHX_KSU_RPKM_table_2_hk2.xlsx")

NewDF <- read.xlsx(xlsxFile = "20250418_RPKMtable_merged_with_pone_match4496_with444.xlsx")

NewDF <- cbind(RPKMtable[, 1:12], NewDF[, 2:22])

NewDF$tif34_down_rank <- order(NewDF$tif34.down.fold, decreasing = TRUE)
NewDF$tif35_up_rank <- order(NewDF$tif35.up.fold, decreasing = TRUE)

subset(NewDF, S4_246 == 1)$tif34_down_rank
order(subset(NewDF, S4_246 == 1)$tif34_down_rank)

NewDF$S4_246_tif34down_rank <- as.integer(NA)
NewDF$S4_246_tif34down_rank[which(NewDF$S4_246 == 1)] <- order(subset(NewDF, S4_246 == 1)$tif34_down_rank)

NewDF$tif34down_top50 <- FALSE
NewDF$tif34down_top50[NewDF$S4_246_tif34down_rank <= 50] <- TRUE

NewDF$tif34down_bottom50 <- FALSE
NewDF$tif34down_bottom50[NewDF$S4_246_tif34down_rank >= (246-49)] <- TRUE

NewDF$tif34down_top50[which(NewDF$S4_246 == 1)]
NewDF$tif34down_bottom50[which(NewDF$S4_246 == 1)]

NewDF$RPS <- FALSE
NewDF$RPL <- FALSE
NewDF$RPS[grep(pattern = "RPS", x = NewDF$Gene.name)] <- TRUE
NewDF$RPL[grep(pattern = "RPL", x = NewDF$Gene.name)] <- TRUE

NewDF$TIF34_TE <- NewDF$RPKM.C / NewDF$RA.t40
NewDF$txn_fold <- NewDF$RA.t40 / NewDF$RA.t0


All <- subset(NewDF, is.na(RA.t40) == FALSE & is.na(RA.t0) == FALSE & is.na(TIF34_TE) == FALSE)
Only246 <- subset(All, S4_246 == 1)
No246 <- subset(All, is.na(S4_246) == TRUE)
Top246 <- subset(All, tif34down_top50 == TRUE)
Bottom246 <- subset(All, tif34down_bottom50 == TRUE)
Middle246 <- subset(All, S4_246 == 1 & tif34down_top50 == FALSE & tif34down_bottom50 == FALSE)
NoRPgenes <- subset(All, RPS == FALSE & RPL == FALSE)
RPSgenes <- subset(All, RPS == TRUE)
RPLgenes <- subset(All, RPL == TRUE)

setwd(HOMEDIR)
setwd("data")
setwd("motifHits")
setwd("XLSX")
write.xlsx(list(No246 = No246, Top246 = Top246, Middle246 = Middle246, Bottom246 = Bottom246, All = All), file = "20250418_S4_246_genes.xlsx")
write.xlsx(list(NoRPgenes = NoRPgenes, RPSgenes = RPSgenes, RPLgenes = RPLgenes, All = All), file = "20250418_RP_genes.xlsx")

setwd("../RData")
save(NewDF, file = "20250418_NewDF.RData")


