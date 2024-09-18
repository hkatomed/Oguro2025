## This file is based on 20240615_Asano_motif_1_bugFixed20240620.txt


## The files below were downloaded from SGD and placed in SGD_20240615
## http://sgd-archive.yeastgenome.org/sequence/S288C_reference/

# orf_dna/orf_coding.fasta.gz		2024-06-13T22:17:19.000Z
# orf_dna/orf_dna.README		2019-10-25T18:56:06.000Z

# SGD_all_ORFs_5prime_UTRs.fsa.zip	2021-04-13T21:06:24.000Z
# SGD_all_ORFs_5prime_UTRs.README	2021-04-13T21:06:24.000Z
# SGD_all_ORFs_3prime_UTRs.fsa.zip	2021-04-13T21:06:25.000Z
# SGD_all_ORFs_3prime_UTRs.README	2021-04-13T21:06:25.000Z

# genome_releases/S288C_reference_genome_R64-5-1_20240529.tgz

# rna_coding.fasta.gz	2024-06-13T22:33:40.000Z
# rna.README		2019-10-25T18:56:13.000Z

# NotFeature.fasta.gz	2024-06-13T22:54:43.000Z
# intergenic.README	2019-10-25T18:55:53.000Z


###################################################
library(Biostrings)
HOMEDIR <- "xxx/yyy/motifDistribution"	# define home directory
setwd(HOMEDIR)
setwd("SGD_20240615")
RNA <- readDNAStringSet(filepath = "rna_coding.fasta.gz")

# > RNA
# DNAStringSet object of length 430:
#       width seq                                     names               
#   [1]   564 GGGCCCTTTCTTCCGTTT...TGGGTTTTCTGGCAAAAA YNCA0001W HRA1 SG...
#   [2]    72 GGGCGTGTGGTCTAGTGG...CAATTCCCAGCTCGCCCC YNCA0002W TRN1 SG...
#   [3]   102 GTGAATGATGAATTTAAT...CTATCACAGTATCTGACG YNCA0003W SNR18 S...
#   [4]    73 GGGCACATGGCGCAGTTG...CGATTCCGGTTGCGTCCA YNCA0004W TGA1 SG...
#   [5]    82 GGTTGTTTGGCCGAGCGG...CGAATCTCTTAGCAACCA YNCA0005W SUP56 S...
#   ...   ... ...
# [426]    75 GCTTTTATAGCTTAGTGG...TTCTCATTAAGGGCAATA YNCQ0023W YNCQ002...
# [427]    74 GTAAATATAATTTAATGG...AATCTTAGTATTTACAAA YNCQ0024C YNCQ002...
# [428]    76 AGGAGATTAGCTTAATTG...ACCCTATATTTCCTAAAT YNCQ0025W YNCQ002...
# [429]    78 TGCAATATGATGTAATTG...CGTATTATTGCTAATAAA YNCQ0026W YNCQ002...
# [430]   483 AGATATTTATATTATTAA...ATTATAATATAATATCCA YNCQ0027W RPM1 SG...

strsplit(temp, split = ", ")[[1]]

# > strsplit(temp, split = ", ")[[1]]
# [1] "YNCM0023C YNCM0023C SGDID:S000006764" "Chr XIII from 420661-420588"         
# [3] "Genome Release 64-5-1"                "reverse complement"                  
# [5] "tRNA_gene"                            "\"Valine tRNA (tRNA-Val)"            
# [7] "predicted by tRNAscan-SE analysis\"" 

info <- strsplit(temp, split = ", ")[[1]]
genomeReleaseIndex <- which(info == "Genome Release 64-5-1")
if(info[genomeReleaseIndex + 1] == "reverse complement"){
	ncRNAtype <- info[genomeReleaseIndex + 2]
}else{
	ncRNAtype <- info[genomeReleaseIndex + 1]
}

ncRNAtypes <- character(length = length(RNA))


for(i in 1:length(RNA)){
	temp <- names(RNA)[i]
	info <- strsplit(temp, split = ", ")[[1]]
	genomeReleaseIndex <- which(info == "Genome Release 64-5-1")
	if(info[genomeReleaseIndex + 1] == "reverse complement"){
		ncRNAtypes[i] <- info[genomeReleaseIndex + 2]
	}else{
		ncRNAtypes[i] <- info[genomeReleaseIndex + 1]
	}
}

ncRNAtypes <- factor(ncRNAtypes)
summary(ncRNAtypes)

# > summary(ncRNAtypes)
#          ncRNA_gene           rRNA_gene         snoRNA_gene 
#                  31                  16                  77 
#          snRNA_gene telomerase_RNA_gene           tRNA_gene 
#                   6                   1                 299 

which(ncRNAtypes == "rRNA_gene")

# > which(ncRNAtypes == "rRNA_gene")
#  [1] 256 257 258 259 260 261 262 263 264 265 266 267 268 269 405 409


whichMito <- logical(length = length(RNA))

for(i in 1:length(RNA)){
	temp <- names(RNA)[i]
	info <- strsplit(temp, split = ", ")[[1]]
	chrName <- strsplit(info[2], split = " ")[[1]][2]
	if(chrName == "Mito") whichMito[i] <- TRUE
}

newNames <- names(RNA)
for(i in 1:length(newNames)){
	newNames[i] <- strsplit(newNames[i], split = " ")[[1]][1]
}

names(RNA) <- newNames

# > RNA
# DNAStringSet object of length 430:
#       width seq                                      names               
#   [1]   564 GGGCCCTTTCTTCCGTTTG...TGGGTTTTCTGGCAAAAA YNCA0001W
#   [2]    72 GGGCGTGTGGTCTAGTGGT...CAATTCCCAGCTCGCCCC YNCA0002W
#   [3]   102 GTGAATGATGAATTTAATT...CTATCACAGTATCTGACG YNCA0003W
#   [4]    73 GGGCACATGGCGCAGTTGG...CGATTCCGGTTGCGTCCA YNCA0004W
#   [5]    82 GGTTGTTTGGCCGAGCGGT...CGAATCTCTTAGCAACCA YNCA0005W
#   ...   ... ...
# [426]    75 GCTTTTATAGCTTAGTGGT...TTCTCATTAAGGGCAATA YNCQ0023W
# [427]    74 GTAAATATAATTTAATGGT...AATCTTAGTATTTACAAA YNCQ0024C
# [428]    76 AGGAGATTAGCTTAATTGG...ACCCTATATTTCCTAAAT YNCQ0025W
# [429]    78 TGCAATATGATGTAATTGG...CGTATTATTGCTAATAAA YNCQ0026W
# [430]   483 AGATATTTATATTATTAAT...ATTATAATATAATATCCA YNCQ0027W


ncRNA <- RNA[which(ncRNAtypes == "ncRNA_gene" & whichMito == FALSE)]
rRNA <- RNA[which(ncRNAtypes == "rRNA_gene" & whichMito == FALSE)]	
snoRNA <- RNA[which(ncRNAtypes == "snoRNA_gene" & whichMito == FALSE)]
snRNA <- RNA[which(ncRNAtypes == "snRNA_gene" & whichMito == FALSE)]
teloRNA <- RNA[which(ncRNAtypes == "telomerase_RNA_gene" & whichMito == FALSE)]
tRNA <- RNA[which(ncRNAtypes == "tRNA_gene" & whichMito == FALSE)]

length(ncRNA)
length(rRNA)
length(snoRNA)
length(snRNA)
length(teloRNA)
length(tRNA)

# > length(ncRNA)
# [1] 30
# > length(rRNA)
# [1] 14
# > length(snoRNA)
# [1] 77
# > length(snRNA)
# [1] 6
# > length(teloRNA)
# [1] 1
# > length(tRNA)
# [1] 275


writeXStringSet(ncRNA, filepath = "ncRNA.fasta.gz", compress = TRUE)
writeXStringSet(rRNA, filepath = "rRNA.fasta.gz", compress = TRUE)
writeXStringSet(snoRNA, filepath = "snoRNA.fasta.gz", compress = TRUE)
writeXStringSet(snRNA, filepath = "snRNA.fasta.gz", compress = TRUE)
writeXStringSet(teloRNA, filepath = "teloRNA.fasta.gz", compress = TRUE)
writeXStringSet(tRNA, filepath = "tRNA.fasta.gz", compress = TRUE)


############################
CDS <- readDNAStringSet(filepath = "orf_coding.fasta.gz")

# > CDS
# DNAStringSet object of length 6039:
#        width seq                                     names               
#    [1]   363 ATGGTCAAATTAACTTCA...TACACTATCGCAAACTAG YAL068C PAU8 SGDI...
#    [2]   228 ATGCCAATTATAGGGGTG...GGAGTCGTATACTGTTAG YAL067W-A YAL067W...
#    [3]  1782 ATGTATTCAATTGTTAAA...GTATCTGATGAAAAATAA YAL067C SEO1 SGDI...
#    [4]   387 ATGAACAGTGCTACCAGT...CTGGCAATCGTATGGTAA YAL065C YAL065C S...
#    [5]   381 ATGGCAGGTGAAGCAGTT...GATTCAGTGCACACATGA YAL064W-B YAL064W...
#    ...   ... ...
# [6035]  1197 ATGAAATTAAAATTATTA...AAATTAAACTTTATTTAA Q0140 VAR1 SGDID:...	# Mito
# [6036]   708 ATGAAAAATATTAAAAAA...GAAACTTTTTTAAAATAA Q0160 SCEI SGDID:...	# Mito
# [6037]   756 ATGTTAGATTTATTAAGA...TGATTAAATGAACAATAA Q0250 COX2 SGDID:...	# Mito
# [6038]  1419 ATGATTAAATGAACAATA...AATATTTTATTAGATTAA Q0255 Q0255 SGDID...	# Mito
# [6039]   810 ATGACACATTTAGAAAGA...TACTGATGAGGAGTCTAA Q0275 COX3 SGDID:...	# Mito


whichMito <- logical(length = length(CDS))

for(i in 1:length(CDS)){
	temp <- names(CDS)[i]
	info <- strsplit(temp, split = ", ")[[1]]
	chrName <- strsplit(info[2], split = " ")[[1]][2]
	if(chrName == "Mito") whichMito[i] <- TRUE
}

CDS <- CDS[whichMito == FALSE]

# > CDS
# DNAStringSet object of length 6020:
#        width seq                                     names               
#    [1]   363 ATGGTCAAATTAACTTCA...TACACTATCGCAAACTAG YAL068C PAU8 SGDI...
#    [2]   228 ATGCCAATTATAGGGGTG...GGAGTCGTATACTGTTAG YAL067W-A YAL067W...
#    [3]  1782 ATGTATTCAATTGTTAAA...GTATCTGATGAAAAATAA YAL067C SEO1 SGDI...
#    [4]   387 ATGAACAGTGCTACCAGT...CTGGCAATCGTATGGTAA YAL065C YAL065C S...
#    [5]   381 ATGGCAGGTGAAGCAGTT...GATTCAGTGCACACATGA YAL064W-B YAL064W...
#    ...   ... ...
# [6016]   393 ATGGTAAGTTTCATAACG...TTGATTGTTAGTGGTTGA YPR200C ARR2 SGDI...
# [6017]  1215 ATGTCAGAAGATCAAAAA...TGGAACAATAGAAATTAA YPR201W ARR3 SGDI...
# [6018]   717 ATGGAAATTGAAAACGAA...GATGCGTACTTTCACTGA YPR202W YPR202W S...
# [6019]   309 ATGCGTACTTTCACTGAC...GTGTGCTGCCCAAGTTGA YPR203W YPR203W S...
# [6020]  3099 ATGGCAGACACACCCTCT...AGAGAAGTTGGAGAGTGA YPR204W YPR204W S...


newNames <- names(CDS)
for(i in 1:length(newNames)){
	newNames[i] <- strsplit(newNames[i], split = " ")[[1]][1]
}

names(CDS) <- newNames

# > CDS
# DNAStringSet object of length 6039:
#        width seq                                     names               
#    [1]   363 ATGGTCAAATTAACTTCA...TACACTATCGCAAACTAG YAL068C
#    [2]   228 ATGCCAATTATAGGGGTG...GGAGTCGTATACTGTTAG YAL067W-A
#    [3]  1782 ATGTATTCAATTGTTAAA...GTATCTGATGAAAAATAA YAL067C
#    [4]   387 ATGAACAGTGCTACCAGT...CTGGCAATCGTATGGTAA YAL065C
#    [5]   381 ATGGCAGGTGAAGCAGTT...GATTCAGTGCACACATGA YAL064W-B
#    ...   ... ...
# [6016]   393 ATGGTAAGTTTCATAACG...TTGATTGTTAGTGGTTGA YPR200C
# [6017]  1215 ATGTCAGAAGATCAAAAA...TGGAACAATAGAAATTAA YPR201W
# [6018]   717 ATGGAAATTGAAAACGAA...GATGCGTACTTTCACTGA YPR202W
# [6019]   309 ATGCGTACTTTCACTGAC...GTGTGCTGCCCAAGTTGA YPR203W
# [6020]  3099 ATGGCAGACACACCCTCT...AGAGAAGTTGGAGAGTGA YPR204W

# > summary(width(CDS))
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#      51     678    1194    1464    1870   14733 

CDSwidth <- width(CDS)

# > head(CDSwidth)
# [1]  363  228 1782  387  381  381


CDS_over101 <- CDS
for(i in 1:length(CDS)){
	if(i%%100 == 0) cat(i, ", ", sep = "")
	if(CDSwidth[i] < 102){
		CDS_over101[[i]] <- DNAString("")
	}
}
names(CDS_over101) <- names(CDS)


CDS_first51 <- CDS
for(i in 1:length(CDS)){
	if(i%%100 == 0) cat(i, ", ", sep = "")
	SEQ <- CDS[[i]]
	if(CDSwidth[i] >= 102){
		CDS_first51[[i]] <- subseq(SEQ, start = 1, end = 51)
	}else{
		CDS_first51[[i]] <- DNAString("")
	}
}
names(CDS_first51) <- names(CDS)

CDS_last51 <- CDS
for(i in 1:length(CDS)){
	if(i%%100 == 0) cat(i, ", ", sep = "")
	SEQ <- CDS[[i]]
	if(CDSwidth[i] >= 102){
		CDS_last51[[i]] <- subseq(SEQ, start = CDSwidth[i] - 50, end = CDSwidth[i])
	}else{
		CDS_last51[[i]] <- DNAString("")
	}
}
names(CDS_last51) <- names(CDS)


CDS_52toLast <- CDS
for(i in 1:length(CDS)){
	if(i%%100 == 0) cat(i, ", ", sep = "")
	SEQ <- CDS[[i]]
	# if(CDSwidth[i] >= 1012){		# バグに気づいた（20240620）
	if(CDSwidth[i] >= 102){
		CDS_52toLast[[i]] <- subseq(SEQ, start = 52, end = CDSwidth[i])
	}else{
		CDS_52toLast[[i]] <- DNAString("")
	}
}
names(CDS_52toLast) <- names(CDS)


CDS_1toLast52 <- CDS
for(i in 1:length(CDS)){
	if(i%%100 == 0) cat(i, ", ", sep = "")
	SEQ <- CDS[[i]]
	if(CDSwidth[i] >= 102){
		CDS_1toLast52[[i]] <- subseq(SEQ, start = 1, end = CDSwidth[i] - 51)
	}else{
		CDS_1toLast52[[i]] <- DNAString("")
	}
}
names(CDS_1toLast52) <- names(CDS)

length(CDS_over101)
length(CDS_first51)
length(CDS_last51)
length(CDS_52toLast)
length(CDS_1toLast52)

# > length(CDS_over101)
# [1] 6020
# > length(CDS_first51)
# [1] 6020
# > length(CDS_last51)
# [1] 6020
# > length(CDS_52toLast)
# [1] 6020
# > length(CDS_1toLast52)
# [1] 6020


writeXStringSet(CDS, filepath = "CDS_total.fasta.gz", compress = TRUE)
writeXStringSet(CDS_over101, filepath = "CDS_over101.fasta.gz", compress = TRUE)
writeXStringSet(CDS_first51, filepath = "CDS_first51.fasta.gz", compress = TRUE)
writeXStringSet(CDS_last51, filepath = "CDS_last51.fasta.gz", compress = TRUE)
writeXStringSet(CDS_52toLast, filepath = "CDS_52toLast.fasta.gz", compress = TRUE)
writeXStringSet(CDS_1toLast52, filepath = "CDS_1toLast52.fasta.gz", compress = TRUE)


#####################################################
command <- "unzip -K SGD_all_ORFs_5prime_UTRs.fsa.zip"
system(command)
command <- "unzip -K SGD_all_ORFs_3prime_UTRs.fsa.zip"
system(command)

UTR5 <- readDNAStringSet(filepath = "SGD_all_ORFs_5prime_UTRs.fsa")
UTR3 <- readDNAStringSet(filepath = "SGD_all_ORFs_3prime_UTRs.fsa")

# > names(UTR5)[9722: 9723]
# [1] "sacCer3_ct_PelechanoonlybasedUTRs_1122_YMR320W_id001_five_prime_UTR range=chrXIII:915443-916746 5'pad=0 3'pad=0 strand=+ repeatMasking=none"
# [2] "sacCer3_ct_PelechanoonlybasedUTRs_1122_YMR320W_id003_five_prime_UTR range=chrXIII:915706-916746 5'pad=0 3'pad=0 strand=+ repeatMasking=none"

# > names(UTR5)[9720: 9721]
# [1] "sacCer3_ct_PelechanoonlybasedUTRs_1122_YMR319C_id004_five_prime_UTR range=chrXIII:914538-914714 5'pad=0 3'pad=0 strand=- repeatMasking=none"
# [2] "sacCer3_ct_PelechanoonlybasedUTRs_1122_YMR319C_id001_five_prime_UTR range=chrXIII:914538-914657 5'pad=0 3'pad=0 strand=- repeatMasking=none"

UTR5.df <- data.frame(names(UTR5))
UTR5.df$geneName <- as.character(NA)
UTR5.df$UTRid <- as.character(NA)
UTR5.df$UTRleft <- as.integer(NA)
UTR5.df$UTRright <- as.integer(NA)
UTR5.df$length <- as.integer(NA)
UTR5.df$strand <- as.character(NA)
for(i in 1:nrow(UTR5.df)){
	if(i%%100 == 0) cat(i, ", ", sep = "")
	temp <- UTR5.df$names.UTR5[i]
	info <- strsplit(temp, split = "_")[[1]]
	UTR5.df$geneName[i] <- info[5]
	UTR5.df$UTRid[i] <- info[6]
	temp2 <- info[9]
	info2 <- strsplit(temp2, split = " ")[[1]]
	UTRrange <- strsplit(info2[2], split = ":")[[1]][2]
	UTRranges <- as.integer(strsplit(UTRrange, split = "-")[[1]])
	UTR5.df$UTRleft[i] <- UTRranges[1]
	UTR5.df$UTRright[i] <- UTRranges[2]
	UTR5.df$length[i] <- abs(UTRranges[2] - UTRranges[1]) + 1
	UTR5.df$strand[i] <- strsplit(info2[5], split = "=")[[1]][2]
}

# > str(UTR5.df)
# 'data.frame':	9723 obs. of  5 variables:
#  $ names.UTR5.: chr  "sacCer3_ct_PelechanoonlybasedUTRs_1122_YAL067C_id001_five_prime_UTR range=chrI:9016-9049 5'pad=0 3'pad=0 strand"| __truncated__ "sacCer3_ct_PelechanoonlybasedUTRs_1122_YAL066W_id001_five_prime_UTR range=chrI:9807-10091 5'pad=0 3'pad=0 stran"| __truncated__ "sacCer3_ct_PelechanoonlybasedUTRs_1122_YAL064W-B_id001_five_prime_UTR range=chrI:11313-12046 5'pad=0 3'pad=0 st"| __truncated__ "sacCer3_ct_PelechanoonlybasedUTRs_1122_YAL063C-A_id001_five_prime_UTR range=chrI:22685-23617 5'pad=0 3'pad=0 st"| __truncated__ ...
#  $ geneName   : chr  "YAL067C" "YAL066W" "YAL064W-B" "YAL063C-A" ...
#  $ UTRid      : chr  "id001" "id001" "id001" "id001" ...
#  $ UTRleft    : int  9016 9807 11313 22685 22685 31527 31529 33345 33359 34983 ...
#  $ UTRright   : int  9049 10091 12046 23617 23219 31567 31567 33448 33448 35155 ...
#  $ length     : num  34 285 734 933 535 41 39 104 90 173 ...
#  $ strand     : chr  "-" "+" "+" "-" ...

# > UTR5.df[1:30, 2:7]
#     geneName UTRid UTRleft UTRright length strand
# 1    YAL067C id001    9016     9049     34      -
# 2    YAL066W id001    9807    10091    285      +
# 3  YAL064W-B id001   11313    12046    734      +
# 4  YAL063C-A id001   22685    23617    933      -
# 5  YAL063C-A id006   22685    23219    535      -
# 6    YAL062W id001   31527    31567     41      +
# 7    YAL062W id003   31529    31567     39      +
# 8    YAL061W id003   33345    33448    104      +
# 9    YAL061W id001   33359    33448     90      +
# 10   YAL060W id001   34983    35155    173      +
# 11   YAL060W id014   35107    35155     49      +
# 12   YAL059W id001   36463    36509     47      +
# 13   YAL059W id002   36466    36509     44      +
# 14   YAL058W id001   37405    37464     60      +
# 15   YAL058W id003   37419    37464     46      +
# 16   YAL056W id001   39196    39259     64      +
# 17   YAL056W id006   39234    39259     26      +
# 18   YAL055W id003   42142    42177     36      +
# 19   YAL055W id001   42167    42177     11      +
# 20   YAL054C id001   45022    45073     52      -
# 21   YAL053W id001   45675    45899    225      +
# 22   YAL053W id003   45728    45899    172      +
# 23   YAL049C id001   52595    54082   1488      -
# 24   YAL049C id002   52595    53580    986      -
# 25 YAL047W-A id001   54209    54584    376      +
# 26 YAL047W-A id002   54221    54584    364      +
# 27   YAL047C id002   56857    56886     30      -	#
# 28   YAL047C id001   56857    56886     30      -	#
# 29 YAL044W-A id001   56955    57518    564      +
# 30 YAL044W-A id004   56973    57518    546      +


UTR5.df$max <- FALSE
for(i in 1:nrow(UTR5.df)){
	if(i%%100 == 0) cat(i, ", ", sep = "")
	geneName <- UTR5.df$geneName[i]
	UTRid <- UTR5.df$UTRid[i]
	UTRlength <- UTR5.df$length[i]
	UTRlengths <- UTR5.df$length[UTR5.df$geneName == geneName]
	if(UTRlength == max(UTRlengths)){
		if(length(which(UTRlengths == max(UTRlengths))) == 1){
			UTR5.df$max[i] <- TRUE
		}
		if(length(which(UTRlengths == max(UTRlengths))) > 1){
			ids <- UTR5.df$UTRid[UTR5.df$geneName == geneName]
			idNumbers <- as.integer(substr(ids, start = 3, stop = 5))
			if(UTRid == ids[which.min(idNumbers)]){
				UTR5.df$max[i] <- TRUE
			}
		}
	}
}

# > UTR5.df[1:30, 2:8]
#     geneName UTRid UTRleft UTRright length strand   max
# 1    YAL067C id001    9016     9049     34      -  TRUE
# 2    YAL066W id001    9807    10091    285      +  TRUE
# 3  YAL064W-B id001   11313    12046    734      +  TRUE
# 4  YAL063C-A id001   22685    23617    933      -  TRUE
# 5  YAL063C-A id006   22685    23219    535      - FALSE
# 6    YAL062W id001   31527    31567     41      +  TRUE
# 7    YAL062W id003   31529    31567     39      + FALSE
# 8    YAL061W id003   33345    33448    104      +  TRUE
# 9    YAL061W id001   33359    33448     90      + FALSE
# 10   YAL060W id001   34983    35155    173      +  TRUE
# 11   YAL060W id014   35107    35155     49      + FALSE
# 12   YAL059W id001   36463    36509     47      +  TRUE
# 13   YAL059W id002   36466    36509     44      + FALSE
# 14   YAL058W id001   37405    37464     60      +  TRUE
# 15   YAL058W id003   37419    37464     46      + FALSE
# 16   YAL056W id001   39196    39259     64      +  TRUE
# 17   YAL056W id006   39234    39259     26      + FALSE
# 18   YAL055W id003   42142    42177     36      +  TRUE
# 19   YAL055W id001   42167    42177     11      + FALSE
# 20   YAL054C id001   45022    45073     52      -  TRUE
# 21   YAL053W id001   45675    45899    225      +  TRUE
# 22   YAL053W id003   45728    45899    172      + FALSE
# 23   YAL049C id001   52595    54082   1488      -  TRUE
# 24   YAL049C id002   52595    53580    986      - FALSE
# 25 YAL047W-A id001   54209    54584    376      +  TRUE
# 26 YAL047W-A id002   54221    54584    364      + FALSE
# 27   YAL047C id002   56857    56886     30      - FALSE	# 
# 28   YAL047C id001   56857    56886     30      -  TRUE	# 
# 29 YAL044W-A id001   56955    57518    564      +  TRUE
# 30 YAL044W-A id004   56973    57518    546      + FALSE

names(UTR5) <- UTR5.df$geneName
UTR5.selected <- UTR5[UTR5.df$max == TRUE]

# > UTR5.selected
# DNAStringSet object of length 5211:
#        width seq                                     names               
#    [1]    34 CAAAGAGAACTACTGCATATATAAATAACATACA      YAL067C
#    [2]   285 GTCAAACCAAATGGTTTT...ATTTCAAAAGCACCATCA YAL066W
#    [3]   734 AAAACACAGATACCTCGA...CCGTGGCCGTCGAAACAA YAL064W-B
#    [4]   933 AACGAGGAAAAATGCCCT...GAGGTCCCACTGTATGTA YAL063C-A
#    [5]    41 GTAGCAACAGTCACCGAA...AAGGTAAAAAGTAAAAAA YAL062W
#    ...   ... ...
# [5207]   343 AGGGTTCGATAACGAAGA...TCAAGAAAAATCTGATTA YMR316W
# [5208]   924 GAGCTTCTGAACTCACTG...AAGAAAGAGAAGAAAAAA YMR316C-B
# [5209]   181 GACAATATAGGAGAAAAA...CAAGAAAAGCCAAAATCA YMR318C
# [5210]   177 GATGCTCTTCATTTAGAA...CATTCGTTAATCAAATTA YMR319C
# [5211]  1304 AAAAAAGAAAAATGGGAC...CACGAAAAGCTCTAAAAA YMR320W



#####################
UTR3.df <- data.frame(names(UTR3))
UTR3.df$geneName <- as.character(NA)
UTR3.df$UTRid <- as.character(NA)
UTR3.df$UTRleft <- as.integer(NA)
UTR3.df$UTRright <- as.integer(NA)
UTR3.df$length <- as.integer(NA)
UTR3.df$strand <- as.character(NA)
for(i in 1:nrow(UTR3.df)){
	if(i%%100 == 0) cat(i, ", ", sep = "")
	temp <- UTR3.df$names.UTR3[i]
	info <- strsplit(temp, split = "_")[[1]]
	UTR3.df$geneName[i] <- info[5]
	UTR3.df$UTRid[i] <- info[6]
	temp2 <- info[9]
	info2 <- strsplit(temp2, split = " ")[[1]]
	UTRrange <- strsplit(info2[2], split = ":")[[1]][2]
	UTRranges <- as.integer(strsplit(UTRrange, split = "-")[[1]])
	UTR3.df$UTRleft[i] <- UTRranges[1]
	UTR3.df$UTRright[i] <- UTRranges[2]
	UTR3.df$length[i] <- abs(UTRranges[2] - UTRranges[1]) + 1
	UTR3.df$strand[i] <- strsplit(info2[5], split = "=")[[1]][2]
}

# > str(UTR3.df)
# 'data.frame':	9723 obs. of  7 variables:
#  $ names.UTR3.: chr  "sacCer3_ct_Pelechanoonlybased3primeUTRs_3950_YAL067C_id001_three_prime_UTR range=chrI:7013-7235 5'pad=0 3'pad=0"| __truncated__ "sacCer3_ct_Pelechanoonlybased3primeUTRs_3950_YAL066W_id001_three_prime_UTR range=chrI:10399-10460 5'pad=0 3'pad"| __truncated__ "sacCer3_ct_Pelechanoonlybased3primeUTRs_3950_YAL064W-B_id001_three_prime_UTR range=chrI:12426-15151 5'pad=0 3'p"| __truncated__ "sacCer3_ct_Pelechanoonlybased3primeUTRs_3950_YAL063C-A_id001_three_prime_UTR range=chrI:22051-22395 5'pad=0 3'p"| __truncated__ ...
#  $ geneName   : chr  "YAL067C" "YAL066W" "YAL064W-B" "YAL063C-A" ...
#  $ UTRid      : chr  "id001" "id001" "id001" "id001" ...
#  $ UTRleft    : int  7013 10399 12426 22051 22363 32940 32940 34701 34701 36303 ...
#  $ UTRright   : int  7235 10460 15151 22395 22395 33039 33028 35113 35018 36433 ...
#  $ length     : num  223 62 2726 345 33 ...
#  $ strand     : chr  "-" "+" "+" "-" ...

# > UTR3.df[1:24, 2:7]
#     geneName UTRid UTRleft UTRright length strand
# 1    YAL067C id001    7013     7235    223      -
# 2    YAL066W id001   10399    10460     62      +
# 3  YAL064W-B id001   12426    15151   2726      +
# 4  YAL063C-A id001   22051    22395    345      -
# 5  YAL063C-A id006   22363    22395     33      -
# 6    YAL062W id001   32940    33039    100      +
# 7    YAL062W id003   32940    33028     89      +
# 8    YAL061W id001   34701    35113    413      +
# 9    YAL061W id003   34701    35018    318      +
# 10   YAL060W id014   36303    36433    131      +
# 11   YAL060W id001   36303    36415    113      +
# 12   YAL059W id002   37147    37344    198      +	#
# 13   YAL059W id001   37147    37344    198      +	#
# 14   YAL058W id003   38972    39163    192      +	#
# 15   YAL058W id001   38972    39163    192      +	#
# 16   YAL056W id006   41901    42040    140      +	#
# 17   YAL056W id001   41901    42040    140      +	#
# 18   YAL054C id001   42713    42881    169      -
# 19   YAL055W id001   42719    42884    166      +
# 20   YAL055W id003   42719    42838    120      +
# 21   YAL053W id001   48250    48312     63      +
# 22   YAL053W id003   48250    48310     61      +
# 23   YAL049C id001   51691    51855    165      -
# 24   YAL049C id002   51719    51855    137      -


UTR3.df$max <- FALSE
for(i in 1:nrow(UTR3.df)){
	if(i%%100 == 0) cat(i, ", ", sep = "")
	geneName <- UTR3.df$geneName[i]
	UTRid <- UTR3.df$UTRid[i]
	UTRlength <- UTR3.df$length[i]
	UTRlengths <- UTR3.df$length[UTR3.df$geneName == geneName]
	if(UTRlength == max(UTRlengths)){
		if(length(which(UTRlengths == max(UTRlengths))) == 1){
			UTR3.df$max[i] <- TRUE
		}
		if(length(which(UTRlengths == max(UTRlengths))) > 1){
			ids <- UTR3.df$UTRid[UTR3.df$geneName == geneName]
			idNumbers <- as.integer(substr(ids, start = 3, stop = 5))
			if(UTRid == ids[which.min(idNumbers)]){
				UTR3.df$max[i] <- TRUE
			}
		}
	}
}


# > UTR3.df[1:24, 2:8]
#     geneName UTRid UTRleft UTRright length strand   max
# 1    YAL067C id001    7013     7235    223      -  TRUE
# 2    YAL066W id001   10399    10460     62      +  TRUE
# 3  YAL064W-B id001   12426    15151   2726      +  TRUE
# 4  YAL063C-A id001   22051    22395    345      -  TRUE
# 5  YAL063C-A id006   22363    22395     33      - FALSE
# 6    YAL062W id001   32940    33039    100      +  TRUE
# 7    YAL062W id003   32940    33028     89      + FALSE
# 8    YAL061W id001   34701    35113    413      +  TRUE
# 9    YAL061W id003   34701    35018    318      + FALSE
# 10   YAL060W id014   36303    36433    131      +  TRUE
# 11   YAL060W id001   36303    36415    113      + FALSE
# 12   YAL059W id002   37147    37344    198      + FALSE	#
# 13   YAL059W id001   37147    37344    198      +  TRUE	#	
# 14   YAL058W id003   38972    39163    192      + FALSE	#
# 15   YAL058W id001   38972    39163    192      +  TRUE	#
# 16   YAL056W id006   41901    42040    140      + FALSE	#
# 17   YAL056W id001   41901    42040    140      +  TRUE	#
# 18   YAL054C id001   42713    42881    169      -  TRUE
# 19   YAL055W id001   42719    42884    166      +  TRUE
# 20   YAL055W id003   42719    42838    120      + FALSE
# 21   YAL053W id001   48250    48312     63      +  TRUE
# 22   YAL053W id003   48250    48310     61      + FALSE
# 23   YAL049C id001   51691    51855    165      -  TRUE
# 24   YAL049C id002   51719    51855    137      - FALSE

names(UTR3) <- UTR3.df$geneName
UTR3.selected <- UTR3[UTR3.df$max == TRUE]

# > UTR3.selected
# DNAStringSet object of length 5211:
#        width seq                                     names               
#    [1]   223 ATACGAGAATAATTTCTC...TCGATAAAGGCTCATCCG YAL067C
#    [2]    62 GAAAATCACAGTACAAAA...GCCTGATATATGTAAGAG YAL066W
#    [3]  2726 AGATTCGCAGGATGTCAC...AAAATCATTGTTTTATCC YAL064W-B
#    [4]   345 ATATGCTGTCCTTAATTT...TAAACCACTTTTATGTTC YAL063C-A
#    [5]   100 GCCGTAAGCGCTATTTTC...GGTTCAGCATTACGGAAG YAL062W
#    ...   ... ...
# [5207]   134 ATTACCGTTTTACGATGT...AACCGTTGAAAAATTTCT YMR315W
# [5208]    60 ACAAAAGTTTGCGTATAA...TAATTTTTAACTTATCAC YMR316W
# [5209]   142 ACAGAATCATTGACAAAT...AAGAAAACATGTTTTCAG YMR316C-B
# [5210]   299 GGTTGTCAAGCTCTTGAT...TACCAAAGCAATAGCGCC YMR318C
# [5211]   351 GCTTCATTGAACATTTTA...ATCACAGCCACTCTCGTC YMR319C

######################################
writeXStringSet(UTR5.selected, filepath = "UTR5.fasta.gz", compress = TRUE)
writeXStringSet(UTR3.selected, filepath = "UTR3.fasta.gz", compress = TRUE)


######################################

nrow(subset(UTR5.df, max == TRUE))
nrow(subset(UTR5.df, UTRid == "id001"))
nrow(subset(UTR5.df, max == TRUE & UTRid == "id001"))
nrow(subset(UTR5.df, max == FALSE & UTRid == "id001"))

# > nrow(subset(UTR5.df, max == TRUE))
# [1] 5211
# > nrow(subset(UTR5.df, UTRid == "id001"))
# [1] 5211
# > nrow(subset(UTR5.df, max == TRUE & UTRid == "id001"))
# [1] 4056
# > nrow(subset(UTR5.df, max == FALSE & UTRid == "id001"))
# [1] 1155

# > subset(UTR5.df, max == FALSE & UTRid == "id001")[1:6, 2:8]
#    geneName UTRid UTRleft UTRright length strand   max
# 9   YAL061W id001   33359    33448     90      + FALSE
# 19  YAL055W id001   42167    42177     11      + FALSE
# 36  YAL043C id001   61052    61068     17      - FALSE
# 38  YAL042W id001   61269    61316     48      + FALSE
# 40  YAL041W id001   62808    62840     33      + FALSE
# 51  YAL036C id001   76152    76203     52      - FALSE

subset(UTR5.df, geneName == "YAL061W")[, 2:8]
subset(UTR5.df, geneName == "YAL043C")[, 2:8]

# > subset(UTR5.df, geneName == "YAL061W")[, 2:8]
#   geneName UTRid UTRleft UTRright length strand   max
# 8  YAL061W id003   33345    33448    104      +  TRUE
# 9  YAL061W id001   33359    33448     90      + FALSE

# > subset(UTR5.df, geneName == "YAL043C")[, 2:8]
#    geneName UTRid UTRleft UTRright length strand   max
# 35  YAL043C id002   61052    61083     32      -  TRUE
# 36  YAL043C id001   61052    61068     17      - FALSE


######################################

# > CDS
# DNAStringSet object of length 6039:
#        width seq                                     names               
#    [1]   363 ATGGTCAAATTAACTTCA...TACACTATCGCAAACTAG YAL068C
#    [2]   228 ATGCCAATTATAGGGGTG...GGAGTCGTATACTGTTAG YAL067W-A
#    [3]  1782 ATGTATTCAATTGTTAAA...GTATCTGATGAAAAATAA YAL067C
#    [4]   387 ATGAACAGTGCTACCAGT...CTGGCAATCGTATGGTAA YAL065C
#    [5]   381 ATGGCAGGTGAAGCAGTT...GATTCAGTGCACACATGA YAL064W-B
#    ...   ... ...
# [6016]   393 ATGGTAAGTTTCATAACG...TTGATTGTTAGTGGTTGA YPR200C
# [6017]  1215 ATGTCAGAAGATCAAAAA...TGGAACAATAGAAATTAA YPR201W
# [6018]   717 ATGGAAATTGAAAACGAA...GATGCGTACTTTCACTGA YPR202W
# [6019]   309 ATGCGTACTTTCACTGAC...GTGTGCTGCCCAAGTTGA YPR203W
# [6020]  3099 ATGGCAGACACACCCTCT...AGAGAAGTTGGAGAGTGA YPR204W

# > UTR5.selected
# DNAStringSet object of length 5211:
#        width seq                                     names               
#    [1]    34 CAAAGAGAACTACTGCATATATAAATAACATACA      YAL067C
#    [2]   285 GTCAAACCAAATGGTTTT...ATTTCAAAAGCACCATCA YAL066W
#    [3]   734 AAAACACAGATACCTCGA...CCGTGGCCGTCGAAACAA YAL064W-B
#    [4]   933 AACGAGGAAAAATGCCCT...GAGGTCCCACTGTATGTA YAL063C-A
#    [5]    41 GTAGCAACAGTCACCGAA...AAGGTAAAAAGTAAAAAA YAL062W
#    ...   ... ...
# [5207]   343 AGGGTTCGATAACGAAGA...TCAAGAAAAATCTGATTA YMR316W
# [5208]   924 GAGCTTCTGAACTCACTG...AAGAAAGAGAAGAAAAAA YMR316C-B
# [5209]   181 GACAATATAGGAGAAAAA...CAAGAAAAGCCAAAATCA YMR318C
# [5210]   177 GATGCTCTTCATTTAGAA...CATTCGTTAATCAAATTA YMR319C
# [5211]  1304 AAAAAAGAAAAATGGGAC...CACGAAAAGCTCTAAAAA YMR320W

# > UTR3.selected
# DNAStringSet object of length 5211:
#        width seq                                     names               
#    [1]   223 ATACGAGAATAATTTCTC...TCGATAAAGGCTCATCCG YAL067C
#    [2]    62 GAAAATCACAGTACAAAA...GCCTGATATATGTAAGAG YAL066W
#    [3]  2726 AGATTCGCAGGATGTCAC...AAAATCATTGTTTTATCC YAL064W-B
#    [4]   345 ATATGCTGTCCTTAATTT...TAAACCACTTTTATGTTC YAL063C-A
#    [5]   100 GCCGTAAGCGCTATTTTC...GGTTCAGCATTACGGAAG YAL062W
#    ...   ... ...
# [5207]   134 ATTACCGTTTTACGATGT...AACCGTTGAAAAATTTCT YMR315W
# [5208]    60 ACAAAAGTTTGCGTATAA...TAATTTTTAACTTATCAC YMR316W
# [5209]   142 ACAGAATCATTGACAAAT...AAGAAAACATGTTTTCAG YMR316C-B
# [5210]   299 GGTTGTCAAGCTCTTGAT...TACCAAAGCAATAGCGCC YMR318C
# [5211]   351 GCTTCATTGAACATTTTA...ATCACAGCCACTCTCGTC YMR319C

transcripts <- CDS
for(i in 1:length(CDS)){
	if(i%%100 == 0) cat(i, ", ", sep = "")
	geneName <- names(CDS)[i]
	transcriptSeq <- CDS[names(CDS) == geneName][[1]]
	if(length(which(names(UTR5.selected) == geneName)) == 1){
		UTR5seq <- UTR5.selected[names(UTR5.selected) == geneName][[1]]
		transcriptSeq <- c(UTR5seq, transcriptSeq)
	}
	if(length(which(names(UTR5.selected) == geneName)) == 1){
		UTR3seq <- UTR3.selected[names(UTR3.selected) == geneName][[1]]
		transcriptSeq <- c(transcriptSeq, UTR3seq)
	}
	if(length(which(names(UTR5.selected) == geneName)) > 1){
		cat(i, "length(which(names(UTR5.selected) == geneName)) > 1", "\n", sep = "")
	}
	if(length(which(names(UTR3.selected) == geneName)) > 1){
		cat(i, "length(which(names(UTR3.selected) == geneName)) > 1", "\n", sep = "")
	}
	transcripts[[i]] <- transcriptSeq
}
names(transcripts) <- names(CDS)

# > transcripts
# DNAStringSet object of length 6039:
#        width seq                                     names               
#    [1]   363 ATGGTCAAATTAACTTCA...TACACTATCGCAAACTAG YAL068C
#    [2]   228 ATGCCAATTATAGGGGTG...GGAGTCGTATACTGTTAG YAL067W-A
#    [3]  2039 CAAAGAGAACTACTGCAT...TCGATAAAGGCTCATCCG YAL067C
#    [4]   387 ATGAACAGTGCTACCAGT...CTGGCAATCGTATGGTAA YAL065C
#    [5]  3841 AAAACACAGATACCTCGA...AAAATCATTGTTTTATCC YAL064W-B
#    ...   ... ...
# [6016]  1989 ATAAGTGGTATTATTCAT...GCAAAACTTTTTAATTCC YPR200C
# [6017]  1816 ACGCTTGCTGGATTGTCA...TATTGCGTATTGTATGTC YPR201W
# [6018]   717 ATGGAAATTGAAAACGAA...GATGCGTACTTTCACTGA YPR202W
# [6019]   309 ATGCGTACTTTCACTGAC...GTGTGCTGCCCAAGTTGA YPR203W
# [6020]  3099 ATGGCAGACACACCCTCT...AGAGAAGTTGGAGAGTGA YPR204W


writeXStringSet(transcripts, filepath = "transcripts.fasta.gz", compress = TRUE)


##############################################
command <- "tar -xvf S288C_reference_genome_R64-5-1_20240529.tgz"
system(command)

setwd("S288C_reference_genome_R64-5-1_20240529")
sc_genome <- readDNAStringSet(filepath = "S288C_reference_sequence_R64-5-1_20240529.fsa.gz")

newNames <- as.character(as.roman(1:16))
newNames[17] <- "MT"

# > newNames
#  [1] "I"    "II"   "III"  "IV"   "V"    "VI"   "VII"  "VIII" "IX"  
# [10] "X"    "XI"   "XII"  "XIII" "XIV"  "XV"   "XVI"  "MT"  

names(sc_genome) <- newNames

sc_genome_woMT <- sc_genome[1:16]

# > sc_genome_woMT
# DNAStringSet object of length 16:
#        width seq                                     names               
#  [1]  230218 CCACACCACACCCACACA...GGGTGTGGTGTGTGTGGG I
#  [2]  813184 AAATAGCCCTCATGTACG...TGTGGTGTGTGGGTGTGT II
#  [3]  316620 CCCACACACCACACCCAC...GTGGGTGTGGTGTGTGTG III
#  [4] 1531933 ACACCACACCCACACCAC...GTAGTAAGTAGCTTTTGG IV
#  [5]  576874 CGTCTCCTCCAAGCCCTG...ATTTTCATTTTTTTTTTT V
#  ...     ... ...
# [12] 1078177 CACACACACACACCACCC...CGTACATGAGGGCTATTT XII
# [13]  924431 CCACACACACACCACACC...TGTGGTGTGTGTGTGGGG XIII
# [14]  784333 CCGGCTTTCTGACCGAAA...TGTGGGTGTGGTGTGGGT XIV
# [15] 1091291 ACACCACACCCACACCAC...TGTGTGGGTGTGGTGTGT XV
# [16]  948066 AAATAGCCCTCATGTACG...TTTAATTTCGGTCAGAAA XVI

setwd("../")
writeXStringSet(sc_genome_woMT, filepath = "sc_genome.fasta.gz", compress = TRUE)


##############################################
Intergenic <- readDNAStringSet(filepath = "NotFeature.fasta.gz")

# > Intergenic
# DNAStringSet object of length 6680:
#        width seq                                     names               
#    [1]  1005 CTTGTGGTAGCAACACTA...TAAAATTGGCGTTTGTCT A:802-1806, Chr I...
#    [2]   310 TGTATTTGTTTTGTTTGT...TATGCGCAGAATGTGGGA A:2170-2479, Chr ...
#    [3]  4527 GGTCTGTAAACTTGTGAA...TGAGAAATTATTCTCGTA A:2708-7234, Chr ...
#    [4]  1074 GTATGTTATTTATATATG...TATTTCAAAAGCACCATC A:9017-10090, Chr...
#    [5]  1165 AAAATCACAGTACAAAAA...TAAGATGATGCCGTGCGT A:10400-11564, Ch...
#    ...   ... ...
# [6676]   370 CTTATTAATTATATAAAA...GGAGTATATATTTTATAA Q:78163-78532, Ch...
# [6677]   604 CTAGATATAATATTATAT...ATATATAACAATAAATTT Q:78609-79212, Ch...
# [6678]  2306 GGCTATAGAATTATATAT...TATTTCAAATATATAAGT Q:80023-82328, Ch...
# [6679]  2434 AGAATAAAAAATAAAAAG...TAAGTATTAATATATAAA Q:82601-85034, Ch...
# [6680]   182 TTAATATATAATATTTAT...TTTTAATTAATTTGAATA Q:85113-85294, Ch...

# > names(Intergenic)[1:3]
# [1] "A:802-1806, Chr I from 802-1806, Genome Release 64-5-1, between TEL01L and YAL068C"     
# [2] "A:2170-2479, Chr I from 2170-2479, Genome Release 64-5-1, between YAL068C and YAL067W-A"
# [3] "A:2708-7234, Chr I from 2708-7234, Genome Release 64-5-1, between YAL067W-A and YAL067C"

newNames <- names(Intergenic)
for(i in 1:length(newNames)){
	temp <- strsplit(newNames[i], split = ", ")[[1]]
	Range <- strsplit(temp[1], split = ":")[[1]][2]
	chrName <- strsplit(temp[2], split = " ")[[1]][2]
	newNames[i] <- paste(chrName, ":", Range, sep = "")
}

names(Intergenic) <- newNames

# > Intergenic
# DNAStringSet object of length 6680:
#        width seq                                     names               
#    [1]  1005 CTTGTGGTAGCAACACTA...TAAAATTGGCGTTTGTCT I:802-1806
#    [2]   310 TGTATTTGTTTTGTTTGT...TATGCGCAGAATGTGGGA I:2170-2479
#    [3]  4527 GGTCTGTAAACTTGTGAA...TGAGAAATTATTCTCGTA I:2708-7234
#    [4]  1074 GTATGTTATTTATATATG...TATTTCAAAAGCACCATC I:9017-10090
#    [5]  1165 AAAATCACAGTACAAAAA...TAAGATGATGCCGTGCGT I:10400-11564
#    ...   ... ...
# [6676]   370 CTTATTAATTATATAAAA...GGAGTATATATTTTATAA Mito:78163-78532
# [6677]   604 CTAGATATAATATTATAT...ATATATAACAATAAATTT Mito:78609-79212
# [6678]  2306 GGCTATAGAATTATATAT...TATTTCAAATATATAAGT Mito:80023-82328
# [6679]  2434 AGAATAAAAAATAAAAAG...TAAGTATTAATATATAAA Mito:82601-85034
# [6680]   182 TTAATATATAATATTTAT...TTTTAATTAATTTGAATA Mito:85113-85294

whichMito <- logical(length = length(newNames))
for(i in 1:length(newNames)){
	if(i%%100 == 0) cat(i, ", ", sep = "")
	chrName <- strsplit(newNames[i], split = ":")[[1]][1]
	if(chrName == "Mito") whichMito[i] <- TRUE
}

nonMito <- whichMito == FALSE
Intergenic_woMT <- Intergenic[nonMito]

# > Intergenic_woMT
# DNAStringSet object of length 6635:
#        width seq                                     names               
#    [1]  1005 CTTGTGGTAGCAACACTA...TAAAATTGGCGTTTGTCT I:802-1806
#    [2]   310 TGTATTTGTTTTGTTTGT...TATGCGCAGAATGTGGGA I:2170-2479
#    [3]  4527 GGTCTGTAAACTTGTGAA...TGAGAAATTATTCTCGTA I:2708-7234
#    [4]  1074 GTATGTTATTTATATATG...TATTTCAAAAGCACCATC I:9017-10090
#    [5]  1165 AAAATCACAGTACAAAAA...TAAGATGATGCCGTGCGT I:10400-11564
#    ...   ... ...
# [6631]   493 CAATACATTGTTTTTGCA...TATCGAAAAGCGGCTTAG XVI:933405-933897
# [6632]  2482 ACAATTTTCCAACGTATA...GTGCTTACAGGAAGAATA XVI:935666-938147
# [6633]   246 TATCTTGGTCACTTATCT...CAATAAAGCTTGAGGAGC XVI:939033-939278
# [6634]   250 TACGCTTGCTGGATTGTC...CATCAGGTTAGTAGAATA XVI:939672-939921
# [6635]   947 TTGTTGACTCACCAAAAA...CAAGACCATGTATGCATA XVI:941137-942083


writeXStringSet(Intergenic_woMT, filepath = "Intergenic.fasta.gz", compress = TRUE)


