## original note: 20250614_AsanoTE_2.txt


library(openxlsx)

HOMEDIR <- ""
setwd(HOMEDIR)
setwd("data/data_plots")
sourceFileName <- "v2-Merged-Yeast_CHX_KSU_RPKM_table_2_hk2_pone2-YM-ka1.xlsx"
RPKMtable <- read.xlsx(xlsxFile = sourceFileName)

NewDF <- RPKMtable

setwd(HOMEDIR)
setwd("data/data_plots")
sourceFileName <- "Table S6.xlsx"
tHSR <- read.xlsx(xlsxFile = sourceFileName)

NewDF$tHSR <- FALSE
for(i in 1:nrow(tHSR)){
	NewDF$tHSR[NewDF$Gene == tHSR$Gene[i]] <- TRUE
}

NewDF$RPS <- FALSE
NewDF$RPL <- FALSE

setwd(HOMEDIR)
setwd("data/data_plots")
RPLgenes <- read.xlsx(xlsxFile = "20250614_RPgenes.xlsx", sheet = "RPL", colNames = FALSE)[, c(1, 15)]
RPSgenes <- read.xlsx(xlsxFile = "20250614_RPgenes.xlsx", sheet = "RPS", colNames = FALSE)[, c(1, 15)]

for(i in 1:nrow(RPLgenes)){
	NewDF$RPL[NewDF$Gene == RPLgenes[i, 1]] <- TRUE
}

for(i in 1:nrow(RPSgenes)){
	NewDF$RPS[NewDF$Gene == RPSgenes[i, 1]] <- TRUE
}

All <- subset(NewDF, is.na(RPKM.C) == FALSE & is.na(Ave.WT.FP) == FALSE & 
			is.na(Ave.WT.mRNA) == FALSE & is.na(TE.DK_WT) == FALSE & 
			is.na(TE.C.RPKM) == FALSE & is.na(NewDF[["Delta.C/DK_WT"]]) == FALSE)
tHSR <- subset(All, tHSR == TRUE)
OtherGenes <- subset(All, RPS == FALSE & RPL == FALSE & tHSR == FALSE) 

RPSgenes <- subset(All, RPS == TRUE)
RPLgenes <- subset(All, RPL == TRUE)

setwd(HOMEDIR)
setwd("data/data_plots")
save(list = c("NewDF", "All", "tHSR", 
		"OtherGenes", "RPSgenes", "RPLgenes"), file = "20250614_NewDFetc.RData")

## tHSR RPS RPL
## 
xmin <- log10(min(OtherGenes$Ave.WT.FP, tHSR$Ave.WT.FP, RPSgenes$Ave.WT.FP, RPLgenes$Ave.WT.FP))
xmax <- log10(max(OtherGenes$Ave.WT.FP, tHSR$Ave.WT.FP, RPSgenes$Ave.WT.FP, RPLgenes$Ave.WT.FP))
xlim <- c(xmin, xmax)

ymin <- log10(min(OtherGenes$RPKM.C, tHSR$RPKM.C, RPSgenes$RPKM.C, RPLgenes$RPKM.C))
ymax <- log10(max(OtherGenes$RPKM.C, tHSR$RPKM.C, RPSgenes$RPKM.C, RPLgenes$RPKM.C))
ylim <- c(ymin, ymax)

plot(x = log10(OtherGenes$Ave.WT.FP), y = log10(OtherGenes$RPKM.C), 
	xlab = "log10(Ave.WT.FP)", ylab = "log10(RPKM.C)", pch = 20,
	col = gray(0, alpha = 0.2), xlim = xlim, ylim = ylim, 
	main = paste("OtherGenes n=", nrow(OtherGenes), sep = ""))

plot(x = log10(tHSR$Ave.WT.FP), y = log10(tHSR$RPKM.C), 
	xlab = "log10(Ave.WT.FP)", ylab = "log10(RPKM.C)", pch = 20,
	col = rgb(0, 0, 1, alpha = 0.5), xlim = xlim, ylim = ylim, 
	main = paste("tHSR n=", nrow(tHSR), sep = ""))

plot(x = log10(RPSgenes$Ave.WT.FP), y = log10(RPSgenes$RPKM.C), 
	xlab = "log10(Ave.WT.FP)", ylab = "log10(RPKM.C)", pch = 20,
	col = rgb(1, 0, 0, alpha = 0.5), xlim = xlim, ylim = ylim, 
	main = paste("RPSgenes n=", nrow(RPSgenes), sep = ""))

plot(x = log10(RPLgenes$Ave.WT.FP), y = log10(RPLgenes$RPKM.C), 
	xlab = "log10(Ave.WT.FP)", ylab = "log10(RPKM.C)", pch = 20,
	col = rgb(0, 0.8, 0.2, alpha = 0.5), xlim = xlim, ylim = ylim, 
	main = paste("RPLgenes n=", nrow(RPLgenes), sep = ""))

plot(x = log10(All$Ave.WT.FP), y = log10(All$RPKM.C), 
	xlab = "log10(Ave.WT.FP)", ylab = "log10(RPKM.C)", pch = 20,
	col = gray(0, alpha = 0.2), xlim = xlim, ylim = ylim, 
	main = paste("All n=", nrow(All), sep = ""))

plot(x = log10(OtherGenes$Ave.WT.FP), y = log10(OtherGenes$RPKM.C), 
	xlab = "log10(Ave.WT.FP)", ylab = "log10(RPKM.C)", pch = 20,
	col = gray(0, alpha = 0.2), xlim = xlim, ylim = ylim, 
	main = paste("All n=", nrow(All), sep = ""))
points(x = log10(tHSR$Ave.WT.FP), y = log10(tHSR$RPKM.C), pch = 20,
	col = rgb(0, 0, 1, alpha = 0.5))
points(x = log10(RPSgenes$Ave.WT.FP), y = log10(RPSgenes$RPKM.C), pch = 20,
	col = rgb(1, 0, 0, alpha = 0.5))
points(x = log10(RPLgenes$Ave.WT.FP), y = log10(RPLgenes$RPKM.C), pch = 20,
	col = rgb(0, 0.8, 0.2, alpha = 0.5))

## 
xmin <- log10(min(OtherGenes$Ave.WT.mRNA, tHSR$Ave.WT.mRNA, RPSgenes$Ave.WT.mRNA, RPLgenes$Ave.WT.mRNA))
xmax <- log10(max(OtherGenes$Ave.WT.mRNA, tHSR$Ave.WT.mRNA, RPSgenes$Ave.WT.mRNA, RPLgenes$Ave.WT.mRNA))
xlim <- c(xmin, xmax)

ymin <- log10(min(OtherGenes$Ave.WT.FP, tHSR$Ave.WT.FP, RPSgenes$Ave.WT.FP, RPLgenes$Ave.WT.FP))
ymax <- log10(max(OtherGenes$Ave.WT.FP, tHSR$Ave.WT.FP, RPSgenes$Ave.WT.FP, RPLgenes$Ave.WT.FP))
ylim <- c(ymin, ymax)

plot(x = log10(OtherGenes$Ave.WT.mRNA), y = log10(OtherGenes$Ave.WT.FP), 
	xlab = "log10(Ave.WT.mRNA)", ylab = "log10(Ave.WT.FP)", pch = 20,
	col = gray(0, alpha = 0.2), xlim = xlim, ylim = ylim, 
	main = paste("OtherGenes n=", nrow(OtherGenes), sep = ""))

plot(x = log10(tHSR$Ave.WT.mRNA), y = log10(tHSR$Ave.WT.FP), 
	xlab = "log10(Ave.WT.mRNA)", ylab = "log10(Ave.WT.FP)", pch = 20,
	col = rgb(0, 0, 1, alpha = 0.5), xlim = xlim, ylim = ylim, 
	main = paste("tHSR n=", nrow(tHSR), sep = ""))

plot(x = log10(RPSgenes$Ave.WT.mRNA), y = log10(RPSgenes$Ave.WT.FP), 
	xlab = "log10(Ave.WT.mRNA)", ylab = "log10(Ave.WT.FP)", pch = 20,
	col = rgb(1, 0, 0, alpha = 0.5), xlim = xlim, ylim = ylim, 
	main = paste("RPSgenes n=", nrow(RPSgenes), sep = ""))

plot(x = log10(RPLgenes$Ave.WT.mRNA), y = log10(RPLgenes$Ave.WT.FP), 
	xlab = "log10(Ave.WT.mRNA)", ylab = "log10(Ave.WT.FP)", pch = 20,
	col = rgb(0, 0.8, 0.2, alpha = 0.5), xlim = xlim, ylim = ylim, 
	main = paste("RPLgenes n=", nrow(RPLgenes), sep = ""))

plot(x = log10(All$Ave.WT.mRNA), log10(All$Ave.WT.FP), 
	xlab = "log10(Ave.WT.mRNA)", ylab = "log10(Ave.WT.FP)", pch = 20,
	col = gray(0, alpha = 0.2), xlim = xlim, ylim = ylim, 
	main = paste("All n=", nrow(All), sep = ""))

plot(x = log10(OtherGenes$Ave.WT.mRNA), y = log10(OtherGenes$Ave.WT.FP), pch = 20,
	col = gray(0, alpha = 0.2), xlim = xlim, ylim = ylim, 
	main = paste("All n=", nrow(All), sep = ""))
points(x = log10(tHSR$Ave.WT.mRNA), y = log10(tHSR$Ave.WT.FP), pch = 20,
	col = rgb(0, 0, 1, alpha = 0.5))
points(x = log10(RPSgenes$Ave.WT.mRNA), y = log10(RPSgenes$Ave.WT.FP), pch = 20,
	col = rgb(1, 0, 0, alpha = 0.5))
points(x = log10(RPLgenes$Ave.WT.mRNA), y = log10(RPLgenes$Ave.WT.FP), pch = 20,
	col = rgb(0, 0.8, 0.2, alpha = 0.5))

## 
xmin <- log10(min(OtherGenes$Ave.WT.mRNA, tHSR$Ave.WT.mRNA, RPSgenes$Ave.WT.mRNA, RPLgenes$Ave.WT.mRNA))
xmax <- log10(max(OtherGenes$Ave.WT.mRNA, tHSR$Ave.WT.mRNA, RPSgenes$Ave.WT.mRNA, RPLgenes$Ave.WT.mRNA))
xlim <- c(xmin, xmax)

ymin <- log10(min(OtherGenes$TE.DK_WT, tHSR$TE.DK_WT, RPSgenes$TE.DK_WT, RPLgenes$TE.DK_WT))
ymax <- log10(max(OtherGenes$TE.DK_WT, tHSR$TE.DK_WT, RPSgenes$TE.DK_WT, RPLgenes$TE.DK_WT))
ylim <- c(ymin, ymax)

plot(x = log10(OtherGenes$Ave.WT.mRNA), y = log10(OtherGenes$TE.DK_WT), 
	xlab = "log10(Ave.WT.mRNA)", ylab = "log10(TE.DK_WT)", pch = 20,
	col = gray(0, alpha = 0.2), xlim = xlim, ylim = ylim, 
	main = paste("OtherGenes n=", nrow(OtherGenes), sep = ""))

plot(x = log10(tHSR$Ave.WT.mRNA), y = log10(tHSR$TE.DK_WT), 
	xlab = "log10(Ave.WT.mRNA)", ylab = "log10(TE.DK_WT)", pch = 20,
	col = rgb(0, 0, 1, alpha = 0.5), xlim = xlim, ylim = ylim, 
	main = paste("tHSR n=", nrow(tHSR), sep = ""))

plot(x = log10(RPSgenes$Ave.WT.mRNA), y = log10(RPSgenes$TE.DK_WT), 
	xlab = "log10(Ave.WT.mRNA)", ylab = "log10(TE.DK_WT)", pch = 20,
	col = rgb(1, 0, 0, alpha = 0.5), xlim = xlim, ylim = ylim, 
	main = paste("RPSgenes n=", nrow(RPSgenes), sep = ""))

plot(x = log10(RPLgenes$Ave.WT.mRNA), y = log10(RPLgenes$TE.DK_WT), 
	xlab = "log10(Ave.WT.mRNA)", ylab = "log10(TE.DK_WT)", pch = 20,
	col = rgb(0, 0.8, 0.2, alpha = 0.5), xlim = xlim, ylim = ylim, 
	main = paste("RPLgenes n=", nrow(RPLgenes), sep = ""))

plot(x = log10(All$Ave.WT.mRNA), y = log10(All$TE.DK_WT), 
	xlab = "log10(Ave.WT.mRNA)", ylab = "log10(TE.DK_WT)", pch = 20,
	col = gray(0, alpha = 0.2), xlim = xlim, ylim = ylim, 
	main = paste("All n=", nrow(All), sep = ""))

plot(x = log10(OtherGenes$Ave.WT.mRNA), y = log10(OtherGenes$TE.DK_WT), pch = 20,
	col = gray(0, alpha = 0.2), xlim = xlim, ylim = ylim, 
	main = paste("All n=", nrow(All), sep = ""))
points(x = log10(tHSR$Ave.WT.mRNA), y = log10(tHSR$TE.DK_WT), pch = 20,
	col = rgb(0, 0, 1, alpha = 0.5))
points(x = log10(RPSgenes$Ave.WT.mRNA), y = log10(RPSgenes$TE.DK_WT), pch = 20,
	col = rgb(1, 0, 0, alpha = 0.5))
points(x = log10(RPLgenes$Ave.WT.mRNA), y = log10(RPLgenes$TE.DK_WT), pch = 20,
	col = rgb(0, 0.8, 0.2, alpha = 0.5))

## 
xmin <- log10(min(OtherGenes$Ave.WT.mRNA, tHSR$Ave.WT.mRNA, RPSgenes$Ave.WT.mRNA, RPLgenes$Ave.WT.mRNA))
xmax <- log10(max(OtherGenes$Ave.WT.mRNA, tHSR$Ave.WT.mRNA, RPSgenes$Ave.WT.mRNA, RPLgenes$Ave.WT.mRNA))
xlim <- c(xmin, xmax)

ymin <- log10(min(OtherGenes$TE.C.RPKM, tHSR$TE.C.RPKM, RPSgenes$TE.C.RPKM, RPLgenes$TE.C.RPKM))
ymax <- log10(max(OtherGenes$TE.C.RPKM, tHSR$TE.C.RPKM, RPSgenes$TE.C.RPKM, RPLgenes$TE.C.RPKM))
ylim <- c(ymin, ymax)

plot(x = log10(OtherGenes$Ave.WT.mRNA), y = log10(OtherGenes$TE.C.RPKM), 
	xlab = "log10(Ave.WT.mRNA)", ylab = "log10(TE.C.RPKM)", pch = 20,
	col = gray(0, alpha = 0.2), xlim = xlim, ylim = ylim, 
	main = paste("OtherGenes n=", nrow(OtherGenes), sep = ""))

plot(x = log10(tHSR$Ave.WT.mRNA), y = log10(tHSR$TE.C.RPKM), 
	xlab = "log10(Ave.WT.mRNA)", ylab = "log10(TE.C.RPKM)", pch = 20,
	col = rgb(0, 0, 1, alpha = 0.5), xlim = xlim, ylim = ylim, 
	main = paste("tHSR n=", nrow(tHSR), sep = ""))

plot(x = log10(RPSgenes$Ave.WT.mRNA), y = log10(RPSgenes$TE.C.RPKM), 
	xlab = "log10(Ave.WT.mRNA)", ylab = "log10(TE.C.RPKM)", pch = 20,
	col = rgb(1, 0, 0, alpha = 0.5), xlim = xlim, ylim = ylim, 
	main = paste("RPSgenes n=", nrow(RPSgenes), sep = ""))

plot(x = log10(RPLgenes$Ave.WT.mRNA), y = log10(RPLgenes$TE.C.RPKM), 
	xlab = "log10(Ave.WT.mRNA)", ylab = "log10(TE.C.RPKM)", pch = 20,
	col = rgb(0, 0.8, 0.2, alpha = 0.5), xlim = xlim, ylim = ylim, 
	main = paste("RPLgenes n=", nrow(RPLgenes), sep = ""))

plot(x = log10(All$Ave.WT.mRNA), y = log10(All$TE.C.RPKM), 
	xlab = "log10(Ave.WT.mRNA)", ylab = "log10(TE.C.RPKM)", pch = 20,
	col = gray(0, alpha = 0.2), xlim = xlim, ylim = ylim, 
	main = paste("All n=", nrow(All), sep = ""))

plot(x = log10(OtherGenes$Ave.WT.mRNA), y = log10(OtherGenes$TE.C.RPKM), 
	xlab = "log10(Ave.WT.mRNA)", ylab = "log10(TE.C.RPKM)", pch = 20,
	col = gray(0, alpha = 0.2), xlim = xlim, ylim = ylim, 
	main = paste("All n=", nrow(All), sep = ""))
points(x = log10(tHSR$Ave.WT.mRNA), y = log10(tHSR$TE.C.RPKM), pch = 20,
	col = rgb(0, 0, 1, alpha = 0.5))
points(x = log10(RPSgenes$Ave.WT.mRNA), y = log10(RPSgenes$TE.C.RPKM), pch = 20,
	col = rgb(1, 0, 0, alpha = 0.5))
points(x = log10(RPLgenes$Ave.WT.mRNA), y = log10(RPLgenes$TE.C.RPKM), pch = 20,
	col = rgb(0, 0.8, 0.2, alpha = 0.5))

setwd(HOMEDIR)
setwd("data/data_plots")
write.xlsx(list(OtherGenes = OtherGenes, tHSR = tHSR, RPSgenes = RPSgenes, RPLgenes = RPLgenes, 
	All = All), file = "20250614_tHSRandRP_genes.xlsx")


## tHSR RPS RPL
## 
xmin <- log10(min(OtherGenes$Ave.WT.FP, tHSR$Ave.WT.FP, RPSgenes$Ave.WT.FP, RPLgenes$Ave.WT.FP))
xmax <- log10(max(OtherGenes$Ave.WT.FP, tHSR$Ave.WT.FP, RPSgenes$Ave.WT.FP, RPLgenes$Ave.WT.FP))
xlim <- c(xmin, xmax)

ymin <- log10(min(OtherGenes$RPKM.C, tHSR$RPKM.C, RPSgenes$RPKM.C, RPLgenes$RPKM.C))
ymax <- log10(max(OtherGenes$RPKM.C, tHSR$RPKM.C, RPSgenes$RPKM.C, RPLgenes$RPKM.C))
ylim <- c(ymin, ymax)

temp <- data.frame(Ave.WT.FP = All$Ave.WT.FP, RPKM.C = All$RPKM.C)
temp.lm <- lm(log10(RPKM.C) ~ log10(Ave.WT.FP), data = temp)
summary(temp.lm)

# > summary(temp.lm)
# 
# Call:
# lm(formula = log10(RPKM.C) ~ log10(Ave.WT.FP), data = temp)
# 
# Residuals:
#     Min      1Q  Median      3Q     Max 
# -2.3740 -0.2656  0.0066  0.2630  3.3534 
# 
# Coefficients:
#                  Estimate Std. Error t value Pr(>|t|)    
# (Intercept)      0.389522   0.016183   24.07   <2e-16 ***
# log10(Ave.WT.FP) 0.786562   0.009159   85.88   <2e-16 ***
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.4426 on 4891 degrees of freedom
# Multiple R-squared:  0.6013,	Adjusted R-squared:  0.6012 
# F-statistic:  7375 on 1 and 4891 DF,  p-value: < 2.2e-16

# intercept: 0.389522
# slope: 0.786562
# R^2: 0.6012

plot(x = log10(OtherGenes$Ave.WT.FP), y = log10(OtherGenes$RPKM.C), 
	xlab = "log 10 (footprint density 30 min in heat), Stanciu et al", 
	ylab = "log 10 (footprint density 3 hr  in heat), This study", pch = 20,
	col = gray(0, alpha = 0.2), xlim = xlim, ylim = ylim, 
	main = paste("All n=", nrow(All), sep = ""))
points(x = log10(tHSR$Ave.WT.FP), y = log10(tHSR$RPKM.C), pch = 20,
	col = rgb(0, 0, 1, alpha = 0.5))
points(x = log10(RPSgenes$Ave.WT.FP), y = log10(RPSgenes$RPKM.C), pch = 20,
	col = rgb(1, 0, 0, alpha = 0.5))
points(x = log10(RPLgenes$Ave.WT.FP), y = log10(RPLgenes$RPKM.C), pch = 20,
	col = rgb(0, 0.8, 0.2, alpha = 0.5))
abline(temp.lm)


## 
xmin <- log10(min(OtherGenes$Ave.WT.mRNA, tHSR$Ave.WT.mRNA, RPSgenes$Ave.WT.mRNA, RPLgenes$Ave.WT.mRNA))
xmax <- log10(max(OtherGenes$Ave.WT.mRNA, tHSR$Ave.WT.mRNA, RPSgenes$Ave.WT.mRNA, RPLgenes$Ave.WT.mRNA))
xlim <- c(xmin, xmax)

ymin <- log10(min(OtherGenes$Ave.WT.FP, tHSR$Ave.WT.FP, RPSgenes$Ave.WT.FP, RPLgenes$Ave.WT.FP))
ymax <- log10(max(OtherGenes$Ave.WT.FP, tHSR$Ave.WT.FP, RPSgenes$Ave.WT.FP, RPLgenes$Ave.WT.FP))
ylim <- c(ymin, ymax)


temp <- data.frame(Ave.WT.mRNA = All$Ave.WT.mRNA, Ave.WT.FP = All$Ave.WT.FP)
temp.lm <- lm(log10(Ave.WT.FP) ~ log10(Ave.WT.mRNA), data = temp)
summary(temp.lm)

# > summary(temp.lm)
# 
# Call:
# lm(formula = log10(Ave.WT.FP) ~ log10(Ave.WT.mRNA), data = temp)
# 
# Residuals:
#      Min       1Q   Median       3Q      Max 
# -2.39777 -0.17163  0.01494  0.19134  2.43951 
# 
# Coefficients:
#                     Estimate Std. Error t value Pr(>|t|)    
# (Intercept)        -0.861237   0.018202  -47.31   <2e-16 ***
# log10(Ave.WT.mRNA)  1.333519   0.009469  140.82   <2e-16 ***
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.3073 on 4891 degrees of freedom
# Multiple R-squared:  0.8022,	Adjusted R-squared:  0.8021 
# F-statistic: 1.983e+04 on 1 and 4891 DF,  p-value: < 2.2e-16

# intercept: -0.861237
# slope: 1.333519
# R^2: 0.8021

plot(x = log10(OtherGenes$Ave.WT.mRNA), y = log10(OtherGenes$Ave.WT.FP), pch = 20,
	col = gray(0, alpha = 0.2), xlim = xlim, ylim = ylim, 
	xlab = "log 10 (mRNA abundance in heat), Stanciu et al",
	ylab = "log 10 (footprint density  in heat), Stanciu et al",
	main = paste("All n=", nrow(All), sep = ""))
points(x = log10(tHSR$Ave.WT.mRNA), y = log10(tHSR$Ave.WT.FP), pch = 20,
	col = rgb(0, 0, 1, alpha = 0.5))
points(x = log10(RPSgenes$Ave.WT.mRNA), y = log10(RPSgenes$Ave.WT.FP), pch = 20,
	col = rgb(1, 0, 0, alpha = 0.5))
points(x = log10(RPLgenes$Ave.WT.mRNA), y = log10(RPLgenes$Ave.WT.FP), pch = 20,
	col = rgb(0, 0.8, 0.2, alpha = 0.5))
abline(temp.lm)


## 
xmin <- log10(min(OtherGenes$Ave.WT.mRNA, tHSR$Ave.WT.mRNA, RPSgenes$Ave.WT.mRNA, RPLgenes$Ave.WT.mRNA))
xmax <- log10(max(OtherGenes$Ave.WT.mRNA, tHSR$Ave.WT.mRNA, RPSgenes$Ave.WT.mRNA, RPLgenes$Ave.WT.mRNA))
xlim <- c(xmin, xmax)

ymin <- log10(min(OtherGenes$TE.DK_WT, tHSR$TE.DK_WT, RPSgenes$TE.DK_WT, RPLgenes$TE.DK_WT))
ymax <- log10(max(OtherGenes$TE.DK_WT, tHSR$TE.DK_WT, RPSgenes$TE.DK_WT, RPLgenes$TE.DK_WT))
ylim <- c(ymin, ymax)

temp <- data.frame(Ave.WT.mRNA = All$Ave.WT.mRNA, TE.DK_WT = All$TE.DK_WT)
temp.lm <- lm(log10(TE.DK_WT) ~ log10(Ave.WT.mRNA), data = temp)
summary(temp.lm)

# > summary(temp.lm)
# 
# Call:
# lm(formula = log10(TE.DK_WT) ~ log10(Ave.WT.mRNA), data = temp)
# 
# Residuals:
#      Min       1Q   Median       3Q      Max 
# -2.39777 -0.17163  0.01494  0.19134  2.43951 
# 
# Coefficients:
#                     Estimate Std. Error t value Pr(>|t|)    
# (Intercept)        -0.861237   0.018202  -47.31   <2e-16 ***
# log10(Ave.WT.mRNA)  0.333519   0.009469   35.22   <2e-16 ***
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.3073 on 4891 degrees of freedom
# Multiple R-squared:  0.2023,	Adjusted R-squared:  0.2022 
# F-statistic:  1240 on 1 and 4891 DF,  p-value: < 2.2e-16

# intercept: -0.861237
# slope: 0.333519
# R^2: 0.2022

plot(x = log10(OtherGenes$Ave.WT.mRNA), y = log10(OtherGenes$TE.DK_WT), pch = 20,
	col = gray(0, alpha = 0.2), xlim = xlim, ylim = ylim, 
	xlab = "log 10 (mRNA abundance in heat), Stanciu et al", 
	ylab = "log 10 (footprint density  in heat), Stanciu et al", 
	main = paste("All n=", nrow(All), sep = ""))
points(x = log10(tHSR$Ave.WT.mRNA), y = log10(tHSR$TE.DK_WT), pch = 20,
	col = rgb(0, 0, 1, alpha = 0.5))
points(x = log10(RPSgenes$Ave.WT.mRNA), y = log10(RPSgenes$TE.DK_WT), pch = 20,
	col = rgb(1, 0, 0, alpha = 0.5))
points(x = log10(RPLgenes$Ave.WT.mRNA), y = log10(RPLgenes$TE.DK_WT), pch = 20,
	col = rgb(0, 0.8, 0.2, alpha = 0.5))
abline(temp.lm)


## 
xmin <- log10(min(OtherGenes$Ave.WT.mRNA, tHSR$Ave.WT.mRNA, RPSgenes$Ave.WT.mRNA, RPLgenes$Ave.WT.mRNA))
xmax <- log10(max(OtherGenes$Ave.WT.mRNA, tHSR$Ave.WT.mRNA, RPSgenes$Ave.WT.mRNA, RPLgenes$Ave.WT.mRNA))
xlim <- c(xmin, xmax)

ymin <- log10(min(OtherGenes$TE.C.RPKM, tHSR$TE.C.RPKM, RPSgenes$TE.C.RPKM, RPLgenes$TE.C.RPKM))
ymax <- log10(max(OtherGenes$TE.C.RPKM, tHSR$TE.C.RPKM, RPSgenes$TE.C.RPKM, RPLgenes$TE.C.RPKM))
ylim <- c(ymin, ymax)

temp <- data.frame(Ave.WT.mRNA = All$Ave.WT.mRNA, TE.C.RPKM = All$TE.C.RPKM)
temp.lm <- lm(log10(TE.C.RPKM) ~ log10(Ave.WT.mRNA), data = temp)
summary(temp.lm)

# > summary(temp.lm)
# 
# Call:
# lm(formula = log10(TE.C.RPKM) ~ log10(Ave.WT.mRNA), data = temp)
# 
# Residuals:
#     Min      1Q  Median      3Q     Max 
# -2.3196 -0.3212  0.0072  0.3271  3.4103 
# 
# Coefficients:
#                    Estimate Std. Error t value Pr(>|t|)    
# (Intercept)        -0.34566    0.02902  -11.91  < 2e-16 ***
# log10(Ave.WT.mRNA)  0.07986    0.01510    5.29 1.28e-07 ***
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.49 on 4891 degrees of freedom
# Multiple R-squared:  0.005688,	Adjusted R-squared:  0.005485 
# F-statistic: 27.98 on 1 and 4891 DF,  p-value: 1.279e-07

# intercept: -0.34566
# slope: 0.07986
# R^2: 0.005485

plot(x = log10(OtherGenes$Ave.WT.mRNA), y = log10(OtherGenes$TE.C.RPKM), pch = 20,
	xlab = "log 10 (mRNA abundance in heat), Stanciu et al", 
	ylab = "log 10 (footprint density  in heat), Stanciu et al", 
	col = gray(0, alpha = 0.2), xlim = xlim, ylim = ylim, 
	main = paste("All n=", nrow(All), sep = ""))
points(x = log10(tHSR$Ave.WT.mRNA), y = log10(tHSR$TE.C.RPKM), pch = 20,
	col = rgb(0, 0, 1, alpha = 0.5))
points(x = log10(RPSgenes$Ave.WT.mRNA), y = log10(RPSgenes$TE.C.RPKM), pch = 20,
	col = rgb(1, 0, 0, alpha = 0.5))
points(x = log10(RPLgenes$Ave.WT.mRNA), y = log10(RPLgenes$TE.C.RPKM), pch = 20,
	col = rgb(0, 0.8, 0.2, alpha = 0.5))
abline(temp.lm)



