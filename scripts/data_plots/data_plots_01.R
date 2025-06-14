## original note: 20250614_AsanoTE_1.txt


library(openxlsx)

HOMEDIR <- ""
setwd(HOMEDIR)
setwd("data/data_plots")
RPKMtable <- read.xlsx(xlsxFile = "Yeast_CHX_KSU_RPKM_table_2_hk2.xlsx")

setwd(HOMEDIR)
setwd("data/data_plots")
NewDF <- read.xlsx(xlsxFile = "20240522_RPKMtable_merged_with_pone_match4496_with444.xlsx")

setwd(HOMEDIR)
setwd("data/data_plots")
Top64 <- readLines(con = "tif34 down top 64.txt")

NewDF$Top64 <- FALSE
for(i in 1:64){
	NewDF$Top64[NewDF$Gene == Top64[i]] <- TRUE
}

NewDF <- cbind(RPKMtable[, 1:12], NewDF[, 2:23])

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

setwd(HOMEDIR)
setwd("data/data_plots")
RPLgenes <- read.xlsx(xlsxFile = "20250614_RPgenes.xlsx", sheet = "RPL", colNames = FALSE)[, c(1, 15)]
RPSgenes <- read.xlsx(xlsxFile = "20250614_RPgenes.xlsx", sheet = "RPS", colNames = FALSE)[, c(1, 15)]

head(RPLgenes, n = 3)
head(RPSgenes, n = 3)

nrow(RPLgenes)
nrow(RPSgenes)

for(i in 1:nrow(RPLgenes)){
	NewDF$RPL[NewDF$Gene == RPLgenes[i, 1]] <- TRUE
}

for(i in 1:nrow(RPSgenes)){
	NewDF$RPS[NewDF$Gene == RPSgenes[i, 1]] <- TRUE
}

NewDF$TIF34_TE <- NewDF$RPKM.C / NewDF$RA.t40
NewDF$txn_fold <- NewDF$RA.t40 / NewDF$RA.t0

All <- subset(NewDF, is.na(RA.t40) == FALSE & is.na(RA.t0) == FALSE & is.na(TIF34_TE) == FALSE)
Only246 <- subset(All, S4_246 == 1)
No246 <- subset(All, is.na(S4_246) == TRUE)
Top246 <- subset(All, tif34down_top50 == TRUE)
Bottom246 <- subset(All, tif34down_bottom50 == TRUE)
Middle246 <- subset(All, S4_246 == 1 & tif34down_top50 == FALSE & tif34down_bottom50 == FALSE)
NoRPno64genes <- subset(All, RPS == FALSE & RPL == FALSE & Top64 == FALSE)

RPSgenes <- subset(All, RPS == TRUE)
RPLgenes <- subset(All, RPL == TRUE)
Top64genes <- subset(All, Top64 == TRUE)

setwd(HOMEDIR)
setwd("data/data_plots")

write.xlsx(list(NoRPno64genes = NoRPno64genes, RPSgenes = RPSgenes, RPLgenes = RPLgenes, 
	Top64genes = Top64genes, All = All), file = "20250614_RP_Top64_genes.xlsx")
save(NewDF, file = "20250614_NewDF.RData")

#############
# Figure 2C, D, E

# 2C
ylim <- xlim <- c(0, max(log10(All$RPKM.C), log10(All$RPKM.TVmean)))

plot(x = NA, y = NA, 
	xlim = xlim, ylim = ylim, 
	xlab = "log10(RPKM.C)", ylab = "log10(RPKM.TVmean)", 
	main = paste("All n=", nrow(All), sep = ""))
# abline(a = 0, b = 1, col = "gray")
points(x = log10(NoRPno64genes$RPKM.C), y = log10(NoRPno64genes$RPKM.TVmean), pch = 20, 
	col = gray(0, alpha = 0.2))
points(x = log10(RPSgenes$RPKM.C), y = log10(RPSgenes$RPKM.TVmean), pch = 20, 
	col = rgb(1, 0, 0, alpha = 0.5))
points(x = log10(RPLgenes$RPKM.C), y = log10(RPLgenes$RPKM.TVmean), pch = 20, 
	col = rgb(0, 0.8, 0.2, alpha = 0.5))
points(x = log10(Top64genes$RPKM.C), y = log10(Top64genes$RPKM.TVmean), pch = 20, 
	col = rgb(0, 0, 1, alpha = 0.5))

y = c(log10(NoRPno64genes$RPKM.TVmean), log10(RPSgenes$RPKM.TVmean), log10(RPLgenes$RPKM.TVmean), log10(Top64genes$RPKM.TVmean))
x = c(log10(NoRPno64genes$RPKM.C), log10(RPSgenes$RPKM.C), log10(RPLgenes$RPKM.C), log10(Top64genes$RPKM.C))

linReg <- lm(y ~ x)

# > summary(linReg)
# 
# Call:
# lm(formula = y ~ x)
# 
# Residuals:
#      Min       1Q   Median       3Q      Max 
# -0.83134 -0.11742 -0.01886  0.10527  0.84580 
# 
# Coefficients:
#             Estimate Std. Error t value Pr(>|t|)    
# (Intercept) 0.258543   0.006825   37.88   <2e-16 ***
# x           0.907181   0.003772  240.47   <2e-16 ***
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.175 on 4459 degrees of freedom
# Multiple R-squared:  0.9284,	Adjusted R-squared:  0.9284 
# F-statistic: 5.783e+04 on 1 and 4459 DF,  p-value: < 2.2e-16

abline(linReg, col = "orange", lwd = 2)



# 2D
ylim <- xlim <- c(0, max(log10(All$RPKM.C), log10(All$RPKM.TT)))

plot(x = NA, y = NA, 
	xlim = xlim, ylim = ylim, 
	xlab = "log10(RPKM.C)", ylab = "log10(RPKM.TT)", 
	main = paste("All n=", nrow(All), sep = ""))
# abline(a = 0, b = 1, col = "gray")
points(x = log10(NoRPno64genes$RPKM.C), y = log10(NoRPno64genes$RPKM.TT), pch = 20, 
	col = gray(0, alpha = 0.2))
points(x = log10(RPSgenes$RPKM.C), y = log10(RPSgenes$RPKM.TT), pch = 20, 
	col = rgb(1, 0, 0, alpha = 0.5))
points(x = log10(RPLgenes$RPKM.C), y = log10(RPLgenes$RPKM.TT), pch = 20, 
	col = rgb(0, 0.8, 0.2, alpha = 0.5))
points(x = log10(Top64genes$RPKM.C), y = log10(Top64genes$RPKM.TT), pch = 20, 
	col = rgb(0, 0, 1, alpha = 0.5))

y = c(log10(NoRPno64genes$RPKM.TT), log10(RPSgenes$RPKM.TT), log10(RPLgenes$RPKM.TT), log10(Top64genes$RPKM.TT))
x = c(log10(NoRPno64genes$RPKM.C), log10(RPSgenes$RPKM.C), log10(RPLgenes$RPKM.C), log10(Top64genes$RPKM.C))

linReg <- lm(y ~ x)

# > summary(linReg)
# 
# Call:
# lm(formula = y ~ x)
# 
# Residuals:
#      Min       1Q   Median       3Q      Max 
# -0.63258 -0.09552 -0.01708  0.08319  1.03288 
# 
# Coefficients:
#             Estimate Std. Error t value Pr(>|t|)    
# (Intercept) 0.130532   0.005390   24.22   <2e-16 ***
# x           0.953299   0.002979  319.99   <2e-16 ***
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.1382 on 4459 degrees of freedom
# Multiple R-squared:  0.9583,	Adjusted R-squared:  0.9583 
# F-statistic: 1.024e+05 on 1 and 4459 DF,  p-value: < 2.2e-16

abline(linReg, col = "orange", lwd = 2)


# TV1 vs TV2
ylim <- xlim <- c(0, max(log10(All$RPKM.TV1), log10(All$RPKM.TV2)))

plot(x = NA, y = NA, 
	xlim = xlim, ylim = ylim, 
	xlab = "log10(RPKM.TV1)", ylab = "log10(RPKM.TV2)", 
	main = paste("All n=", nrow(All), sep = ""))
# abline(a = 0, b = 1, col = "gray")
points(x = log10(NoRPno64genes$RPKM.TV1), y = log10(NoRPno64genes$RPKM.TV2), pch = 20, 
	col = gray(0, alpha = 0.2))
points(x = log10(RPSgenes$RPKM.TV1), y = log10(RPSgenes$RPKM.TV2), pch = 20, 
	col = rgb(1, 0, 0, alpha = 0.5))
points(x = log10(RPLgenes$RPKM.TV1), y = log10(RPLgenes$RPKM.TV2), pch = 20, 
	col = rgb(0, 0.8, 0.2, alpha = 0.5))
points(x = log10(Top64genes$RPKM.TV1), y = log10(Top64genes$RPKM.TV2), pch = 20, 
	col = rgb(0, 0, 1, alpha = 0.5))

y = c(log10(NoRPno64genes$RPKM.TV2), log10(RPSgenes$RPKM.TV2), log10(RPLgenes$RPKM.TV2), log10(Top64genes$RPKM.TV2))
x = c(log10(NoRPno64genes$RPKM.TV1), log10(RPSgenes$RPKM.TV1), log10(RPLgenes$RPKM.TV1), log10(Top64genes$RPKM.TV1))

linReg <- lm(y ~ x)

# > summary(linReg)
# 
# Call:
# lm(formula = y ~ x)
# 
# Residuals:
#      Min       1Q   Median       3Q      Max 
# -0.62279 -0.05354  0.00388  0.05869  1.13115 
# 
# Coefficients:
#              Estimate Std. Error t value Pr(>|t|)    
# (Intercept) -0.040479   0.004205  -9.626   <2e-16 ***
# x            1.008652   0.002215 455.420   <2e-16 ***
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.09619 on 4459 degrees of freedom
# Multiple R-squared:  0.979,	Adjusted R-squared:  0.9789 
# F-statistic: 2.074e+05 on 1 and 4459 DF,  p-value: < 2.2e-16

abline(linReg, col = "orange", lwd = 2)


# 2E
xlim <- c(0, max(All$tif34.down.fold))
ylim <- c(0, max(All$tif35.up.fold))

plot(x = NA, y = NA, 
	xlim = xlim, ylim = ylim, 
	xlab = "tif34.down.fold)", ylab = "tif35.up.fold)", 
	main = paste("All n=", nrow(All), sep = ""))
# abline(a = 0, b = 1, col = "gray")
points(x = NoRPno64genes$tif34.down.fold, y = NoRPno64genes$tif35.up.fold, pch = 20, 
	col = gray(0, alpha = 0.2))
points(x = RPSgenes$tif34.down.fold, y = RPSgenes$tif35.up.fold, pch = 20, 
	col = rgb(1, 0, 0, alpha = 0.5))
points(x = RPLgenes$tif34.down.fold, y = RPLgenes$tif35.up.fold, pch = 20, 
	col = rgb(0, 0.8, 0.2, alpha = 0.5))
points(x = Top64genes$tif34.down.fold, y = Top64genes$tif35.up.fold, pch = 20, 
	col = rgb(0, 0, 1, alpha = 0.5))

y = c(NoRPno64genes$tif35.up.fold, RPSgenes$tif35.up.fold, RPLgenes$tif35.up.fold, Top64genes$tif35.up.fold)
x = c(NoRPno64genes$tif34.down.fold, RPSgenes$tif34.down.fold, RPLgenes$tif34.down.fold, Top64genes$tif34.down.fold)

linReg <- lm(y ~ x)

# > summary(linReg)
# 
# Call:
# lm(formula = y ~ x)
# 
# Residuals:
#     Min      1Q  Median      3Q     Max 
# -0.7903 -0.0999 -0.0110  0.0791  3.5806 
# 
# Coefficients:
#             Estimate Std. Error t value Pr(>|t|)    
# (Intercept) 0.592431   0.006489   91.30   <2e-16 ***
# x           0.373254   0.006980   53.47   <2e-16 ***
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.1657 on 4459 degrees of freedom
# Multiple R-squared:  0.3907,	Adjusted R-squared:  0.3906 
# F-statistic:  2859 on 1 and 4459 DF,  p-value: < 2.2e-16

abline(linReg, col = "orange", lwd = 2)


###################
## linear regression line 

# 2C
ylim <- xlim <- c(0, max(log10(All$RPKM.C), log10(All$RPKM.TVmean)))

# 20250512
y = c(log10(NoRPno64genes$RPKM.TVmean), log10(RPSgenes$RPKM.TVmean), log10(RPLgenes$RPKM.TVmean), log10(Top64genes$RPKM.TVmean))
x = c(log10(NoRPno64genes$RPKM.C), log10(RPSgenes$RPKM.C), log10(RPLgenes$RPKM.C), log10(Top64genes$RPKM.C))

linReg <- lm(y ~ x)

plot(x = NA, y = NA, 
	xlim = xlim, ylim = ylim, 
	xlab = "log10(RPKM.C)", ylab = "log10(RPKM.TVmean)", 
	main = paste("All n=", nrow(All), sep = ""))
# abline(a = 0, b = 1, col = "gray")
abline(linReg, col = "gray", lwd = 2)
points(x = log10(NoRPno64genes$RPKM.C), y = log10(NoRPno64genes$RPKM.TVmean), pch = 20, 
	col = gray(0, alpha = 0.2))
points(x = log10(RPSgenes$RPKM.C), y = log10(RPSgenes$RPKM.TVmean), pch = 20, 
	col = rgb(1, 0, 0, alpha = 0.5))
points(x = log10(RPLgenes$RPKM.C), y = log10(RPLgenes$RPKM.TVmean), pch = 20, 
	col = rgb(0, 0.8, 0.2, alpha = 0.5))
points(x = log10(Top64genes$RPKM.C), y = log10(Top64genes$RPKM.TVmean), pch = 20, 
	col = rgb(0, 0, 1, alpha = 0.5))



# 2D
ylim <- xlim <- c(0, max(log10(All$RPKM.C), log10(All$RPKM.TT)))

# 20250512
y = c(log10(NoRPno64genes$RPKM.TT), log10(RPSgenes$RPKM.TT), log10(RPLgenes$RPKM.TT), log10(Top64genes$RPKM.TT))
x = c(log10(NoRPno64genes$RPKM.C), log10(RPSgenes$RPKM.C), log10(RPLgenes$RPKM.C), log10(Top64genes$RPKM.C))

linReg <- lm(y ~ x)

plot(x = NA, y = NA, 
	xlim = xlim, ylim = ylim, 
	xlab = "log10(RPKM.C)", ylab = "log10(RPKM.TT)", 
	main = paste("All n=", nrow(All), sep = ""))
# abline(a = 0, b = 1, col = "gray")
abline(linReg, col = "gray", lwd = 2)
points(x = log10(NoRPno64genes$RPKM.C), y = log10(NoRPno64genes$RPKM.TT), pch = 20, 
	col = gray(0, alpha = 0.2))
points(x = log10(RPSgenes$RPKM.C), y = log10(RPSgenes$RPKM.TT), pch = 20, 
	col = rgb(1, 0, 0, alpha = 0.5))
points(x = log10(RPLgenes$RPKM.C), y = log10(RPLgenes$RPKM.TT), pch = 20, 
	col = rgb(0, 0.8, 0.2, alpha = 0.5))
points(x = log10(Top64genes$RPKM.C), y = log10(Top64genes$RPKM.TT), pch = 20, 
	col = rgb(0, 0, 1, alpha = 0.5))


# TV1 vs TV2
ylim <- xlim <- c(0, max(log10(All$RPKM.TV1), log10(All$RPKM.TV2)))

# 20250512
y = c(log10(NoRPno64genes$RPKM.TV2), log10(RPSgenes$RPKM.TV2), log10(RPLgenes$RPKM.TV2), log10(Top64genes$RPKM.TV2))
x = c(log10(NoRPno64genes$RPKM.TV1), log10(RPSgenes$RPKM.TV1), log10(RPLgenes$RPKM.TV1), log10(Top64genes$RPKM.TV1))

linReg <- lm(y ~ x)

plot(x = NA, y = NA, 
	xlim = xlim, ylim = ylim, 
	xlab = "log10(RPKM.TV1)", ylab = "log10(RPKM.TV2)", 
	main = paste("All n=", nrow(All), sep = ""))
# abline(a = 0, b = 1, col = "gray")
abline(linReg, col = "gray", lwd = 2)
points(x = log10(NoRPno64genes$RPKM.TV1), y = log10(NoRPno64genes$RPKM.TV2), pch = 20, 
	col = gray(0, alpha = 0.2))
points(x = log10(RPSgenes$RPKM.TV1), y = log10(RPSgenes$RPKM.TV2), pch = 20, 
	col = rgb(1, 0, 0, alpha = 0.5))
points(x = log10(RPLgenes$RPKM.TV1), y = log10(RPLgenes$RPKM.TV2), pch = 20, 
	col = rgb(0, 0.8, 0.2, alpha = 0.5))
points(x = log10(Top64genes$RPKM.TV1), y = log10(Top64genes$RPKM.TV2), pch = 20, 
	col = rgb(0, 0, 1, alpha = 0.5))


# 2E
xlim <- c(0, max(All$tif34.down.fold))
ylim <- c(0, max(All$tif35.up.fold))

# 20250512
y = c(NoRPno64genes$tif35.up.fold, RPSgenes$tif35.up.fold, RPLgenes$tif35.up.fold, Top64genes$tif35.up.fold)
x = c(NoRPno64genes$tif34.down.fold, RPSgenes$tif34.down.fold, RPLgenes$tif34.down.fold, Top64genes$tif34.down.fold)

linReg <- lm(y ~ x)

plot(x = NA, y = NA, 
	xlim = xlim, ylim = ylim, 
	xlab = "tif34.down.fold)", ylab = "tif35.up.fold)", 
	main = paste("All n=", nrow(All), sep = ""))
# abline(a = 0, b = 1, col = "gray")
abline(linReg, col = "gray", lwd = 2)
points(x = NoRPno64genes$tif34.down.fold, y = NoRPno64genes$tif35.up.fold, pch = 20, 
	col = gray(0, alpha = 0.2))
points(x = RPSgenes$tif34.down.fold, y = RPSgenes$tif35.up.fold, pch = 20, 
	col = rgb(1, 0, 0, alpha = 0.5))
points(x = RPLgenes$tif34.down.fold, y = RPLgenes$tif35.up.fold, pch = 20, 
	col = rgb(0, 0.8, 0.2, alpha = 0.5))
points(x = Top64genes$tif34.down.fold, y = Top64genes$tif35.up.fold, pch = 20, 
	col = rgb(0, 0, 1, alpha = 0.5))

