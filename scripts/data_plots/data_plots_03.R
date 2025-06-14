## original note: 20250614_AsanoTE_3.txt


library(openxlsx)

HOMEDIR <- ""
setwd(HOMEDIR)
setwd("data/data_plots")
load(file = "20250614_NewDFetc.RData")

rm(list = c("All", "OtherGenes", "RPLgenes", "RPSgenes", "tHSR"))

All <- subset(NewDF, is.na(tif34.down.fold) == FALSE & 
		is.na(NewDF[["RA40/0"]]) == FALSE & 
		NewDF[["RA40/0"]] > 0)
tHSR <- subset(All, tHSR == TRUE)
OtherGenes <- subset(All, RPS == FALSE & RPL == FALSE & tHSR == FALSE) 

RPSgenes <- subset(All, RPS == TRUE)
RPLgenes <- subset(All, RPL == TRUE)

save(list = c("NewDF", "All", "tHSR", 
		"OtherGenes", "RPSgenes", "RPLgenes"), file = "20250614_NewDFetc.RData")


## tHSR RPS RPL
## 
xmin <- log10(min(OtherGenes$tif34.down.fold, tHSR$tif34.down.fold, RPSgenes$tif34.down.fold, RPLgenes$tif34.down.fold))
xmax <- log10(max(OtherGenes$tif34.down.fold, tHSR$tif34.down.fold, RPSgenes$tif34.down.fold, RPLgenes$tif34.down.fold))
xlim <- c(xmin, xmax)

ymin <- log10(min(OtherGenes[["RA40/0"]], tHSR[["RA40/0"]], RPSgenes[["RA40/0"]], RPLgenes[["RA40/0"]]))
ymax <- log10(max(OtherGenes[["RA40/0"]], tHSR[["RA40/0"]], RPSgenes[["RA40/0"]], RPLgenes[["RA40/0"]]))
ylim <- c(ymin, ymax)

plot(x = log10(OtherGenes$tif34.down.fold), y = log10(OtherGenes[["RA40/0"]]), 
	xlab = "log10(tif34.down.fold)", ylab = "log10(RA40/0)", pch = 20,
	col = gray(0, alpha = 0.2), xlim = xlim, ylim = ylim, 
	main = paste("OtherGenes n=", nrow(OtherGenes), sep = ""))

plot(x = log10(tHSR$tif34.down.fold), y = log10(tHSR[["RA40/0"]]), 
	xlab = "log10(tif34.down.fold)", ylab = "log10(RA40/0)", pch = 20,
	col = rgb(0, 0, 1, alpha = 0.5), xlim = xlim, ylim = ylim, 
	main = paste("tHSR n=", nrow(tHSR), sep = ""))

plot(x = log10(RPSgenes$tif34.down.fold), y = log10(RPSgenes[["RA40/0"]]), 
	xlab = "log10(tif34.down.fold)", ylab = "log10(RA40/0)", pch = 20,
	col = rgb(1, 0, 0, alpha = 0.5), xlim = xlim, ylim = ylim, 
	main = paste("RPSgenes n=", nrow(RPSgenes), sep = ""))

plot(x = log10(RPLgenes$tif34.down.fold), y = log10(RPLgenes[["RA40/0"]]), 
	xlab = "log10(tif34.down.fold)", ylab = "log10(RA40/0)", pch = 20,
	col = rgb(0, 0.8, 0.2, alpha = 0.5), xlim = xlim, ylim = ylim, 
	main = paste("RPLgenes n=", nrow(RPLgenes), sep = ""))

plot(x = log10(All$tif34.down.fold), y = log10(All[["RA40/0"]]), 
	xlab = "log10(tif34.down.fold)", ylab = "log10(RA40/0)", pch = 20,
	col = gray(0, alpha = 0.2), xlim = xlim, ylim = ylim, 
	main = paste("All n=", nrow(All), sep = ""))

plot(x = log10(OtherGenes$tif34.down.fold), y = log10(OtherGenes[["RA40/0"]]), 
	xlab = "log10(tif34.down.fold)", ylab = "log10(RA40/0)", pch = 20,
	col = gray(0, alpha = 0.2), xlim = xlim, ylim = ylim, 
	main = paste("All n=", nrow(All), sep = ""))
points(x = log10(tHSR$tif34.down.fold), y = log10(tHSR[["RA40/0"]]), pch = 20,
	col = rgb(0, 0, 1, alpha = 0.5))
points(x = log10(RPSgenes$tif34.down.fold), y = log10(RPSgenes[["RA40/0"]]), pch = 20,
	col = rgb(1, 0, 0, alpha = 0.5))
points(x = log10(RPLgenes$tif34.down.fold), y = log10(RPLgenes[["RA40/0"]]), pch = 20,
	col = rgb(0, 0.8, 0.2, alpha = 0.5))

temp <- data.frame(tif34.down.fold = All$tif34.down.fold, RPKM.C = All[["RA40/0"]])
temp.lm <- lm(log10(RPKM.C) ~ log10(tif34.down.fold), data = temp)
summary(temp.lm)

# > summary(temp.lm)
# 
# Call:
# lm(formula = log10(RPKM.C) ~ log10(tif34.down.fold), data = temp)
# 
# Residuals:
#      Min       1Q   Median       3Q      Max 
# -1.61229 -0.21022  0.02978  0.23586  2.37100 
# 
# Coefficients:
#                         Estimate Std. Error t value Pr(>|t|)    
# (Intercept)             0.125955   0.007008   17.97   <2e-16 ***
# log10(tif34.down.fold) -0.349717   0.032864  -10.64   <2e-16 ***
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.4093 on 4459 degrees of freedom
# Multiple R-squared:  0.02477,	Adjusted R-squared:  0.02455 
# F-statistic: 113.2 on 1 and 4459 DF,  p-value: < 2.2e-16

# intercept: 0.125955
# slope: -0.349717
# R^2: 0.02455

plot(x = log10(OtherGenes$tif34.down.fold), y = log10(OtherGenes[["RA40/0"]]), 
	xlab = "log10(tif34.down.fold)", ylab = "log10(RA40/0)", pch = 20,
	col = gray(0, alpha = 0.2), xlim = xlim, ylim = ylim, 
	main = paste("All n=", nrow(All), sep = ""))
points(x = log10(tHSR$tif34.down.fold), y = log10(tHSR[["RA40/0"]]), pch = 20,
	col = rgb(0, 0, 1, alpha = 0.5))
points(x = log10(RPSgenes$tif34.down.fold), y = log10(RPSgenes[["RA40/0"]]), pch = 20,
	col = rgb(1, 0, 0, alpha = 0.5))
points(x = log10(RPLgenes$tif34.down.fold), y = log10(RPLgenes[["RA40/0"]]), pch = 20,
	col = rgb(0, 0.8, 0.2, alpha = 0.5))
abline(temp.lm)

