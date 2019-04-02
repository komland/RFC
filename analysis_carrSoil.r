rm(list = ls())

library(MASS)
library(lattice)

gr <- (sqrt(5) + 1)/2

## identify the latest RDS data file available
(availableFiles <- sort(grep("dat4c.RDS",
                             list.files(), value = TRUE)))
dat4c <- readRDS(file = availableFiles[length(availableFiles)])

## no longer need to deal with the comma-separated strings!
set3 <- grep("soil3.xrf", names(dat4c), value = TRUE)
set3[nchar(set3) >= 13]
set3 <- set3[nchar(set3) < 13]
(set3 <- sub("soil3\\.xrf\\.", "", set3))

set9 <- grep("soil9.xrf", names(dat4c), value = TRUE)
set9[nchar(set9) >= 13]
set9 <- set9[nchar(set9) < 13]
(set9 <- sub("soil9\\.xrf\\.", "", set9))

setC <- grep("carrots.xrf", names(dat4c), value = TRUE)
setC[nchar(setC) >= 15]
setC <- setC[nchar(setC) < 15]
(setC <- sub("carrots\\.xrf\\.", "", setC)) # a shorter list

identical(set3, set9)
setdiff(set3, setC)
setdiff(setC, set3)

## consensus list
setA <- intersect(set3, setC)

## Rhodium 0 in all soil and carrot samples
foo <- data.frame(xtabs(~ soil9.xrf.Rh + soil3.xrf.Rh + carrots.xrf.Rh,
                        data = dat4c, addNA = TRUE))
foo <- foo[foo$Freq > 0,]
foo[order(foo$Freq, decreasing = TRUE),]
setA <- setdiff(setA, c("Rh"))

## limit to records with those elements
has3 <- apply(dat4c[,paste("soil3.xrf.", setA, sep = "")], 1, function(f){
    !all(is.na(f))})
unique(dat4c[!has3, paste("soil3.xrf.", setA, sep = "")])
has9 <- apply(dat4c[,paste("soil9.xrf.", setA, sep = "")], 1, function(f){
    !all(is.na(f))})
unique(dat4c[!has9, paste("soil9.xrf.", setA, sep = "")])
hasC <- apply(dat4c[,paste("carrots.xrf.", setA, sep = "")], 1, function(f){
    !all(is.na(f))})
unique(dat4c[!hasC, paste("carrots.xrf.", setA, sep = "")])
## n.b. when one element was missing, all were
data.frame(xtabs(~ hasC + has9 + has3))[,c(3:1,4)]
## 117 records with shallow, deep, and carrot
##  11 with shallow and carrot (but not deep)
## 396 with just carrot
##  18 with just soil, both shallow and deep

## establish deep-shallow correlation using has9 and has3
xtabs(~ has9 + has3) # n = 117 + 18 = 135
nrow(dat4cSoil <- dat4c[has9 & has3,])
## study shallow-carrot correlation using has3 and hasC
xtabs(~ hasC + has3) # n = 117 + 11 = 128
nrow(dat4cSlCr <- dat4c[has3 & hasC,])

for(i in setA){
    Y <- dat4cSoil[paste(c("soil9.xrf.", "soil3.xrf."), i, sep = "")]
    print(list(i, xtabs(~ (Y < 0) + (Y > 0))))}
## OK to use log transformation for these:
allPos <- c("Mg", "Al", "Si", "P", "S", "K", "Mn", "Fe", "Cu", "Zn")
setdiff(setA, allPos)

## four negative values for Ca
i <- "Ca"
Y <- dat4cSoil[paste(c("soil9.xrf.", "soil3.xrf."), i, sep = "")]
eqscplot(x = Y[,1], y = Y[,2], main = i); abline(h = 0); abline(v = 0)
Y[Y[1] <= 0 | Y[2] <= 0,]
## substitute half minimum positive value
min(Y[Y > 0])
## call it 100
Y[Y < 0] <- 100
eqscplot(x = log10(Y[,1]), y = log10(Y[,2]), main = i); abline(h = 0); abline(v = 0)

## seven negative values for Ni (one where both shallow and deep were negative)
i <- "Ni"
Y <- dat4cSoil[paste(c("soil9.xrf.", "soil3.xrf."), i, sep = "")]
eqscplot(x = Y[,1], y = Y[,2], main = i); abline(h = 0); abline(v = 0)
Y[Y[1] <= 0 | Y[2] <= 0,]
## substitute half minimum positive value
min(Y[Y > 0])
## call it 0.1
Y[Y < 0] <- 0.1
eqscplot(x = log10(Y[,1]), y = log10(Y[,2]), main = i); abline(h = 0); abline(v = 0)
## that's too violent; try ...
Y[Y < 1] <- 1
eqscplot(x = log10(Y[,1]), y = log10(Y[,2]), main = i); abline(h = 0); abline(v = 0)
## still too violent; try ...
Y[Y < 10] <- 10
eqscplot(x = log10(Y[,1]), y = log10(Y[,2]), main = i); abline(h = 0); abline(v = 0)

## two (extreme outlier) negative values for Pb
i <- "Pb"
Y <- dat4cSoil[paste(c("soil9.xrf.", "soil3.xrf."), i, sep = "")]
eqscplot(x = Y[,1], y = Y[,2], main = i); abline(h = 0); abline(v = 0)
Y[Y[1] <= 0 | Y[2] <= 0,]
## for Pb, I think the best approach is to omit those three records
omit <- apply(Y, 1, function(f) !all(is.finite(f) & f > 0))
Y[omit,]
Y <- Y[!omit,]
eqscplot(x = log10(Y[,1]), y = log10(Y[,2]), main = i); abline(h = 0); abline(v = 0)

## two NAs for each of As, Se, and Mo
i <- "As"
Y <- dat4cSoil[paste(c("soil9.xrf.", "soil3.xrf."), i, sep = "")]
eqscplot(x = Y[,1], y = Y[,2], main = i); abline(h = 0); abline(v = 0)
## here the real problem are the two 5000s
wickedHigh <- apply(Y, 1, function(f) !all(is.finite(f) & f < 5000))
eqscplot(x = Y[!wickedHigh,1], y = Y[!wickedHigh,2], main = i)

i <- "Se"
Y <- dat4cSoil[paste(c("soil9.xrf.", "soil3.xrf."), i, sep = "")]
eqscplot(x = Y[,1], y = Y[,2], main = i); abline(h = 0); abline(v = 0)
## again, real problem is wicked high outliers; set it up to see if they're the same
wickedHig2 <- apply(Y, 1, function(f) !all(is.finite(f) & f < 6000))
table(wickedHigh, wickedHig2) # yes, same readouts
eqscplot(x = Y[!wickedHigh,1], y = Y[!wickedHigh,2], main = i)

i <- "Mo"
Y <- dat4cSoil[paste(c("soil9.xrf.", "soil3.xrf."), i, sep = "")]
eqscplot(x = Y[,1], y = Y[,2], main = i); abline(h = 0); abline(v = 0)
## again ...
wickedHig3 <- apply(Y, 1, function(f) !all(is.finite(f) & f < 100))
table(wickedHigh, wickedHig3) # again, same readouts
eqscplot(x = Y[!wickedHigh,1], y = Y[!wickedHigh,2], main = i)

## are those problematic readouts problematic in other elements, too?
dat4cSoil[wickedHigh,
          c("Sample_ID", "xlLine",
            as.vector(outer(c("soil9.xrf.", "soil3.xrf."), setA, paste, sep="")))]
## yes, for Pb, anyway

## Na another story altogether: lots of negative values!
i <- "Na"
Y <- dat4cSoil[paste(c("soil9.xrf.", "soil3.xrf."), i, sep = "")]
eqscplot(x = Y[,1], y = Y[,2], main = i); abline(h = 0); abline(v = 0)
## let's just glance ahead at soil3 & carrot ...
Y <- dat4cSoil[paste(c("soil3.xrf.", "carrots.xrf."), i, sep = "")]
## not eqscplot!
plot(x = Y[,1], y = Y[,2], main = i, xlab = names(Y)[1], ylab = names(Y)[2])
abline(h = 0); abline(v = 0)
## which is just crazy: the familiar inflation of variance in the carrot ...
## but what does -8,000 ppm mean?
plot(x = Y[,1], y = Y[,2], main = i, xlab = names(Y)[1], ylab = names(Y)[2])
abline(h = 0); abline(v = 0)
## forget sodium
setA <- setdiff(setA, c("Na"))

if(exists("smryElements")) rm(smryElements)
bgn <- Sys.time()
## beyond that, I think coercion to logs will just omit negative readings
for(i in setA){

    plotWidth <- 7.5
    ## png(filename = paste("output/corrCarr_", i, ".png", sep = ""),
    ##     width = plotWidth, height = plotWidth/2, units = "in", res = 600,
    ##     pointsize = 10)
    par(mfrow = c(1, 2), mar = c(4.1, 4.1, 0.6, 0.6))

    ## deep-shallow correlation
    Y <- dat4cSoil[paste(c("soil9.xrf.", "soil3.xrf."), i, sep = "")]
    Y <- Y[apply(Y, 1, function(f) !any(is.na(f)) & all(f > 0)),]
    if(i == "As") Y <- Y[Y[2] < 5000,]
    if(i == "Se") Y <- Y[Y[2] < 6000,]
    if(i == "Mo") Y <- Y[Y[2] < 199,]
    yDeep <- log2(Y[1]); yShlw <- log2(Y[2])
    rDpSh <- as.numeric(cor(yDeep, yShlw, use = "complete.obs"))
    vDeep <- as.numeric(var(yDeep))
    vShlw <- as.numeric(var(yShlw))
    mDeep <- as.numeric(sapply(yDeep, median))
    mShlw <- as.numeric(sapply(yShlw, median))
    axsLim <- range(c(yDeep, yShlw), na.rm = TRUE)
    if(axsLim[2] - axsLim[1] < 2){
        axsCtr <- round(median(axsLim))
        axsLim <- axsCtr + c(-1,1)}
    axsTcks <- seq(floor(axsLim[1]), ceiling(axsLim[2]), by = 1)

    eqscplot(x = as.numeric(yDeep[,1]), y = as.numeric(yShlw[,1]),
             xlim = axsLim, ylim = axsLim,
             xaxt = "n", xlab = paste(i, "(ppm), 15-30 cm"),
             yaxt = "n", ylab = paste(i, "(ppm), 0-15 cm"))
    axis(1, at = axsTcks, labels = 2^(axsTcks))
    abline(v = axsTcks, col = gray(0.9))
    axis(2, at = axsTcks, labels = 2^(axsTcks))
    abline(h = axsTcks, col = gray(0.9))
    abline(0, 1, lty = 2)
    text(x = axsLim[1], y = axsLim[2],
         labels = paste("r =", signif(rDpSh, digits = 2)),
         adj = c(-0.2, 1.2))
    points(x = as.numeric(yDeep[,1]), y = as.numeric(yShlw[,1]), pch = 16)
    box()

    ## then shallow-carrot correlation
    Y <- dat4cSoil[paste(c("soil3.xrf.", "carrots.xrf."), i, sep = "")]
    Y <- Y[apply(Y, 1, function(f) !any(is.na(f)) & all(f > 0)),]
    if(i == "As") Y <- Y[Y[1] < 5000,]
    if(i == "Ni") Y <- Y[Y[1] > 10,]
    if(i == "P") Y <- Y[Y[2] < 20000,]
    if(i == "Se") Y <- Y[Y[1] < 6000,]
    yShlw <- log2(Y[1]); yCrrt <- log2(Y[2])
    rShCr <- as.numeric(cor(yShlw, yCrrt, use = "complete.obs"))
    vCrrt <- as.numeric(var(yCrrt))
    mCrrt <- as.numeric(sapply(yCrrt, median))
    axsLim <- range(c(yShlw, yCrrt), na.rm = TRUE)
    if(axsLim[2] - axsLim[1] < 2){
        axsCtr <- round(median(axsLim))
        axsLim <- axsCtr + c(-1,1)}
    axsTcks <- seq(floor(axsLim[1]), ceiling(axsLim[2]), by = 1)

    eqscplot(x = as.numeric(yShlw[,1]), y = as.numeric(yCrrt[,1]),
             xlim = axsLim, ylim = axsLim,
             xaxt = "n", xlab = paste(i, "(ppm), 0-15 cm"),
             yaxt = "n", ylab = paste(i, "(ppm), Carrot"))
    axis(1, at = axsTcks, labels = 2^(axsTcks))
    abline(v = axsTcks, col = gray(0.9))
    axis(2, at = axsTcks, labels = 2^(axsTcks))
    abline(h = axsTcks, col = gray(0.9))
    abline(0, 1, lty = 2)
    text(x = axsLim[1], y = axsLim[2],
         labels = paste("r =", signif(rShCr, digits = 2)),
         adj = c(-0.2, 1.2))
    points(x = as.numeric(yShlw[,1]), y = as.numeric(yCrrt[,1]), pch = 16)
    box()
    ## dev.off()

    dfInt <- data.frame(element = i,
                        mDeep = mDeep,
                        mShlw = mShlw,
                        mCrrt = mCrrt,
                        rDpSh = rDpSh,
                        rShCr = rShCr,
                        vDeep = vDeep,
                        vShlw = vShlw,
                        vCrrt = vCrrt)
    if(!exists("smryElements")){
        smryElements <- dfInt
    }else{
        smryElements <- rbind(smryElements, dfInt)
    }
}
Sys.time() - bgn

smryElements$soilMratio <- exp(smryElements$mShlw - smryElements$mDeep)
## they're all about the same
smryElements[order(smryElements$soilMratio),
             c("element", "mShlw", "mDeep", "soilMratio")]
1/min(smryElements$soilMratio)
## phosphorus, calcium, lead, sulfur augmented

smryElements$crrtMratio <- exp(smryElements$mCrrt - smryElements$mShlw)
smryElements$crrtMrClss <- round(log10(smryElements$crrtMratio))
smryElements[order(smryElements$crrtMratio),
             c("element", "mCrrt", "mShlw", "crrtMratio", "crrtMrClss")]
1/min(smryElements$crrtMratio)
formatC(sort((smryElements$crrtMratio)))
formatC(sort((1/smryElements$crrtMratio)))
## concentrated 96x! Si
## concentrated 3-4x: K, P
## eschewed by factors of about 10,000! Fe, Al
## eschewed by factors of about 1,000: Pb, As
## eschewed by factors of about 100: Mo, Mn, Cu
## eschewed by factors of about 10: Ni, Zn
## held at about soil concentration: Mg, Se, Ca, S

smryElements$soilVratio <- smryElements$vShlw/smryElements$vDeep
## they're all about the same
smryElements[order(smryElements$soilVratio),
             c("element", "rDpSh", "vDeep", "vShlw", "soilVratio")]
1/min(smryElements$soilVratio)

smryElements$crrtVratio <- smryElements$vCrrt/smryElements$vShlw
smryElements$crrtVrClss <- ifelse(smryElements$crrtVratio < 0.5, "less",
                           ifelse(smryElements$crrtVratio > 2, "more",
                                  "similar"))
## they're wildly different!
smryElements[order(smryElements$crrtVratio),
             c("element", "rShCr", "vShlw", "vCrrt", "crrtVratio", "crrtVrClss")]
## 6 less than half: Al, Se, Ca, Si, Cu, P
## 5 more than double: Mg, Mo, Pb, Fe, K
## 5 similar: As, Mn, Zn, Ni, S

xtabs(~ crrtMrClss + crrtVrClss, smryElements)
## regulated and concentrated: Si, P
smryElements[smryElements$crrtMrClss > 0 &
             smryElements$crrtVrClss == "less",]
## regulated at about ambient level: Ca, Se
smryElements[smryElements$crrtMrClss == 0 &
             smryElements$crrtVrClss == "less",]
## regulated and eschewed or sequestered: Al, Cu
smryElements[smryElements$crrtMrClss < 0 &
             smryElements$crrtVrClss == "less",]

## unregulated: S
smryElements[smryElements$crrtMrClss == 0 &
             smryElements$crrtVrClss == "similar",]
## eschewed: Mn, Ni, Zn, As
smryElements[smryElements$crrtMrClss < 0 &
             smryElements$crrtVrClss == "similar",]
## concentrated proportionally (an empty set)
smryElements[smryElements$crrtMrClss > 0 &
             smryElements$crrtVrClss == "similar",]

## weird #1 - sometimes sequestered?: Fe, Pb, Mo, ...
smryElements[smryElements$crrtMrClss < 0 &
             smryElements$crrtVrClss == "more",]
## ... Mg
smryElements[smryElements$crrtMrClss == 0 &
             smryElements$crrtVrClss == "more",]
## concentrated willy-nilly: K
smryElements[smryElements$crrtMrClss > 0 &
             smryElements$crrtVrClss == "more",]
