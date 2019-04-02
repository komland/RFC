rm(list = ls())

library(MASS)
library(lattice)

gr <- (sqrt(5) + 1)/2

## identify the latest RDS data file available
(availableFiles <- sort(grep("dat4s.RDS",
                             list.files(), value = TRUE)))
dat4s <- readRDS(file = availableFiles[length(availableFiles)])

## late revision - see analysis_spinach_20190204.r
dat4s$omit <- dat4s$VARIETY %in% c("Malabar", "Molokhia (egyptian spinach)",
                                   "Mustard Spinach", "New Zealand",
                                   "Sioux or Aztec")
xtabs(~ VARIETY + omit, dat4s, addNA = TRUE)
xtabs(~ omit, dat4s, addNA = TRUE)
nrow(dat4s <- dat4s[!dat4s$omit,])

## no longer need to deal with the comma-separated strings!
set3 <- grep("soil3.xrf", names(dat4s), value = TRUE)
set3[nchar(set3) >= 13]
set3 <- set3[nchar(set3) < 13]
(set3 <- sub("soil3\\.xrf\\.", "", set3))

set9 <- grep("soil9.xrf", names(dat4s), value = TRUE)
set9[nchar(set9) >= 13]
set9 <- set9[nchar(set9) < 13]
(set9 <- sub("soil9\\.xrf\\.", "", set9))

setS <- grep("spinach.xrf", names(dat4s), value = TRUE)
setS[nchar(setS) >= 15]
setS <- setS[nchar(setS) < 15]
(setS <- sub("spinach\\.xrf\\.", "", setS)) # a shorter list

identical(set3, set9)
setdiff(set3, setS)
setdiff(setS, set3)

## consensus list
setA <- intersect(set3, setS)

## Rhodium 0 in all soil and spinach samples
sort(unique(matrix(as.matrix(
    dat4s[,c("soil9.xrf.Rh", "soil3.xrf.Rh", "spinach.xrf.Rh")]))))
setA <- setdiff(setA, c("Rh"))

## Sodium centered around 0 - has it been transformed already?
## view the eqscplots without log transformation
par(mfrow = c(1, 2), mar = c(4.1, 4.1, 0.6, 0.6))
i <- "Na"
## deep-shallow correlation
Y <- dat4s[paste(c("soil9.xrf.", "soil3.xrf."), i, sep = "")]
yDeep <- Y[1]; yShlw <- Y[2]
rDpSh <- as.numeric(cor(yDeep, yShlw, use = "complete.obs"))
axsLim <- range(c(yDeep, yShlw), na.rm = TRUE)
axsTcks <- seq(floor(axsLim[1]), ceiling(axsLim[2]), by = 1)
axsLim <- range(c(yDeep, yShlw, axsTcks), na.rm = TRUE)

eqscplot(x = as.numeric(yDeep[,1]), y = as.numeric(yShlw[,1]),
         xlim = axsLim, ylim = axsLim,
         xaxt = "n", xlab = paste(i, "(ppm), 15-30 cm"),
         yaxt = "n", ylab = paste(i, "(ppm), 0-15 cm"))
axis(1, at = axsTcks, labels = axsTcks)
abline(v = axsTcks, col = gray(0.9))
axis(2, at = axsTcks, labels = axsTcks)
abline(h = axsTcks, col = gray(0.9))
abline(0, 1, lty = 2)
text(x = axsLim[1], y = axsLim[2],
     labels = paste("r =", signif(rDpSh, digits = 2)),
     adj = c(-0.2, 1.2))
points(x = as.numeric(yDeep[,1]), y = as.numeric(yShlw[,1]), pch = 16)
box()

## then shallow-spinach correlation
Y <- dat4s[paste(c("soil3.xrf.", "spinach.xrf."), i, sep = "")]
yShlw <- Y[1]; ySpch <- Y[2]
rShCr <- as.numeric(cor(yShlw, ySpch, use = "complete.obs"))
## forget about eqscplot!
range(yShlw); range(ySpch)
plot(x = as.numeric(yShlw[,1]), y = as.numeric(ySpch[,1]),
     xlab = paste(i, "(ppm), 0-15 cm"),
     ylab = paste(i, "(ppm), Spinach"),
     pch = 16)
text(x = min(yShlw), y = max(ySpch),
     labels = paste("r =", signif(rShCr, digits = 2)),
     adj = c(-0.2, 1.2))

## so also exclude Na
setA <- setdiff(setA, c("Na"))

## limit to records with those elements
has3 <- apply(dat4s[,paste("soil3.xrf.", setA, sep = "")], 1, function(f){
    !all(is.na(f))})
unique(dat4s[!has3, paste("soil3.xrf.", setA, sep = "")])
has9 <- apply(dat4s[,paste("soil9.xrf.", setA, sep = "")], 1, function(f){
    !all(is.na(f))})
unique(dat4s[!has9, paste("soil9.xrf.", setA, sep = "")])
hasS <- apply(dat4s[,paste("spinach.xrf.", setA, sep = "")], 1, function(f){
    !all(is.na(f))})
unique(dat4s[!hasS, paste("spinach.xrf.", setA, sep = "")])
## n.b. when one element was missing, all were
data.frame(xtabs(~ hasS + has9 + has3))[,c(3:1,4)]
## 14 records with shallow, deep, and spinach
##  5 with shallow and spinach (but not deep)
## 96 with just spinach (which have not data for this analysis)
## 20 with none of the three (ditto)

## establish deep-shallow correlation using has9 and has3
xtabs(~ has9 + has3) # n = 14
nrow(dat4sSoil <- dat4s[has9 & has3,])
## study shallow-spinach correlation using has3 and hasS
xtabs(~ hasS + has3) # n = 14 + 5 = 19
nrow(dat4sSlCr <- dat4s[has3 & hasS,])

if(exists("smryElements")) rm(smryElements)
bgn <- Sys.time()
## beyond that, I think coercion to logs will just omit negative readings
for(i in setA){
    plotWidth <- 7.5
    ## png(filename = paste("output/corrSpch_", i, ".png", sep = ""),
    ##     width = plotWidth, height = plotWidth/2, units = "in", res = 600,
    ##     pointsize = 10)
    par(mfrow = c(1, 2), mar = c(4.1, 4.1, 0.6, 0.6))

    ## deep-shallow correlation
    Y <- dat4sSoil[paste(c("soil9.xrf.", "soil3.xrf."), i, sep = "")]
    Y <- Y[apply(Y, 1, function(f) !any(is.na(f)) & all(f > 0)),]
    ## if(i == "As") Y <- Y[Y[2] < 5000,]
    ## if(i == "Se") Y <- Y[Y[2] < 6000,]
    ## if(i == "Mo") Y <- Y[Y[2] < 199,]
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

    ## then shallow-spinach correlation
    Y <- dat4sSoil[paste(c("soil3.xrf.", "spinach.xrf."), i, sep = "")]
    ## kludge to fix column 2 for iron, which comes in as character, for some reason
    Y[,2] <- as.numeric(Y[,2])
    Y <- Y[apply(Y, 1, function(f) !any(is.na(f)) & all(f > 0)),]
    ## if(i == "As") Y <- Y[Y[1] < 5000,]
    ## if(i == "Ni") Y <- Y[Y[1] > 10,]
    ## if(i == "P") Y <- Y[Y[2] < 20000,]
    ## if(i == "Se") Y <- Y[Y[1] < 6000,]
    yShlw <- log2(Y[1]); ySpch <- log2(Y[2])
    rShCr <- as.numeric(cor(yShlw, ySpch, use = "complete.obs"))
    vSpch <- as.numeric(var(ySpch))
    mSpch <- as.numeric(sapply(ySpch, median))
    axsLim <- range(c(yShlw, ySpch), na.rm = TRUE)
    if(axsLim[2] - axsLim[1] < 2){
        axsCtr <- round(median(axsLim))
        axsLim <- axsCtr + c(-1,1)}
    axsTcks <- seq(floor(axsLim[1]), ceiling(axsLim[2]), by = 1)

    eqscplot(x = as.numeric(yShlw[,1]), y = as.numeric(ySpch[,1]),
             xlim = axsLim, ylim = axsLim,
             xaxt = "n", xlab = paste(i, "(ppm), 0-15 cm"),
             yaxt = "n", ylab = paste(i, "(ppm), Spinach"))
    axis(1, at = axsTcks, labels = 2^(axsTcks))
    abline(v = axsTcks, col = gray(0.9))
    axis(2, at = axsTcks, labels = 2^(axsTcks))
    abline(h = axsTcks, col = gray(0.9))
    abline(0, 1, lty = 2)
    text(x = axsLim[1], y = axsLim[2],
         labels = paste("r =", signif(rShCr, digits = 2)),
         adj = c(-0.2, 1.2))
    points(x = as.numeric(yShlw[,1]), y = as.numeric(ySpch[,1]), pch = 16)
    box()
    ## dev.off()

    dfInt <- data.frame(element = i,
                        mDeep = mDeep,
                        mShlw = mShlw,
                        mSpch = mSpch,
                        rDpSh = rDpSh,
                        rShCr = rShCr,
                        vDeep = vDeep,
                        vShlw = vShlw,
                        vSpch = vSpch)
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

smryElements$spchMratio <- exp(smryElements$mSpch - smryElements$mShlw)
smryElements$spchMrClss <- round(log10(smryElements$spchMratio))
## they're all about the same
smryElements[order(smryElements$spchMratio),
             c("element", "mSpch", "mShlw", "spchMratio", "spchMrClss")]
1/min(smryElements$spchMratio)
formatC(sort((smryElements$spchMratio)))
formatC(sort((1/smryElements$spchMratio)))
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

smryElements$spchVratio <- smryElements$vSpch/smryElements$vShlw
smryElements$spchVrClss <- ifelse(smryElements$spchVratio < 0.5, "less",
                           ifelse(smryElements$spchVratio > 2, "more",
                                  "similar"))
## they're wildly different!
smryElements[order(smryElements$spchVratio),
             c("element", "rShCr", "vShlw", "vSpch", "spchVratio", "spchVrClss")]
## 6 more than double: Al, Se, Ca, Si, Cu, P - none with strong correlation
## 5 less than half: Mg, Mo, Pb, Fe, K
## 5 similar: As, Mn, Zn, Ni, S

xtabs(~ spchMrClss + spchVrClss, smryElements)
## regulated and concentrated: Si, P
smryElements[smryElements$spchMrClss > 0 &
             smryElements$spchVrClss == "less",]
## regulated at about ambient level: Ca, Se
smryElements[smryElements$spchMrClss == 0 &
             smryElements$spchVrClss == "less",]
## regulated and eschewed or sequestered: Al, Cu
smryElements[smryElements$spchMrClss < 0 &
             smryElements$spchVrClss == "less",]

## unregulated: S
smryElements[smryElements$spchMrClss == 0 &
             smryElements$spchVrClss == "similar",]
## eschewed: Mn, Ni, Zn
smryElements[smryElements$spchMrClss < 0 &
             smryElements$spchVrClss == "similar",]
## concentrated proportionally
smryElements[smryElements$spchMrClss > 0 &
             smryElements$spchVrClss == "similar",]

## weird #1 - sometimes sequestered?: Mg, Fe, As, Pb, Mo
smryElements[smryElements$spchMrClss < 0 &
             smryElements$spchVrClss == "more",]
smryElements[smryElements$spchMrClss == 0 &
             smryElements$spchVrClss == "more",]
## concentrated willy-nilly: K
smryElements[smryElements$spchMrClss > 0 &
             smryElements$spchVrClss == "more",]
