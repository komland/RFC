rm(list = ls())

bgn0 <- Sys.time()

library(MASS)
library(mgcv)

## identify the latest RDS data file available
(availableFiles <- sort(grep("dat4c.RDS",
                             list.files(), value = TRUE)))
dat4c <- readRDS(file = availableFiles[length(availableFiles)])

## carrot quality
(crrtQVars <- grep("carrots.veg_results.data.", names(dat4c)))
names(dat4c)[crrtQVars]
names(dat4c) <- sub("carrots.veg_results.data", "crrtQ", names(dat4c))
names(dat4c)[crrtQVars]

## carrot scan
(crrtScanVars <- grep("carrotscan", names(dat4c)))
names(dat4c)[crrtScanVars]
names(dat4c) <- sub("carrots.carrotscan.median", "crrtScan", names(dat4c))
names(dat4c)[crrtScanVars]

## supernatant scan
(sprnScanVars <- grep("supernatant", names(dat4c)))
names(dat4c)[sprnScanVars]
names(dat4c) <- sub("carrots.supernatant.data.median", "sprnScan", names(dat4c))
names(dat4c)[sprnScanVars]

foo <- data.frame(
    xtabs(~ apply(dat4c[,crrtQVars], 1, function(f) all(is.finite(f))) +
              apply(dat4c[,crrtScanVars], 1, function(f) all(is.finite(f))) +
              apply(dat4c[,sprnScanVars], 1, function(f) all(is.finite(f)))))
names(foo) <- sub("...1..function.f..all.is.finite.f...", "",
                  sub("apply.dat4c...", "", names(foo)))
foo
## 3 with no data:
dat4c[apply(dat4c[,c(crrtQVars, crrtScanVars, sprnScanVars)], 1, function(f){
  all(!is.finite(f))}),
  c(crrtQVars, crrtScanVars, sprnScanVars)]
## only losing 85 partially populated records
## (25 with carrot scan but that's all, 9 with sprn scan but that's all)

## 491 with all 3 sets of variables fully populated
nrow(dat4c4 <- dat4c[apply(dat4c[,c(crrtQVars, crrtScanVars, sprnScanVars)],
                           1, function(f) all(is.finite(f))),])
## 415 of them orange-only
nrow(dat4cOrange <- dat4c4[!is.na(dat4c4$carrots.carrot_color) &
                             dat4c4$carrots.carrot_color == "orange",])

## PCA to condense the four crrtQVars into a single response
## look at distributions
sort(dat4cOrange$crrtQ.polyphenolsMgGae100gFw, na.last = TRUE)
sort(dat4cOrange$crrtQ.proteinMgPer100g, na.last = TRUE)
sort(dat4cOrange$crrtQ.antioxidentsFrap, na.last = TRUE)
sort(dat4cOrange$crrtQ.moistureContent, na.last = TRUE)

par(mfrow = c(2, 2))
qqnorm(dat4cOrange$crrtQ.polyphenolsMgGae100gFw)
qqline(dat4cOrange$crrtQ.polyphenolsMgGae100gFw)
qqnorm(dat4cOrange$crrtQ.proteinMgPer100g)
qqline(dat4cOrange$crrtQ.proteinMgPer100g)
qqnorm(dat4cOrange$crrtQ.antioxidentsFrap)
qqline(dat4cOrange$crrtQ.antioxidentsFrap)
qqnorm(dat4cOrange$crrtQ.moistureContent)
qqline(dat4cOrange$crrtQ.moistureContent)

## polyphenols = -208 is truly problematic; two other negative readings, too
dat4cOrange <- dat4cOrange[dat4cOrange$crrtQ.polyphenolsMgGae100gFw > 0,]
## return to line 51

## polyphenols, protein, antioxidants appear lognormal
qqnorm(log(dat4cOrange$crrtQ.polyphenolsMgGae100gFw))
qqline(log(dat4cOrange$crrtQ.polyphenolsMgGae100gFw))
qqnorm(log(dat4cOrange$crrtQ.proteinMgPer100g))
qqline(log(dat4cOrange$crrtQ.proteinMgPer100g))
qqnorm(log(dat4cOrange$crrtQ.antioxidentsFrap))
qqline(log(dat4cOrange$crrtQ.antioxidentsFrap))
## but not moisture ... but I'm not sure I'll use it anyway
qqnorm(dat4cOrange$crrtQ.moistureContent)
qqline(dat4cOrange$crrtQ.moistureContent)

grep("crrtQ\\.", names(dat4cOrange), value = TRUE)
crrtQ.fm <- princomp(
    log(as.matrix(unname(dat4cOrange[,grep("crrtQ\\.", names(dat4cOrange))]))),
    cor = TRUE)
summary(crrtQ.fm)
loadings(crrtQ.fm)
crrtQ.pc <- predict(crrtQ.fm)
eqscplot(crrtQ.pc[,1:2])
## two records high on PC1, very low on PC2:
dat4cOrange[crrtQ.pc[,1] > 5 & crrtQ.pc[,2] < -6, crrtQVars]
## they are the very low moisture samples

## use PC1 as indicator of carrot quality
dat4cOrange$crrtQ <- crrtQ.pc[,1]
## which, certainly, is a continuous indicator
qqnorm(dat4cOrange$crrtQ); qqline(dat4cOrange$crrtQ)
range(dat4cOrange$crrtQ)
hist(dat4cOrange$crrtQ, breaks = seq(-5, 6.5, 0.25))
## symmetrical/normal

## but let's get some other ones omitting moisture then protein
crrtQ.fm2 <- princomp(
    log(as.matrix(unname(dat4cOrange[,c("crrtQ.polyphenolsMgGae100gFw",
                                        "crrtQ.proteinMgPer100g",
                                        "crrtQ.antioxidentsFrap")]))),
    cor = TRUE)
summary(crrtQ.fm2)
loadings(crrtQ.fm2)
crrtQ.pc2 <- predict(crrtQ.fm2)
eqscplot(crrtQ.pc2[,1:2])
## use PC1 as indicator of carrot quality
dat4cOrange$crrtQ2 <- crrtQ.pc2[,1]
## which, certainly, is a continuous indicator
qqnorm(dat4cOrange$crrtQ2); qqline(dat4cOrange$crrtQ2)
range(dat4cOrange$crrtQ2)
hist(dat4cOrange$crrtQ2, breaks = seq(-5, 6.5, 0.25))
## quite nicely correlated with first estimate ... noting two outliers
plot(dat4cOrange$crrtQ, dat4cOrange$crrtQ2); abline(c(0, 1))

## and without protein
crrtQ.fm3 <- princomp(
    log(as.matrix(unname(dat4cOrange[,c("crrtQ.polyphenolsMgGae100gFw",
                                        "crrtQ.antioxidentsFrap")]))),
    cor = TRUE)
biplot(crrtQ.fm3)
summary(crrtQ.fm3)
loadings(crrtQ.fm3)
crrtQ.pc3 <- predict(crrtQ.fm3)
eqscplot(crrtQ.pc3[,1:2])
## use PC1 as indicator of carrot quality
dat4cOrange$crrtQ3 <- crrtQ.pc3[,1]
## which, certainly, is a continuous indicator
qqnorm(dat4cOrange$crrtQ3); qqline(dat4cOrange$crrtQ3)
range(dat4cOrange$crrtQ3)
hist(dat4cOrange$crrtQ3, breaks = seq(-5, 6.5, 0.25))
## quite nicely correlated with first estimate ... noting two outliers
plot(dat4cOrange$crrtQ, dat4cOrange$crrtQ2); abline(c(0, 1))
plot(dat4cOrange$crrtQ, dat4cOrange$crrtQ3); abline(c(0, 1))
plot(dat4cOrange$crrtQ2, dat4cOrange$crrtQ3); abline(c(0, 1))
## that case with about -4 with protein and about 0 without jumps out
dat4cOrange[dat4cOrange$crrtQ < -3 &
            dat4cOrange$crrtQ2 < - 3 &
            dat4cOrange$crrtQ3 > 0,
            c("crrtQ.polyphenolsMgGae100gFw",
              "crrtQ.proteinMgPer100g",
              "crrtQ.antioxidentsFrap")]
## the possible LLOQ on protein ... so I think even crrtQ2 is suspect
x <- log(dat4cOrange$crrtQ.antioxidentsFrap)
y <- log(dat4cOrange$crrtQ.polyphenolsMgGae100gFw)
range(exp(c(x, y)))
tcks <- c(seq(0.1, 0.5, 0.1), 0.7,
          seq(1, 5, 1), 7,
          seq(10, 50, 10), 70,
          seq(100, 500, 100), 700)

## png(filename = "output/qualityCorr.png",
##     width = 6, height = 6, units = "in", res = 600)
par(las = 1)
plot(x = x, xlab = "Antioxidants", xaxt = "n",
     y = y, ylab = "Polyphenols", yaxt = "n",
     abline(coef(lm(y ~ x))))
axis(1, at = log(tcks), labels = formatC(tcks))
axis(2, at = log(tcks), labels = formatC(tcks))
## dev.off()

## multiple regression to attempt to model quality by scan
summary(fm1 <- lm(crrtQ3 ~ 1, dat4cOrange))

fmCrrtScan <- stepAIC(fm1,
                  scope = list(upper = ~ crrtScan_365 *
                                   crrtScan_385 *
                                   crrtScan_450 *
                                   crrtScan_500 *
                                   crrtScan_587 *
                                   crrtScan_632 *
                                   crrtScan_850 *
                                   crrtScan_880 *
                                   crrtScan_940))
summary(fmCrrtScan)
## 450 (violet-blue) comes in first, then 880 (infrared), then 587 (yellow-orange)
## eventually also 632 (orange-red) and 850 (infrared)
## plus some interactions with 450
par(mfrow = c(2, 2))
plot(fmCrrtScan) # obs 329 with low crrtQ3 not influential

fmSprnScan <- stepAIC(fm1,
                      scope = list(upper = ~ sprnScan_365 *
                                   sprnScan_385 *
                                   sprnScan_450 *
                                   sprnScan_500 *
                                   sprnScan_587 *
                                   sprnScan_632 *
                                   sprnScan_850 *
                                   sprnScan_880 *
                                   sprnScan_940))
## 365 first (ultraviolet), then 500 (blue-green), then 587 (yellow-orange)
## also 632 (orange-red), 385 (ultraviolet), 850 and 880 (infrared)
## leaving out only 450 (which was strongest indicator in direct scans) and 940
## no interactions
summary(fmSprnScan)
plot(fmSprnScan) # obs 329 still not influential
dev.off()

## game: reserve ~ 10% of the data to try to predict G/F/P crrtQ3
0.9 * nrow(dat4cOrange)
## make it 372 to leave 40 to predict making nice decimals (multiples of 0.025)
traiN <- 372; nrow(dat4cOrange) - traiN

## quality classes - cutoffs at 1/4 and 1/2 quantiles of normal approximation of dist
## (noting mean(dat4cOrange$crrtQ3) = 0 due to PCA)
pCutoffs <- c(1/4, 1/2)
(qCutoffs <- qnorm(pCutoffs) * sd(dat4cOrange$crrtQ3))
QaxsPts <- qnorm(c(pCutoffs[1]/2, median(pCutoffs), (pCutoffs[2]+1)/2)) *
    sd(dat4cOrange$crrtQ3)
dat4cOrange$crrtQ3class <- factor(ifelse(dat4cOrange$crrtQ3 > qCutoffs[2], "Good",
                            ifelse(dat4cOrange$crrtQ3 > qCutoffs[1], "Fair",
                                   "Poor")), levels = c("Good", "Fair", "Poor"))
xtabs(~ crrtQ3class, dat4cOrange)
by(dat4cOrange, dat4cOrange$crrtQ3class, function(f) range(f$crrtQ3))

bgn1 <- Sys.time()
if(exists("smryCrrtScanConfTabl")){
    rm(smryCrrtScanConfTabl)}
#### begin loop
for(i in 1:10){#1000){
    ## training (traiN = 372) and test (n = 40) sets
    selectedRecords <- sample(1:nrow(dat4cOrange), size = traiN)
    selectedRecords <- dat4cOrange$Sample_ID[sort(selectedRecords)]
    dat4.tr <- dat4cOrange[dat4cOrange$Sample_ID %in% selectedRecords,]
    dat4.ts <- dat4cOrange[!dat4cOrange$Sample_ID %in% selectedRecords,]

    fm1 <- lm(crrtQ3 ~ 1, dat4.tr)
    fmCrrtScan <- stepAIC(fm1,
                          scope = list(upper = ~ crrtScan_365 *
                                           crrtScan_385 *
                                           crrtScan_450 *
                                           crrtScan_500 *
                                           crrtScan_587 *
                                           crrtScan_632 *
                                           crrtScan_850 *
                                           crrtScan_880 *
                                           crrtScan_940))
    dat4.ts$pred.crrtQ3 <- predict(fmCrrtScan, newdata = dat4.ts)

    dat4.ts$pred.crrtQ3class <-
        factor(ifelse(dat4.ts$pred.crrtQ3 > qCutoffs[2], "Good",
               ifelse(dat4.ts$pred.crrtQ3 > qCutoffs[1], "Fair",
                      "Poor")),
               levels = c("Good", "Fair", "Poor"))
    confMatr <- xtabs(~ crrtQ3class + pred.crrtQ3class, dat4.ts)
    confMatr <- data.frame(confMatr)
    confTabl <-
        data.frame(truPos = confMatr$Freq[confMatr$pred.crrtQ3class == "Good" &
                                          confMatr$crrtQ3class == "Good"],
                   truNeg = sum(confMatr$Freq[confMatr$pred.crrtQ3class != "Good" &
                                              confMatr$crrtQ3class != "Good"]),
                   flsPos = sum(confMatr$Freq[confMatr$pred.crrtQ3class == "Good" &
                                              confMatr$crrtQ3class != "Good"]),
                   flsNeg = sum(confMatr$Freq[confMatr$pred.crrtQ3class != "Good" &
                                              confMatr$crrtQ3class == "Good"]))

    ## tst <- sum(confTabl[grep("tru", names(confTabl))])
    ## if(tst >= 32 | tst <= 20){
    ##     print(paste("Hey! Check this one out!  ",
    ##                 paste(unname(confTabl), collapse = "-")))
    ##     png(filename = "output/currentImage.png",
    ##         width = 6, height = 6, units = "in", res = 600)
    ##     par(mar = c(4.1, 4.1, 0.6, 0.6))
    ##     eqscplot(x = dat4.ts$pred.crrtQ3, xlab = "Predicted Quality", xaxt = "n",
    ##              y = dat4.ts$crrtQ3, ylab = "Measured Quality", yaxt = "n")
    ##     abline(0, 1, lty = 2)
    ##     abline(v = qCutoffs)
    ##     abline(h = qCutoffs)
    ##     axis(1, at = QaxsPts, labels = rev(levels(dat4cOrange$crrtQ3class)),
    ##          tick = FALSE)
    ##     axis(2, at = QaxsPts, labels = rev(levels(dat4cOrange$crrtQ3class)),
    ##          tick = FALSE)
    ##     points(x = dat4.ts$pred.crrtQ3, y = dat4.ts$crrtQ3, pch = 16)
    ##     dev.off()
    ##     scan()} ## to enable renaming currentImage

    if(!exists("smryCrrtScanConfTabl")){
        smryCrrtScanConfTabl <- confTabl
    }else{
        smryCrrtScanConfTabl <- rbind(smryCrrtScanConfTabl, confTabl)
    }
}
#### end loop
ndg1 <- Sys.time()
ndg1 - bgn1

smryCrrtScanConfTabl
sapply(smryCrrtScanConfTabl, mean)/40

## same for supernatant
bgn2 <- Sys.time()
if(exists("smrySprnScanConfTabl")){
    rm(smrySprnScanConfTabl)}
#### begin loop
for(i in 1:10){#1000){
    ## training (traiN = 372) and test (n = 40) sets
    selectedRecords <- sample(1:nrow(dat4cOrange), size = traiN)
    selectedRecords <- dat4cOrange$Sample_ID[sort(selectedRecords)]
    dat4.tr <- dat4cOrange[dat4cOrange$Sample_ID %in% selectedRecords,]
    dat4.ts <- dat4cOrange[!dat4cOrange$Sample_ID %in% selectedRecords,]

    fm1 <- lm(crrtQ3 ~ 1, dat4.tr)
    fmSprnScan <- stepAIC(fm1,
                          scope = list(upper = ~ sprnScan_365 *
                                           sprnScan_385 *
                                           sprnScan_450 *
                                           sprnScan_500 *
                                           sprnScan_587 *
                                           sprnScan_632 *
                                           sprnScan_850 *
                                           sprnScan_880 *
                                           sprnScan_940))
    dat4.ts$pred.crrtQ3 <- predict(fmSprnScan, newdata = dat4.ts)

    dat4.ts$pred.crrtQ3class <-
        factor(ifelse(dat4.ts$pred.crrtQ3 > qCutoffs[2], "Good",
               ifelse(dat4.ts$pred.crrtQ3 > qCutoffs[1], "Fair",
                      "Poor")),
               levels = c("Good", "Fair", "Poor"))
    confMatr <- xtabs(~ crrtQ3class + pred.crrtQ3class, dat4.ts)
    confMatr <- data.frame(confMatr)
    confTabl <-
        data.frame(truPos = confMatr$Freq[confMatr$pred.crrtQ3class == "Good" &
                                          confMatr$crrtQ3class == "Good"],
                   truNeg = sum(confMatr$Freq[confMatr$pred.crrtQ3class != "Good" &
                                              confMatr$crrtQ3class != "Good"]),
                   flsPos = sum(confMatr$Freq[confMatr$pred.crrtQ3class == "Good" &
                                              confMatr$crrtQ3class != "Good"]),
                   flsNeg = sum(confMatr$Freq[confMatr$pred.crrtQ3class != "Good" &
                                              confMatr$crrtQ3class == "Good"]))

    ## tst <- sum(confTabl[grep("tru", names(confTabl))])
    ## if(tst >= 35 | tst <= 22){
    ##     print(paste("Hey! Check this one out!  ",
    ##                 paste(unname(confTabl), collapse = "-")))
    ##     png(filename = "output/currentImage.png",
    ##         width = 6, height = 6, units = "in", res = 600)
    ##     par(mar = c(4.1, 4.1, 0.6, 0.6))
    ##     eqscplot(x = dat4.ts$pred.crrtQ3, xlab = "Predicted Quality", xaxt = "n",
    ##              y = dat4.ts$crrtQ3, ylab = "Measured Quality", yaxt = "n")
    ##     abline(0, 1, lty = 2)
    ##     abline(v = qCutoffs)
    ##     abline(h = qCutoffs)
    ##     axis(1, at = QaxsPts, labels = rev(levels(dat4cOrange$crrtQ3class)),
    ##          tick = FALSE)
    ##     axis(2, at = QaxsPts, labels = rev(levels(dat4cOrange$crrtQ3class)),
    ##          tick = FALSE)
    ##     points(x = dat4.ts$pred.crrtQ3, y = dat4.ts$crrtQ3, pch = 16)
    ##     dev.off()
    ##     scan()} ## to enable renaming currentImage

    if(!exists("smrySprnScanConfTabl")){
        smrySprnScanConfTabl <- confTabl
    }else{
        smrySprnScanConfTabl <- rbind(smrySprnScanConfTabl, confTabl)
    }
}
#### end loop
ndg2 <- Sys.time()
ndg2 - bgn2

smrySprnScanConfTabl
sapply(smrySprnScanConfTabl, mean)/40

## ## more flexible spline models
## smCrrtScan <- gam(crrtQ3 ~ s(crrtScan_365) +
##                     s(crrtScan_385) +
##                     s(crrtScan_450) +
##                     s(crrtScan_500) +
##                     s(crrtScan_587) +
##                     s(crrtScan_632) +
##                     s(crrtScan_850) +
##                     s(crrtScan_880) +
##                     s(crrtScan_940),
##                   data = dat4cOrange)
## par(mfrow = c(3, 3))
## plot(smCrrtScan, residuals = TRUE)
## summary(smCrrtScan)

## ## surely that's overfitted
## ## backwards elimination, starting with 940
## smCrrtScan2 <- update(smCrrtScan, formula =  ~ s(crrtScan_365) +
##                                       s(crrtScan_385) +
##                                       s(crrtScan_450) +
##                                       s(crrtScan_500) +
##                                       s(crrtScan_587) +
##                                       s(crrtScan_632) +
##                                       s(crrtScan_850) +
##                                       s(crrtScan_880))
## summary(smCrrtScan2)
## anova(smCrrtScan, smCrrtScan2, test = "F")
## AIC(smCrrtScan, smCrrtScan2)
## par(mfrow = c(3, 3))
## plot(smCrrtScan2, residuals = TRUE)

## ## eliminate 365
## smCrrtScan3 <- update(smCrrtScan, formula =  ~ s(crrtScan_385) +
##                                       s(crrtScan_450) +
##                                       s(crrtScan_500) +
##                                       s(crrtScan_587) +
##                                       s(crrtScan_632) +
##                                       s(crrtScan_850) +
##                                       s(crrtScan_880))
## summary(smCrrtScan3)
## anova(smCrrtScan3, smCrrtScan2, test = "F")
## AIC(smCrrtScan, smCrrtScan2, smCrrtScan3)
## par(mfrow = c(3, 3))
## plot(smCrrtScan3, residuals = TRUE)
## ## AIC going up a little

## ## eliminate 632
## smCrrtScan4 <- update(smCrrtScan, formula =  ~ s(crrtScan_385) +
##                                       s(crrtScan_450) +
##                                       s(crrtScan_500) +
##                                       s(crrtScan_587) +
##                                       s(crrtScan_850) +
##                                       s(crrtScan_880))
## summary(smCrrtScan4)
## anova(smCrrtScan4, smCrrtScan3, test = "F")
## AIC(smCrrtScan, smCrrtScan2, smCrrtScan3, smCrrtScan4)
## ## AIC holding steady
## par(mfrow = c(2, 3))
## plot(smCrrtScan4, residuals = TRUE)

## ## eliminate 587
## smCrrtScan5 <- update(smCrrtScan, formula =  ~ s(crrtScan_385) +
##                                       s(crrtScan_450) +
##                                       s(crrtScan_500) +
##                                       s(crrtScan_850) +
##                                       s(crrtScan_880))
## summary(smCrrtScan5)
## anova(smCrrtScan5, smCrrtScan4, test = "F")
## AIC(smCrrtScan, smCrrtScan2, smCrrtScan3, smCrrtScan4, smCrrtScan5)
## ## AIC dropping
## par(mfrow = c(2, 3))
## plot(smCrrtScan5, residuals = TRUE)

## ## eliminate 850
## smCrrtScan6 <- update(smCrrtScan, formula =  ~ s(crrtScan_385) +
##                                       s(crrtScan_450) +
##                                       s(crrtScan_500) +
##                                       s(crrtScan_880))
## summary(smCrrtScan6)
## anova(smCrrtScan6, smCrrtScan5, test = "F")
## AIC(smCrrtScan, smCrrtScan2, smCrrtScan3, smCrrtScan4, smCrrtScan5, smCrrtScan6)
## ## AIC holding steady
## par(mfrow = c(2, 2))
## plot(smCrrtScan6, residuals = TRUE)

## ## eliminate 450
## smCrrtScan7 <- update(smCrrtScan, formula =  ~ s(crrtScan_385) +
##                                       s(crrtScan_500) +
##                                       s(crrtScan_880))
## summary(smCrrtScan7)
## anova(smCrrtScan7, smCrrtScan6, test = "F")
## AIC(smCrrtScan, smCrrtScan2, smCrrtScan3, smCrrtScan4, smCrrtScan5, smCrrtScan6, smCrrtScan7)
## ## 450 was useful

## ## 880 needn't be a smooth function
## smCrrtScan6b <- update(smCrrtScan, formula =  ~ s(crrtScan_385) +
##                                        s(crrtScan_450) +
##                                        s(crrtScan_500) +
##                                        crrtScan_880)
## summary(smCrrtScan6b)
## anova(smCrrtScan6b, smCrrtScan6, test = "F")
## AIC(smCrrtScan, smCrrtScan2, smCrrtScan3, smCrrtScan4, smCrrtScan5, smCrrtScan6, smCrrtScan6b, smCrrtScan7)
## ## latter just as good
## par(mfrow = c(2, 2))
## plot(smCrrtScan6b, residuals = TRUE)

## dat4cOrange$fitted.crrtQ3 <- smCrrtScan6b$fitted.values
## dat4cOrange$fitted.crrtQ3class <- factor(ifelse(dat4cOrange$fitted.crrtQ3 > qCutoffs[2], "Good",
##                                         ifelse(dat4cOrange$fitted.crrtQ3 > qCutoffs[1], "Fair",
##                                                "Poor")), levels = c("Good", "Fair", "Poor"))
## xtabs(~ crrtQ3class + fitted.crrtQ3class, dat4cOrange)

## par(mfrow = c(1, 1))
## eqscplot(smCrrtScan$fitted.values, dat4cOrange$crrtQ3)
## abline(h = 0, lty = 2); abline(v = 0, lty = 2); abline(c(0, 1))
## cor(smCrrtScan$fitted.values, dat4cOrange$crrtQ3)

## plot(smCrrtScan6b$fitted.values, smCrrtScan6b$residuals)
## abline(h = 0)

## ## supernatant
## smSprnScan <- gam(crrtQ3 ~ s(sprnScan_365) +
##                     s(sprnScan_385) +
##                     s(sprnScan_450) +
##                     s(sprnScan_500) +
##                     s(sprnScan_587) +
##                     s(sprnScan_632) +
##                     s(sprnScan_850) +
##                     s(sprnScan_880) +
##                     s(sprnScan_940),
##                   data = dat4cOrange)
## summary(smSprnScan)
## par(mfrow = c(3, 3))
## plot(smSprnScan, residuals = TRUE)

## ## eliminate 940
## smSprnScan2 <- update(smSprnScan, formula = crrtQ3 ~ s(sprnScan_365) +
##                                       s(sprnScan_385) +
##                                       s(sprnScan_450) +
##                                       s(sprnScan_500) +
##                                       s(sprnScan_587) +
##                                       s(sprnScan_632) +
##                                       s(sprnScan_850) +
##                                       s(sprnScan_880))
## summary(smSprnScan2)
## anova(smSprnScan2, smSprnScan, test = "F")
## AIC(smSprnScan2, smSprnScan)
## ## latter just as good
## par(mfrow = c(3, 3))
## plot(smSprnScan2, residuals = TRUE)

## ## eliminate 385
## smSprnScan3 <- update(smSprnScan, formula = crrtQ3 ~ s(sprnScan_365) +
##                                       s(sprnScan_450) +
##                                       s(sprnScan_500) +
##                                       s(sprnScan_587) +
##                                       s(sprnScan_632) +
##                                       s(sprnScan_850) +
##                                       s(sprnScan_880))
## summary(smSprnScan3)
## anova(smSprnScan3, smSprnScan2, test = "F")
## AIC(smSprnScan3, smSprnScan2, smSprnScan)
## ## just getting better
## par(mfrow = c(3, 3))
## plot(smSprnScan3, residuals = TRUE)

## ## eliminate 385
## smSprnScan3 <- update(smSprnScan, formula = crrtQ3 ~ s(sprnScan_365) +
##                                       s(sprnScan_450) +
##                                       s(sprnScan_500) +
##                                       s(sprnScan_587) +
##                                       s(sprnScan_632) +
##                                       s(sprnScan_850) +
##                                       s(sprnScan_880))
## summary(smSprnScan3)
## anova(smSprnScan3, smSprnScan2, test = "F")
## AIC(smSprnScan3, smSprnScan2, smSprnScan)
## ## just getting better
## par(mfrow = c(3, 3))
## plot(smSprnScan3, residuals = TRUE)




## dat4cOrange$fitt2d.crrtQ3 <- smSprnScan$fitted.values
## dat4cOrange$fitt2d.crrtQ3class <- factor(ifelse(dat4cOrange$fitt2d.crrtQ3 > qCutoffs[2], "Good",
##                                                ifelse(dat4cOrange$fitt2d.crrtQ3 > qCutoffs[1], "Fair",
##                                                       "Poor")), levels = c("Good", "Fair", "Poor"))
## xtabs(~ crrtQ3class + fitt2d.crrtQ3class, dat4cOrange)

## eqscplot(smSprnScan$fitted.values, dat4cOrange$crrtQ3)
## abline(h = 0, lty = 2); abline(v = 0, lty = 2)
## abline(c(0, 1))
## cor(smSprnScan$fitted.values, dat4cOrange$crrtQ3)

## plot(smSprnScan$fitted.values, smSprnScan$residuals)
## abline(h = 0)

## ## simulations using smooth models
## bgn3 <- Sys.time()
## ## if(exists("smr2CrrtScanConfTabl")){
## ##     rm(smr2CrrtScanConfTabl)}
## #### begin loop
## for(i in 1:1000){
##     ## training (traiN = 372) and test (n = 40) sets
##     selectedRecords <- sample(1:nrow(dat4cOrange), size = traiN)
##     selectedRecords <- dat4cOrange$Sample_ID[sort(selectedRecords)]
##     dat4.tr <- dat4cOrange[dat4cOrange$Sample_ID %in% selectedRecords,]
##     dat4.ts <- dat4cOrange[!dat4cOrange$Sample_ID %in% selectedRecords,]

##     smCrrtScan <- gam(crrtQ3 ~ s(crrtScan_365) +
##                           s(crrtScan_385) +
##                           s(crrtScan_450) +
##                           s(crrtScan_500) +
##                           s(crrtScan_587) +
##                           s(crrtScan_632) +
##                           s(crrtScan_850) +
##                           s(crrtScan_880) +
##                           s(crrtScan_940),
##                       data = dat4.tr)
##     dat4.ts$pred.crrtQ3 <- predict(smCrrtScan, newdata = dat4.ts)

##     dat4.ts$pred.crrtQ3class <-
##         factor(ifelse(dat4.ts$pred.crrtQ3 > qCutoffs[2], "Good",
##                ifelse(dat4.ts$pred.crrtQ3 > qCutoffs[1], "Fair",
##                       "Poor")),
##                levels = c("Good", "Fair", "Poor"))
##     confMatr <- xtabs(~ crrtQ3class + pred.crrtQ3class, dat4.ts)
##     confMatr <- data.frame(confMatr)
##     confTabl <-
##         data.frame(truPos = confMatr$Freq[confMatr$pred.crrtQ3class == "Good" &
##                                           confMatr$crrtQ3class == "Good"],
##                    truNeg = sum(confMatr$Freq[confMatr$pred.crrtQ3class != "Good" &
##                                               confMatr$crrtQ3class != "Good"]),
##                    flsPos = sum(confMatr$Freq[confMatr$pred.crrtQ3class == "Good" &
##                                               confMatr$crrtQ3class != "Good"]),
##                    flsNeg = sum(confMatr$Freq[confMatr$pred.crrtQ3class != "Good" &
##                                               confMatr$crrtQ3class == "Good"]))

##     if(!exists("smr2CrrtScanConfTabl")){
##         smr2CrrtScanConfTabl <- confTabl
##     }else{
##         smr2CrrtScanConfTabl <- rbind(smr2CrrtScanConfTabl, confTabl)
##     }
## }
## #### end loop
## ndg3 <- Sys.time()
## ndg3 - bgn3

## smr2CrrtScanConfTabl
## sapply(smr2CrrtScanConfTabl, mean)/40

## ## same for supernatant
## bgn4 <- Sys.time()
## ## if(exists("smr2SprnScanConfTabl")){
## ##     rm(smr2SprnScanConfTabl)}
## #### begin loop
## for(i in 1:1000){
##     ## training (traiN = 372) and test (n = 40) sets
##     selectedRecords <- sample(1:nrow(dat4cOrange), size = traiN)
##     selectedRecords <- dat4cOrange$Sample_ID[sort(selectedRecords)]
##     dat4.tr <- dat4cOrange[dat4cOrange$Sample_ID %in% selectedRecords,]
##     dat4.ts <- dat4cOrange[!dat4cOrange$Sample_ID %in% selectedRecords,]

##     smSprnScan <- gam(crrtQ3 ~ s(sprnScan_365) +
##                     s(sprnScan_385) +
##                     s(sprnScan_450) +
##                     s(sprnScan_500) +
##                     s(sprnScan_587) +
##                     s(sprnScan_632) +
##                     s(sprnScan_850) +
##                     s(sprnScan_880) +
##                     s(sprnScan_940),
##                   data =  dat4.tr)
##     dat4.ts$pred.crrtQ3 <- predict(smSprnScan, newdata = dat4.ts)

##     dat4.ts$pred.crrtQ3class <-
##         factor(ifelse(dat4.ts$pred.crrtQ3 > qCutoffs[2], "Good",
##                ifelse(dat4.ts$pred.crrtQ3 > qCutoffs[1], "Fair",
##                       "Poor")),
##                levels = c("Good", "Fair", "Poor"))
##     confMatr <- xtabs(~ crrtQ3class + pred.crrtQ3class, dat4.ts)
##     confMatr <- data.frame(confMatr)
##     confTabl <-
##         data.frame(truPos = confMatr$Freq[confMatr$pred.crrtQ3class == "Good" &
##                                           confMatr$crrtQ3class == "Good"],
##                    truNeg = sum(confMatr$Freq[confMatr$pred.crrtQ3class != "Good" &
##                                               confMatr$crrtQ3class != "Good"]),
##                    flsPos = sum(confMatr$Freq[confMatr$pred.crrtQ3class == "Good" &
##                                               confMatr$crrtQ3class != "Good"]),
##                    flsNeg = sum(confMatr$Freq[confMatr$pred.crrtQ3class != "Good" &
##                                               confMatr$crrtQ3class == "Good"]))

##     if(!exists("smr2SprnScanConfTabl")){
##         smr2SprnScanConfTabl <- confTabl
##     }else{
##         smr2SprnScanConfTabl <- rbind(smr2SprnScanConfTabl, confTabl)
##     }
## }
## #### end loop
## ndg4 <- Sys.time()
## ndg4 - bgn4

## smr2SprnScanConfTabl
## sapply(smr2SprnScanConfTabl, mean)/40

## include some metadata
## organic, whether certified or not
xtabs(~ organic + certified_organic, dat4cOrange)
## irrigated, whether scheduled or not
xtabs(~ irrigation + scheduled_irrigation, dat4cOrange)
## amendments, cover crops

## get fm1 back
summary(fm1 <- lm(crrtQ3 ~ 1, dat4cOrange))

fmCrrtScanM <- stepAIC(fm1,
                       scope = list(lower = ~ 1,
                                    upper = ~ crrtScan_365 *
                                        crrtScan_385 *
                                        crrtScan_450 *
                                        crrtScan_500 *
                                        crrtScan_587 *
                                        crrtScan_632 *
                                        crrtScan_850 *
                                        crrtScan_880 *
                                        crrtScan_940 *
                                        organic * irrigation *
                                        biological_amendments * cover_crops))
summary(fmCrrtScanM)

## in immitation of the other variables included in Dan's random forest models
## week, state, diameter, color ... all orange
boxplot(crrtQ3 ~ Weeks, data = dat4cOrange, varwidth = TRUE)
abline(h = median(dat4cOrange$crrtQ3))
abline(h = quantile(dat4cOrange$crrtQ3, probs = c(0.25, 0.75)), lty = 2)
## weeks 29, 31 were unusual
dat4cOr2 <- dat4cOrange[is.finite(dat4cOrange$Weeks) & dat4cOrange$Weeks >= 32,]
boxplot(crrtQ3 ~ Weeks, data = dat4cOr2, varwidth = TRUE)
abline(h = median(dat4cOr2$crrtQ3))
abline(h = quantile(dat4cOr2$crrtQ3, probs = c(0.25, 0.75)), lty = 2)
## week should be a random effect
summary(fm1b <- lm(crrtQ3 ~ 1, dat4cOr2))
fm1c <- lme(fixed = crrtQ3 ~ 1, data = dat4cOr2, random = ~ 1 | as.factor(Weeks),
            method = "ML")
summary(fm1c)
AIC(fm1b, fm1c)
VarCorr(fm1c)

## state
boxplot(crrtQ3 ~ STATE,
        data = dat4cOrange[!is.na(dat4cOrange$STATE) &
                           dat4cOrange$STATE %in% c("CT", "ME"),],
        varwidth = TRUE)

xtabs(~ STATE, dat4cOrange, addNA = TRUE)
prop.table(xtabs(~ STATE, dat4cOrange, addNA = TRUE))
## I don't like the lack of balance: > 40% for CT, another 36% from MI & NY
dat4cOr2 <- dat4cOrange[!is.na(dat4cOrange$STATE) & dat4cOrange$STATE %in% c("CT", "MI", "NY"),]
boxplot(crrtQ3 ~ STATE, data = dat4cOr2, varwidth = TRUE)
abline(h = median(dat4cOr2$crrtQ3))
abline(h = quantile(dat4cOr2$crrtQ3, probs = c(0.25, 0.75)), lty = 2)
## state should be a random effect
summary(fm1b <- lm(crrtQ3 ~ 1, dat4cOr2))
fm1c <- lme(fixed = crrtQ3 ~ 1, data = dat4cOr2, random = ~ 1 | as.factor(STATE),
            method = "ML")
summary(fm1c)
AIC(fm1b, fm1c)
11.961-7.272
VarCorr(fm1c)
0.05796113/(0.05796113+1.37039143)
ranef(fm1c)

## diameter as a potential covariate
nrow(dat4cOrange)
nrow(dat4cOr2 <- dat4cOrange[is.finite(dat4cOrange$carrots.carrot_diameter),])
summary(fm1b <- lm(crrtQ3 ~ 1, dat4cOr2))
fmCrrtScanM2 <- stepAIC(fm1b,
                       scope = list(lower = ~ 1,
                                    upper = ~ crrtScan_365 *
                                        crrtScan_385 *
                                        crrtScan_450 *
                                        crrtScan_500 *
                                        crrtScan_587 *
                                        crrtScan_632 *
                                        crrtScan_850 *
                                        crrtScan_880 *
                                        crrtScan_940))
fmCrrtScanM3 <- stepAIC(fmCrrtScanM2,
                        scope = list(lower = ~ 1,
                                     upper = ~ crrtScan_365 *
                                         crrtScan_385 *
                                         crrtScan_450 *
                                         crrtScan_500 *
                                         crrtScan_587 *
                                         crrtScan_632 *
                                         crrtScan_850 *
                                         crrtScan_880 *
                                         crrtScan_940 *
                                         carrots.carrot_diameter))
AIC(fmCrrtScanM2, fmCrrtScanM3)
anova(fmCrrtScanM2, fmCrrtScanM3)
summary(fmCrrtScanM2)
summary(fmCrrtScanM3)

fmSprnScanM2 <- stepAIC(fm1b,
                       scope = list(lower = ~ 1,
                                    upper = ~ sprnScan_365 *
                                        sprnScan_385 *
                                        sprnScan_450 *
                                        sprnScan_500 *
                                        sprnScan_587 *
                                        sprnScan_632 *
                                        sprnScan_850 *
                                        sprnScan_880 *
                                        sprnScan_940))
fmSprnScanM3 <- stepAIC(fmSprnScanM2,
                        scope = list(lower = ~ 1,
                                     upper = ~ sprnScan_365 *
                                         sprnScan_385 *
                                         sprnScan_450 *
                                         sprnScan_500 *
                                         sprnScan_587 *
                                         sprnScan_632 *
                                         sprnScan_850 *
                                         sprnScan_880 *
                                         sprnScan_940 *
                                         carrots.carrot_diameter))
AIC(fmSprnScanM2, fmSprnScanM3)
anova(fmSprnScanM2, fmSprnScanM3)
summary(fmSprnScanM2)
summary(fmSprnScanM3)

