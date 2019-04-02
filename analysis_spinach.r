rm(list = ls())

bgn0 <- Sys.time()

library(MASS)
library(mgcv)
library(tree)
library(randomForest)

gr <- (1 + sqrt(5))/2

## identify the latest RDS data file available
(availableFiles <- sort(grep("dat4s.RDS",
                             list.files(), value = TRUE)))
dat4s <- readRDS(file = availableFiles[length(availableFiles)])

## late revision: some of these "varieties" are far from Spinacia oleracea!
data.frame(xtabs(~ VARIETY, dat4s))
## omit:
##   Malabar (perennial vine, Basella alba, Basellaceae)
##   Molokhia (shrub, Corchorus olitorius, Malvaceae)
##   New Zealand (Tetragonia tetragonoides, Aizoaceae)
##   Mustard (Brassica rapa var. perviridis, Brassicaceae)
##   Sious or Aztec (Chenopodium nuttalliae, Chenopodioideae)
## retain:
##   "baby" - but be on guard for associated patterns
##   Emperor, Flamingo (https://www.johnnyseeds.com/growers-library/vegetables/spinachprogram.html)
dat4s$omit <- dat4s$VARIETY %in% c("Malabar", "Molokhia (egyptian spinach)",
                                   "Mustard Spinach", "New Zealand",
                                   "Sioux or Aztec")
xtabs(~ VARIETY + omit, dat4s, addNA = TRUE)

## scatterplots of polyphenols and protein over antioxidants
par(mfrow = c(1, 2))
plot(x = log(dat4s$spinach.veg_results.data.antioxidentsFrap),
     y = log(dat4s$spinach.veg_results.data.polyphenolsMgGae100gFw))
points(x = log(dat4s$spinach.veg_results.data.antioxidentsFrap[dat4s$omit]),
       y = log(dat4s$spinach.veg_results.data.polyphenolsMgGae100gFw[dat4s$omit]),
       pch = 8, col = "red")
## two in upper right are "Molokhia (egyptian spinach)" and will be omitted
plot(x = log(dat4s$spinach.veg_results.data.antioxidentsFrap),
     y = dat4s$spinach.veg_results.data.proteinMgPer100g)
points(x = log(dat4s$spinach.veg_results.data.antioxidentsFrap[dat4s$omit]),
       y = dat4s$spinach.veg_results.data.proteinMgPer100g[dat4s$omit],
       pch = 8, col = "red")
xtabs(~ VARIETY + (dat4s$spinach.veg_results.data.proteinMgPer100g > 10),
      dat4s, addNA = TRUE)
xtabs(~ omit + (dat4s$spinach.veg_results.data.proteinMgPer100g > 10),
      dat4s, addNA = TRUE)

## png(filename = "output/spchProtein.png",
##     width = 3 * gr, height = 3, units = "in", res = 600)
par(mar = c(5.1, 4.1, 1.1, 1.1))
hist(dat4s$spinach.veg_results.data.proteinMgPer100g, breaks = 0:27, main = "",
     xlab = "Protein (Mg per 100g)",
     yaxs = "i", ylim = c(0, 40),
     las = 1)
box()
## dev.off()

## spinach quality - polyphenols and antioxidants only
(spchQVars <- grep("spinach.veg_results.data.", names(dat4s)))
names(dat4s)[spchQVars]
spchQVars <- spchQVars[1:2]
names(dat4s) <- sub("spinach.veg_results.data", "spchQ", names(dat4s))
names(dat4s)[spchQVars]

## spinach scan
(leafScanVars <- grep("spinachscan", names(dat4s)))
names(dat4s)[leafScanVars]
names(dat4s) <- sub("spinach.spinachscan.median", "leafScan", names(dat4s))
names(dat4s)[leafScanVars]

## supernatant scan
(xtrcScanVars <- grep("supernatant", names(dat4s)))
names(dat4s)[xtrcScanVars]
names(dat4s) <- sub("spinach.supernatant.data.median", "xtrcScan", names(dat4s))
names(dat4s)[xtrcScanVars]

foo <- data.frame(
    xtabs(~ dat4s$omit +
              apply(dat4s[,spchQVars], 1, function(f) all(is.finite(f))) +
              apply(dat4s[,leafScanVars], 1, function(f) all(is.finite(f))) +
              apply(dat4s[,xtrcScanVars], 1, function(f) all(is.finite(f)))))
names(foo) <- sub("Vars...1..function.f..all.is.finite.f...", "",
                  sub("dat4s.", "",
                      sub("apply.dat4s...", "", names(foo))))
foo <- foo[foo$Freq > 0,]
foo[order(foo$omit, foo$Freq),]
## will use 94 records (not non-spinach and have quality and scan data)
## losing a grand total of 59 records:
## 18 because they aren't spinach - 5 lack other information, too
##  8 because they don't have antioxidant _and_ polyphenols (2 have one or other)
## 20 because all they have is leaf scan data
## 13 because they're missing quality and both scans
nrow(dat4s4 <- dat4s[!dat4s$omit &
                     apply(dat4s[,c(spchQVars, leafScanVars, xtrcScanVars)],
                           1, function(f) all(is.finite(f))),])

## identify three classes based on antioxidants >= 300 and polyphenols >= 30:
dat4s4$Qclass <- factor(ifelse(dat4s4$spchQ.antioxidentsFrap >= 300 &
                               dat4s4$spchQ.polyphenolsMgGae100gFw >= 30, "highQ",
                        ifelse(dat4s4$spchQ.antioxidentsFrap < 300 &
                               dat4s4$spchQ.polyphenolsMgGae100gFw >= 30, "lowAO",
                        ifelse(dat4s4$spchQ.antioxidentsFrap >= 300 &
                               dat4s4$spchQ.polyphenolsMgGae100gFw < 30, "lowPP",
                               NA))),
                        levels = c("highQ", "lowAO", "lowPP"))
xtabs(~ Qclass, dat4s4, addNA = TRUE)
xtabs(~ VARIETY + Qclass, dat4s4, addNA = TRUE)
by(dat4s4, dat4s4$Qclass, function(f){
    return(list(aoRange = range(f$spchQ.antioxidentsFrap),
                ppRange = range(f$spchQ.polyphenolsMgGae100gFw)))})

tcks <- c(seq(1, 5, 1), 7,
          seq(10, 50, 10), 70,
          seq(100, 500, 100), 700,
          seq(1000, 5000, 1000))

## png(filename = "output/spchQCor2.png",
##     width = 6, height = 6, units = "in", res = 600)
par(las = 1, mar = c(5.1, 4.1, 1.1, 1.1))
eqscplot(x = log(dat4s4$spchQ.antioxidentsFrap), xlab = "Antioxidants", xaxt = "n",
     y = log(dat4s4$spchQ.polyphenolsMgGae100gFw), ylab = "Polyphenols", yaxt = "n",
     type = "n")
axis(1, at = log(tcks), labels = formatC(tcks))
axis(2, at = log(tcks), labels = formatC(tcks))
abline(v = log(c(300, 1000)), lty = 2)
abline(h = log(c(30, 100)), lty = 2)
points(x = log(dat4s4$spchQ.antioxidentsFrap[dat4s4$Qclass == "highQ"]),
       y = log(dat4s4$spchQ.polyphenolsMgGae100gFw[dat4s4$Qclass == "highQ"]),
       pch = 1)
points(x = log(dat4s4$spchQ.antioxidentsFrap[dat4s4$Qclass == "lowAO"]),
       y = log(dat4s4$spchQ.polyphenolsMgGae100gFw[dat4s4$Qclass == "lowAO"]),
       pch = 2)
points(x = log(dat4s4$spchQ.antioxidentsFrap[dat4s4$Qclass == "lowPP"]),
       y = log(dat4s4$spchQ.polyphenolsMgGae100gFw[dat4s4$Qclass == "lowPP"]),
       pch = 6)
## dev.off()

## now: can those classes be identified by spectral scans?
## leaves
leafScanTree <- tree(Qclass ~ leafScan_365 +
                         leafScan_385 +
                         leafScan_450 +
                         leafScan_500 +
                         leafScan_587 +
                         leafScan_632 +
                         leafScan_850 +
                         leafScan_880 +
                         leafScan_940,
                     data = dat4s4)
summary(leafScanTree)
plot(leafScanTree)
text(leafScanTree, all = TRUE, cex = 0.5)

dev.off()

## png(filename = "output/leafScanTree.png",
##     width = 6, height = 6, units = "in", res = 600)
tree.screens()
plot(leafScanTree)
text(leafScanTree, all = TRUE, cex = 0.5)
tile.tree(leafScanTree, dat4s4$Qclass)
close.screen(all = TRUE)
## dev.off()

leafScanTre2 <- prune.tree(leafScanTree)
plot(leafScanTre2, order = "increasing")

## extracts
xtrcScanTree <- tree(Qclass ~ xtrcScan_365 +
                         xtrcScan_385 +
                         xtrcScan_450 +
                         xtrcScan_500 +
                         xtrcScan_587 +
                         xtrcScan_632 +
                         xtrcScan_850 +
                         xtrcScan_880 +
                         xtrcScan_940,
                     data = dat4s4)
summary(xtrcScanTree)
plot(xtrcScanTree)
text(xtrcScanTree, all = TRUE, cex = 0.5)

dev.off()

## png(filename = "output/xtrcScanTree.png",
##     width = 6, height = 6, units = "in", res = 600)
tree.screens()
plot(xtrcScanTree)
text(xtrcScanTree, all = TRUE, cex = 0.5)
tile.tree(xtrcScanTree, dat4s4$Qclass)
close.screen(all = TRUE)
## dev.off()

xtrcScanTre2 <- prune.tree(xtrcScanTree)
plot(xtrcScanTre2, order = "increasing")

## scoring simulations
## reserve ~ 10% of the data to try to predict G/F/P spchQ
c(0.1, 0.9) * nrow(dat4s4)
## make it 84 to leave 10 to predict making nice decimals
(traiN <- 84); nrow(dat4s4) - traiN

bgn1 <- Sys.time()
if(exists("smryLeafScanScores")){
    rm(smryLeafScanScores)}
#### begin loop
for(i in 1:10){#1000
    ## training (traiN = 84) and test (n = 10) sets
    selectedRecords <- sample(1:nrow(dat4s4), size = traiN)
    selectedRecords <- dat4s4$Sample_ID[sort(selectedRecords)]
    dat4.tr <- dat4s4[dat4s4$Sample_ID %in% selectedRecords,]
    dat4.ts <- dat4s4[!dat4s4$Sample_ID %in% selectedRecords,]

    leafScanTree <- tree(Qclass ~ leafScan_365 +
                         leafScan_385 +
                         leafScan_450 +
                         leafScan_500 +
                         leafScan_587 +
                         leafScan_632 +
                         leafScan_850 +
                         leafScan_880 +
                         leafScan_940,
                     data = dat4.tr)
    dat4.ts$pred.Qclass <- predict(leafScanTree, newdata = dat4.ts, type = "class")

    confMatr <- as.matrix(xtabs(~ Qclass + pred.Qclass, dat4.ts))
    leafScanScores <- data.frame(truPos = confMatr[1,1]/sum(confMatr[1,]),
                                 flsPos = confMatr[1,1]/sum(confMatr[,1]),
                                 ttlCor = sum(diag(confMatr)))
    if(!exists("smryLeafScanScores")){
        smryLeafScanScores <- leafScanScores
    }else{
        smryLeafScanScores <- rbind(smryLeafScanScores, leafScanScores)
    }
}
#### end loop
ndg1 <- Sys.time()
ndg1 - bgn1

head(smryLeafScanScores); tail(smryLeafScanScores)
mean(smryLeafScanScores$ttlCor)/10
table(is.na(smryLeafScanScores$truPos), is.na(smryLeafScanScores$flsPos))
mean(smryLeafScanScores$truPos, na.rm = TRUE) # prob id + | +
mean(smryLeafScanScores$flsPos, na.rm = TRUE) # prob id + | -

## same for supernatant
bgn2 <- Sys.time()
if(exists("smryXtrcScanScores")){
    rm(smryXtrcScanScores)}
#### begin loop
for(i in 1:10){#1000
    ## training (traiN = 84) and test (n = 10) sets
    selectedRecords <- sample(1:nrow(dat4s4), size = traiN)
    selectedRecords <- dat4s4$Sample_ID[sort(selectedRecords)]
    dat4.tr <- dat4s4[dat4s4$Sample_ID %in% selectedRecords,]
    dat4.ts <- dat4s4[!dat4s4$Sample_ID %in% selectedRecords,]

    xtrcScanTree <- tree(Qclass ~ xtrcScan_365 +
                         xtrcScan_385 +
                         xtrcScan_450 +
                         xtrcScan_500 +
                         xtrcScan_587 +
                         xtrcScan_632 +
                         xtrcScan_850 +
                         xtrcScan_880 +
                         xtrcScan_940,
                     data = dat4.tr)
    dat4.ts$pred.Qclass <- predict(xtrcScanTree, newdata = dat4.ts, type = "class")

    confMatr <- as.matrix(xtabs(~ Qclass + pred.Qclass, dat4.ts))
    xtrcScanScores <- data.frame(truPos = confMatr[1,1]/sum(confMatr[1,]),
                                 flsPos = confMatr[1,1]/sum(confMatr[,1]),
                                 ttlCor = sum(diag(confMatr)))
    if(!exists("smryXtrcScanScores")){
        smryXtrcScanScores <- xtrcScanScores
    }else{
        smryXtrcScanScores <- rbind(smryXtrcScanScores, xtrcScanScores)
    }
}
#### end loop
ndg2 <- Sys.time()
ndg2 - bgn2

head(smryXtrcScanScores); tail(smryXtrcScanScores)
mean(smryXtrcScanScores$ttlCor)/10
table(is.na(smryXtrcScanScores$truPos), is.na(smryXtrcScanScores$flsPos))
mean(smryXtrcScanScores$truPos, na.rm = TRUE) # prob id + | +
mean(smryXtrcScanScores$flsPos, na.rm = TRUE) # prob id + | -

## try random forests
## leaves
leafScanRfst <- randomForest(Qclass ~ leafScan_365 +
                                 leafScan_385 +
                                 leafScan_450 +
                                 leafScan_500 +
                                 leafScan_587 +
                                 leafScan_632 +
                                 leafScan_850 +
                                 leafScan_880 +
                                 leafScan_940,
                             data = dat4s4)
print(leafScanRfst)
t(t(round(importance(leafScanRfst), 1)[order(importance(leafScanRfst)[,1], decreasing = TRUE),]))
varImpPlot(leafScanRfst)

## scoring simulations
bgn3 <- Sys.time()
if(exists("smryLeafScanScorRF")){
    rm(smryLeafScanScorRF)}
#### begin loop
for(i in 1:10){#1000
    ## training (traiN = 84) and test (n = 10) sets
    selectedRecords <- sample(1:nrow(dat4s4), size = traiN)
    selectedRecords <- dat4s4$Sample_ID[sort(selectedRecords)]
    dat4.tr <- dat4s4[dat4s4$Sample_ID %in% selectedRecords,]
    dat4.ts <- dat4s4[!dat4s4$Sample_ID %in% selectedRecords,]

    leafScanRfst <- randomForest(Qclass ~ leafScan_365 +
                                     leafScan_385 +
                                     leafScan_450 +
                                     leafScan_500 +
                                     leafScan_587 +
                                     leafScan_632 +
                                     leafScan_850 +
                                     leafScan_880 +
                                     leafScan_940,
                                 data = dat4.tr)
    dat4.ts$pred.Qclass <- predict(leafScanRfst, newdata = dat4.ts, type = "class")

    confMatr <- as.matrix(xtabs(~ Qclass + pred.Qclass, dat4.ts))
    xtrcScanScores <- data.frame(truPos = confMatr[1,1]/sum(confMatr[1,]),
                                 flsPos = confMatr[1,1]/sum(confMatr[,1]),
                                 ttlCor = sum(diag(confMatr)))
    if(!exists("smryLeafScanScorRF")){
        smryLeafScanScorRF <- xtrcScanScores
    }else{
        smryLeafScanScorRF <- rbind(smryLeafScanScorRF, xtrcScanScores)
    }
}
#### end loop
ndg3 <- Sys.time()
ndg3 - bgn3

head(smryLeafScanScorRF); tail(smryLeafScanScorRF)
mean(smryLeafScanScorRF$ttlCor)/10
table(is.na(smryLeafScanScorRF$truPos), is.na(smryLeafScanScorRF$flsPos))
mean(smryLeafScanScorRF$truPos, na.rm = TRUE) # prob id + | +
mean(smryLeafScanScorRF$flsPos, na.rm = TRUE) # prob id + | -

## extracts
xtrcScanRfst <- randomForest(Qclass ~ xtrcScan_365 +
                                 xtrcScan_385 +
                                 xtrcScan_450 +
                                 xtrcScan_500 +
                                 xtrcScan_587 +
                                 xtrcScan_632 +
                                 xtrcScan_850 +
                                 xtrcScan_880 +
                                 xtrcScan_940,
                             data = dat4s4)
print(xtrcScanRfst)
t(t(round(importance(xtrcScanRfst), 1)[order(importance(xtrcScanRfst)[,1], decreasing = TRUE),]))
varImpPlot(xtrcScanRfst)

## scoring simulations
bgn4 <- Sys.time()
if(exists("smryXtrcScanScorRF")){
    rm(smryXtrcScanScorRF)}
#### begin loop
for(i in 1:10){#1000
    ## training (traiN = 84) and test (n = 10) sets
    selectedRecords <- sample(1:nrow(dat4s4), size = traiN)
    selectedRecords <- dat4s4$Sample_ID[sort(selectedRecords)]
    dat4.tr <- dat4s4[dat4s4$Sample_ID %in% selectedRecords,]
    dat4.ts <- dat4s4[!dat4s4$Sample_ID %in% selectedRecords,]

    xtrcScanRfst <- randomForest(Qclass ~ xtrcScan_365 +
                                     xtrcScan_385 +
                                     xtrcScan_450 +
                                     xtrcScan_500 +
                                     xtrcScan_587 +
                                     xtrcScan_632 +
                                     xtrcScan_850 +
                                     xtrcScan_880 +
                                     xtrcScan_940,
                                 data = dat4.tr)
    dat4.ts$pred.Qclass <- predict(xtrcScanRfst, newdata = dat4.ts, type = "class")

    confMatr <- as.matrix(xtabs(~ Qclass + pred.Qclass, dat4.ts))
    xtrcScanScores <- data.frame(truPos = confMatr[1,1]/sum(confMatr[1,]),
                                 flsPos = confMatr[1,1]/sum(confMatr[,1]),
                                 ttlCor = sum(diag(confMatr)))
    if(!exists("smryXtrcScanScorRF")){
        smryXtrcScanScorRF <- xtrcScanScores
    }else{
        smryXtrcScanScorRF <- rbind(smryXtrcScanScorRF, xtrcScanScores)
    }
}
#### end loop
ndg4 <- Sys.time()
ndg4 - bgn4

head(smryXtrcScanScorRF); tail(smryXtrcScanScorRF)
mean(smryXtrcScanScorRF$ttlCor)/10
table(is.na(smryXtrcScanScorRF$truPos), is.na(smryXtrcScanScorRF$flsPos))
mean(smryXtrcScanScorRF$truPos, na.rm = TRUE) # prob id + | +
mean(smryXtrcScanScorRF$flsPos, na.rm = TRUE) # prob id + | -

ndg0 <- Sys.time()
ndg0 - bgn0
