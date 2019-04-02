rm(list = ls())

library(lattice)

## identify the latest RDS data file available
(availableFiles <- sort(grep("dat4c.RDS",
                             list.files(), value = TRUE)))
dat4c <- readRDS(file = availableFiles[length(availableFiles)])

## model Mg in carrot by Mg available in soil and soil respiration
## identify and select records that are fully populated
fullRecord <- is.finite(dat4c$carrots.xrf.Mg) &
    is.finite(dat4c$soil3.xrf.Mg) &
    is.finite(dat4c$soil9.xrf.Mg) &
    is.finite(dat4c$soil3.display_respiration.data.ugc_gsoil) &
    is.finite(dat4c$soil9.display_respiration.data.ugc_gsoil)
table(fullRecord, useNA = "always")
dat4Mg <- dat4c[fullRecord,]
## center and rescale the predictors

## respiration
par(mfrow = c(2, 2))
qqnorm(dat4Mg$soil9.display_respiration.data.ugc_gsoil)
qqline(dat4Mg$soil9.display_respiration.data.ugc_gsoil)
qqnorm(dat4Mg$soil3.display_respiration.data.ugc_gsoil)
qqline(dat4Mg$soil3.display_respiration.data.ugc_gsoil)
plot(x = dat4Mg$soil9.display_respiration.data.ugc_gsoil,
     y = dat4Mg$soil3.display_respiration.data.ugc_gsoil)
abline(c(0, 1)); abline(h = 0); abline(v = 0)
## look at that high reading
sort(dat4Mg$soil9.display_respiration.data.ugc_gsoil)
sort(dat4Mg$soil3.display_respiration.data.ugc_gsoil)
## mostly respiration in shallow soil > deep
dat4Mg$shDpRatio <- dat4Mg$soil3.display_respiration.data.ugc_gsoil /
    dat4Mg$soil9.display_respiration.data.ugc_gsoil
sort(dat4Mg$shDpRatio)
## the negative reading and the very high reading are troubling
dat4Mg[dat4Mg$shDpRatio <= 0,
       c("soil9.display_respiration.data.ugc_gsoil",
         "soil3.display_respiration.data.ugc_gsoil",
         "shDpRatio")]
dat4Mg[dat4Mg$shDpRatio > 0 & dat4Mg$shDpRatio <= 0.2,
       c("soil9.display_respiration.data.ugc_gsoil",
         "soil3.display_respiration.data.ugc_gsoil",
         "shDpRatio")]
## the astronomical ratio is troubling, too
dat4Mg[dat4Mg$shDpRatio >= 4,
       c("soil9.display_respiration.data.ugc_gsoil",
         "soil3.display_respiration.data.ugc_gsoil",
         "shDpRatio")]
## omitting those would be tantamount to accepting respiration readings b/w 1 and 50
dat4Mg[dat4Mg$soil9.display_respiration.data.ugc_gsoil <= 1 |
       dat4Mg$soil9.display_respiration.data.ugc_gsoil >= 50 |
       dat4Mg$soil3.display_respiration.data.ugc_gsoil <= 1 |
       dat4Mg$soil3.display_respiration.data.ugc_gsoil >= 50,
       c("soil9.display_respiration.data.ugc_gsoil",
         "soil3.display_respiration.data.ugc_gsoil",
         "shDpRatio")]
## but note that all three readings outside that range are from deep soil ...
## move on just don't use deep respiration
dat4Mg$x.soilResp <- c(scale(dat4Mg$soil3.display_respiration.data.ugc_gsoil))
par(mfrow = c(1, 1))
qqnorm(dat4Mg$x.soilResp); qqline(dat4Mg$x.soilResp)

## soil Mg
par(mfrow = c(1, 2))
qqnorm(dat4Mg$soil3.xrf.Mg); qqline(dat4Mg$soil3.xrf.Mg)
qqnorm(dat4Mg$soil9.xrf.Mg); qqline(dat4Mg$soil9.xrf.Mg)
## just use shallow here, too
dat4Mg$x.soilMg <- c(scale(dat4Mg$soil3.xrf.Mg))
par(mfrow = c(1, 1))
qqnorm(dat4Mg$x.soilMg); qqline(dat4Mg$x.soilMg)

par(mfrow = c(1, 2))
plot(x = dat4Mg$x.soilResp,
     y = dat4Mg$carrots.xrf.Mg)
plot(x = dat4Mg$x.soilMg,
     y = dat4Mg$carrots.xrf.Mg)

crrtMg.fm <- lm(carrots.xrf.Mg ~ x.soilResp * x.soilMg,
                data = dat4Mg)
par(mfrow = c(2, 2))
plot(crrtMg.fm) # extreme outlier 513 not influential
dat4Mg[resid(crrtMg.fm) >= 10000,
       c("soil3.display_respiration.data.ugc_gsoil", "x.soilResp",
         "soil3.xrf.Mg", "x.soilMg",
         "carrots.xrf.Mg")]
sort(dat4Mg$carrots.xrf.Mg)
## but it is quite extreme
dat4Mg$crrtMg.outlier <- dat4Mg$carrots.xrf.Mg >= 10000
par(mfrow = c(1, 2))
qqnorm(dat4Mg$carrots.xrf.Mg); qqline(dat4Mg$carrots.xrf.Mg)
qqnorm(dat4Mg$carrots.xrf.Mg[!dat4Mg$crrtMg.outlier])
qqline(dat4Mg$carrots.xrf.Mg[!dat4Mg$crrtMg.outlier])
## model fit without it
crrtMg.fm2 <- lm(carrots.xrf.Mg ~ x.soilResp * x.soilMg,
                data = dat4Mg[!dat4Mg$crrtMg.outlier,])
par(mfrow = c(2, 2))
plot(crrtMg.fm2)
summary(crrtMg.fm)
summary(crrtMg.fm2)

levelplot(carrots.xrf.Mg ~ x.soilResp * x.soilMg,
          data = dat4Mg[!dat4Mg$crrtMg.outlier,])

## plot the fitted model's predictions
dat5Mg <- expand.grid(x.soilResp = seq(min(dat4Mg$x.soilResp),
                                       max(dat4Mg$x.soilResp),
                                       length.out = 1001),
                      x.soilMg = seq(min(dat4Mg$x.soilMg),
                                     max(dat4Mg$x.soilMg),
                                     length.out = 1001))
dat5Mg$pred.crrtMg <- predict(crrtMg.fm2, newdata = dat5Mg)
contourplot(pred.crrtMg ~  x.soilResp * x.soilMg, data = dat5Mg)

## try Si, P (regulated, up); Si first ...
## identify and select records that are fully populated
fullRecor2 <- is.finite(dat4c$carrots.xrf.Si) &
    is.finite(dat4c$soil3.xrf.Si) &
    is.finite(dat4c$soil9.xrf.Si) &
    is.finite(dat4c$soil3.display_respiration.data.ugc_gsoil) &
    is.finite(dat4c$soil9.display_respiration.data.ugc_gsoil)
table(fullRecord, fullRecor2, useNA = "always") ## same ones
dat4Si <- dat4c[fullRecor2,]

## center and rescale the predictors
## respiration
dat4Si$x.soilResp <- c(scale(dat4Si$soil3.display_respiration.data.ugc_gsoil))

## soil Si
qqnorm(dat4Si$soil3.xrf.Si); qqline(dat4Si$soil3.xrf.Si)
dat4Si$x.soilSi <- c(scale(dat4Si$soil3.xrf.Si))

crrtSi.fm <- lm(carrots.xrf.Si ~ x.soilResp * x.soilSi,
                data = dat4Si)
par(mfrow = c(2, 2))
plot(crrtSi.fm) # odd but perhaps indicative that 513 is the outlier again
sort(dat4Si$carrots.xrf.Si)
dat4Si[dimnames(dat4Si)[[1]] == "513",
       c("carrots.xrf.Mg", "carrots.xrf.Si")]
dat4Si$crrtSi.outlier <- dat4Si$carrots.xrf.Si >= 700
crrtSi.fm2 <- lm(carrots.xrf.Si ~ x.soilResp * x.soilSi,
                data = dat4Si[!dat4Si$crrtSi.outlier,])
par(mfrow = c(2, 2))
plot(crrtSi.fm2)
summary(crrtSi.fm)
summary(crrtSi.fm2)

levelplot(carrots.xrf.Si ~ x.soilResp * x.soilSi,
          data = dat4Si[!dat4Si$crrtSi.outlier,])

## plot the fitted model's predictions
dat5Si <- expand.grid(x.soilResp = seq(min(dat4Si$x.soilResp),
                                       max(dat4Si$x.soilResp),
                                       length.out = 1001),
                      x.soilSi = seq(min(dat4Si$x.soilSi),
                                     max(dat4Si$x.soilSi),
                                     length.out = 1001))
dat5Si$pred.crrtSi <- predict(crrtSi.fm2, newdata = dat5Si)
contourplot(pred.crrtSi ~  x.soilResp * x.soilSi, data = dat5Si)

## ... try P (regulated, up)
## identify and select records that are fully populated
fullRecor2 <- is.finite(dat4c$carrots.xrf.P) &
    is.finite(dat4c$soil3.xrf.P) &
    is.finite(dat4c$soil9.xrf.P) &
    is.finite(dat4c$soil3.display_respiration.data.ugc_gsoil) &
    is.finite(dat4c$soil9.display_respiration.data.ugc_gsoil)
table(fullRecord, fullRecor2, useNA = "always") ## same ones
dat4P <- dat4c[fullRecor2,]

## center and rescale the predictors
## respiration
dat4P$x.soilResp <- c(scale(dat4P$soil3.display_respiration.data.ugc_gsoil))

## soil P
qqnorm(dat4P$soil3.xrf.P); qqline(dat4P$soil3.xrf.P)
dat4P$x.soilP <- c(scale(dat4P$soil3.xrf.P))

crrtP.fm <- lm(carrots.xrf.P ~ x.soilResp * x.soilP,
                data = dat4P)
par(mfrow = c(2, 2))
plot(crrtP.fm) # 513 again
sort(dat4P$carrots.xrf.P)
dat4P[dimnames(dat4P)[[1]] == "513",
       c("carrots.xrf.Mg", "carrots.xrf.Si", "carrots.xrf.P")]
dat4P$crrtP.outlier <- dat4P$carrots.xrf.P >= 11300
crrtP.fm2 <- lm(carrots.xrf.P ~ x.soilResp * x.soilP,
                data = dat4P[!dat4P$crrtP.outlier,])
par(mfrow = c(2, 2))
plot(crrtP.fm2)
summary(crrtP.fm)
summary(crrtP.fm2)

levelplot(carrots.xrf.P ~ x.soilResp * x.soilP,
          data = dat4P[!dat4P$crrtP.outlier,])

## plot the fitted model's predictions
dat5P <- expand.grid(x.soilResp = seq(min(dat4P$x.soilResp),
                                       max(dat4P$x.soilResp),
                                       length.out = 1001),
                      x.soilP = seq(min(dat4P$x.soilP),
                                     max(dat4P$x.soilP),
                                     length.out = 1001))
dat5P$pred.crrtP <- predict(crrtP.fm2, newdata = dat5P)

contourplot(pred.crrtP ~  x.soilResp * x.soilP, data = dat5P)

## try Ca, Se (regulated, ambient); Ca first ...
## identify and select records that are fully populated
fullRecor2 <- is.finite(dat4c$carrots.xrf.Ca) &
    is.finite(dat4c$soil3.xrf.Ca) &
    is.finite(dat4c$soil9.xrf.Ca) &
    is.finite(dat4c$soil3.display_respiration.data.ugc_gsoil) &
    is.finite(dat4c$soil9.display_respiration.data.ugc_gsoil)
table(fullRecord, fullRecor2, useNA = "always") ## same ones
dat4Ca <- dat4c[fullRecor2,]

## center and rescale the predictors
## respiration
dat4Ca$x.soilResp <- c(scale(dat4Ca$soil3.display_respiration.data.ugc_gsoil))

## soil Ca
qqnorm(dat4Ca$soil3.xrf.Ca); qqline(dat4Ca$soil3.xrf.Ca)
dat4Ca$x.soilCa <- c(scale(dat4Ca$soil3.xrf.Ca))

crrtCa.fm <- lm(carrots.xrf.Ca ~ x.soilResp * x.soilCa,
                data = dat4Ca)
par(mfrow = c(2, 2))
plot(crrtCa.fm) # 513 right up there with 493
sort(dat4Ca$carrots.xrf.Ca)
dat4Ca[dimnames(dat4Ca)[[1]] %in% c("513", "493"),
       c("carrots.xrf.Mg", "carrots.xrf.Si", "carrots.xrf.P", "carrots.xrf.Ca")]
## close as an outlier, actually larger as a datum
dat4Ca$crrtCa.outlier <- dat4Ca$carrots.xrf.Ca >= 10511
crrtCa.fm2 <- lm(carrots.xrf.Ca ~ x.soilResp * x.soilCa,
                data = dat4Ca[!dat4Ca$crrtCa.outlier,])
par(mfrow = c(2, 2))
plot(crrtCa.fm2)
summary(crrtCa.fm)
summary(crrtCa.fm2)

levelplot(carrots.xrf.Ca ~ x.soilResp * x.soilCa,
          data = dat4Ca[!dat4Ca$crrtCa.outlier,])

## plot the fitted model's predictions
dat5Ca <- expand.grid(x.soilResp = seq(min(dat4Ca$x.soilResp),
                                       max(dat4Ca$x.soilResp),
                                       length.out = 1001),
                      x.soilCa = seq(min(dat4Ca$x.soilCa),
                                     max(dat4Ca$x.soilCa),
                                     length.out = 1001))
dat5Ca$pred.crrtCa <- predict(crrtCa.fm2, newdata = dat5Ca)
contourplot(pred.crrtCa ~  x.soilResp * x.soilCa, data = dat5Ca)

## ... try Se (regulated, ambient)
## identify and select records that are fully populated
fullRecor2 <- is.finite(dat4c$carrots.xrf.Se) &
    is.finite(dat4c$soil3.xrf.Se) &
    is.finite(dat4c$soil9.xrf.Se) &
    is.finite(dat4c$soil3.display_respiration.data.ugc_gsoil) &
    is.finite(dat4c$soil9.display_respiration.data.ugc_gsoil)
table(fullRecord, fullRecor2, useNA = "always") ## 6 with something missing ...
## ... make sure it isn't soil9.xrf.Se
dat4c[fullRecord & !fullRecor2,
      c("carrots.xrf.Se", "soil3.xrf.Se", "soil9.xrf.Se",
        "soil3.display_respiration.data.ugc_gsoil",
        "soil9.display_respiration.data.ugc_gsoil")]
## right - missing from carrot
dat4Se <- dat4c[fullRecor2,]

## center and rescale the predictors
## respiration
dat4Se$x.soilResp <- c(scale(dat4Se$soil3.display_respiration.data.ugc_gsoil))

## soil Se
qqnorm(dat4Se$soil3.xrf.Se); qqline(dat4Se$soil3.xrf.Se)
sort(dat4Se$soil3.xrf.Se)
dat4Se <- dat4Se[dat4Se$soil3.xrf.Se < 5,]
qqnorm(dat4Se$soil3.xrf.Se); qqline(dat4Se$soil3.xrf.Se)
dat4Se$x.soilSe <- c(scale(dat4Se$soil3.xrf.Se))

crrtSe.fm <- lm(carrots.xrf.Se ~ x.soilResp * x.soilSe,
                data = dat4Se)
par(mfrow = c(2, 2))
plot(crrtSe.fm) # hurray! 513 nowhere to be seen
sort(dat4Se$carrots.xrf.Se)
dat4Se[dimnames(dat4Se)[[1]] == "513",
       c("carrots.xrf.Mg", "carrots.xrf.Si", "carrots.xrf.P", "carrots.xrf.Ca",
         "carrots.xrf.Se")]
summary(crrtSe.fm)

levelplot(carrots.xrf.Se ~ x.soilResp * x.soilSe,
          data = dat4Se)

## plot the fitted model's predictions
dat5Se <- expand.grid(x.soilResp = seq(min(dat4Se$x.soilResp),
                                       max(dat4Se$x.soilResp),
                                       length.out = 1001),
                      x.soilSe = seq(min(dat4Se$x.soilSe),
                                     max(dat4Se$x.soilSe),
                                     length.out = 1001))
dat5Se$pred.crrtSe <- predict(crrtSe.fm, newdata = dat5Se)
contourplot(pred.crrtSe ~  x.soilResp * x.soilSe, data = dat5Se)
