rm(list = ls())

library(MASS)

## identify the latest RDS data file available
(availableFiles <- sort(grep("soil-food-combined-clean.RDS",
                             list.files(), value = TRUE)))
dat4 <- readRDS(file = availableFiles[length(availableFiles)])

## I have noticed that sample_type doesn't exactly correspond to populated
## carrot/spinach fields
## also there are records that are carrot+soil or spinach+soil so dummy variables

dat4$carrot <- as.numeric(
    apply(dat4[,grep("carrot", names(dat4), ignore.case = TRUE)],
          1,
          function(f) !all(is.na(f) | f == "")))
xtabs(~ carrot, dat4)

dat4$spinach <- as.numeric(
    apply(dat4[,grep("spinach", names(dat4), ignore.case = TRUE)],
          1,
          function(f) !all(is.na(f) | f == "")))
xtabs(~ spinach, dat4)

dat4$soil <- as.numeric(
    apply(dat4[,grep("soil", names(dat4), ignore.case = TRUE)],
          1,
          function(f) !all(is.na(f) | f == "")))
xtabs(~ soil, dat4)
data.frame(xtabs(~ soil + spinach + carrot, dat4))[8:1,c(3:1,4)]
## 422 just-carrot
## 151 carrot-and-soil
## 110 just-spinach
##  37 spinach-and-soil
## also:
##   6 carrot _and_ spinach ... really?!
##   3 just-soil
##   3 with no data about carrots, spinach, or soil

## candidate predictors as dummy variables
(foo <-
     data.frame(sort(table(matrix(unlist(strsplit(dat4$farm_practice, " ")))[,1]),
                     decreasing = TRUE))
)
for(i in as.character(foo$Var1)){
    dat4[i] <- as.numeric(grepl(i, dat4$farm_practice))
}
## review a random sample
dat4[sample((1:nrow(dat4)), size = 10),
      c("farm_practice", names(dat4)[length(names(dat4)) - ((nrow(foo)-1):0)])]
## look at the "nones"
dat4[dat4$none == 1,
       c("farm_practice", names(dat4)[length(names(dat4)) - ((nrow(foo)-1):0)])]
## (but that isn't really a predictor - 0 under all others)
## and the other seldom-used description
dat4[dat4$biodynamic == 1,
      c("farm_practice", names(dat4)[length(names(dat4)) - ((nrow(foo)-1):0)])]

## UPSHOT:
## focus on samples with vegetable-and-soil
nrow(dat4)
nrow(dat4c <- dat4[dat4$carrot == 1,])
nrow(dat4s <- dat4[dat4$spinach == 1,])

## drop uninformative columns ...
(nonC <- c(names(dat4c)[sapply(dat4c, function(f){length(unique(f)) == 1})],
           grep("spinach", names(dat4c), ignore.case = TRUE, value = TRUE)))
dat4c <- dat4c[,setdiff(names(dat4), nonC)]

(nonS <- c(names(dat4s)[sapply(dat4s, function(f){length(unique(f)) == 1})],
           grep("carrot", names(dat4s), ignore.case = TRUE, value = TRUE)))
dat4s <- dat4s[,setdiff(names(dat4), nonS)]

saveRDS(dat4c, file = paste(gsub("-", "", Sys.Date()),
                            "dat4c.RDS",
                            sep = ""))

saveRDS(dat4s, file = paste(gsub("-", "", Sys.Date()),
                            "dat4s.RDS",
                            sep = ""))
