rm(list = ls())

library(readxl)

## update-able filename
(fileName <- grep(".xls", list.files(recursive = TRUE), value = TRUE))
dim(dat4 <- as.data.frame(read_xlsx(fileName)))
## add line number to permit re-ordering and reference back to Excel file
dat4$xlLine <- 2:(nrow(dat4)+1)

## there are apparent duplicate Sample_IDs
IDfreq <- data.frame(xtabs(~ Sample_ID, data = dat4))
xtabs(~ Freq, IDfreq)

## review the apparent duplicates
IDfreq[IDfreq$Freq == 2,]

foo <- dat4[dat4$Sample_ID == "630",]#"189"#"619"
## verify: if they're both populated, they're identical
table(sapply(foo, function(f){
    (!is.na(f[1]) & !is.na(f[2]) & identical(f[1], f[2])) |
        is.na(f[1]) |
        is.na(f[2])
}))
## both populated, not identical - only the field I added
table(sapply(foo, function(f){!is.na(f[1]) & !is.na(f[2]) & f[2] != f[1]}))
names(foo)[sapply(foo, function(f){!is.na(f[1]) & !is.na(f[2]) & f[2] != f[1]})]
## neither NA, identical
table(sapply(foo, function(f){!is.na(f[1]) & !is.na(f[2]) & f[2] == f[1]}))
## both NA
table(sapply(foo, function(f){all(is.na(f))}))
## earlier record is populated and later one is NA
table(sapply(foo, function(f){!is.na(f[1]) & is.na(f[2])}))
## vice versa
table(sapply(foo, function(f){is.na(f[1]) & !is.na(f[2])}))
foo[,sapply(foo, function(f){is.na(f[1]) & !is.na(f[2])})]

## for each pair, keep earlier record
## where earlier NA and later populated, move data up
for(i in IDfreq$Sample_ID[IDfreq$Freq == 2]){
    focalRows <- which(dat4$Sample_ID == i)
    for(j in 1:(ncol(dat4) - 1)){# last column is my xlLine
        if(is.na(dat4[focalRows[1], j]) & !is.na(dat4[focalRows[2], j])){
            dat4[focalRows[1], j] <- dat4[focalRows[2], j]}
    }
    dat4 <- dat4[-focalRows[2],]
}

saveRDS(dat4, file = paste(gsub("-", "", Sys.Date()),
                           "soil-food-combined-clean.RDS",
                           sep = ""))
