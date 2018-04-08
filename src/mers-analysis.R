#!/usr/bin/Rscript

cmd <- read.csv("../data/cauchemez-mers-data.csv", skip = 1, colClasses = list(region = "character"))
names(cmd)[1] <- "date"
dates <- lubridate::mdy(cmd$date)
cmd$date <- NULL
row.names(cmd) <- dates


regs <- c("Ar Riyad", "Al Jawf", "Ash Sharqiyah", "Makkah", "Tabuk", "Aseer",
          "Al Madinah", "Najran", "Al Bahah", "Al Hudud ash Shamaliyah", "Al Quassim")

colify <- function(regname, mat){
    regname <- gsub(" ", "\\\\.", regname)
    pattern <- paste0("^", regname, "[\\.0-9]{0,3}$")
    grep(pattern, names(mat))
}

col_inds <- lapply(regs, colify, mat = cmd)

stopifnot(all(diff(sort(unlist(col_inds))) == 1))
stopifnot(max(unlist(col_inds)) == ncol(cmd))

tmpf <- function(x, mat){
    rowSums(mat[, x, drop = FALSE])
}
reg_cases <- lapply(col_inds, tmpf, mat = cmd)

mdr <- do.call(cbind, reg_cases)
colnames(mdr) <- regs

#matplot(dates, mdr, type = 's')
#matplot(dates, apply(mdr, 2, cumsum), type = 's')

acores <- apply(mdr, 2, acf)

hist(colSums(mdr)) ## two clusters
is_big <- which(colSums(mdr) > 50) ## cluster with more cases

matplot(seq_len(nrow(mdr)), mdr[, is_big], type = 'l')

acf(mdr[seq(1, 300), is_big])
acf(mdr[-seq(1, 300), is_big]) ## so with the large hospital outbreaks, ac is higher


## more recent data

don <- read.csv("../data/grepped-don.csv")

library(lubridate)
library(dplyr)

don %>% filter(Reporting.country == "Saudi Arabia") %>%
  mutate(onset_date = ymd(Date.of.symptoms.onset..yyyy.mm.dd.)) %>%
    select(onset_date) -> odate

D <- data.frame(table(odate))
D$odate <- as.Date(D$odate)
D <- padr::pad(D)
D$Freq[is.na(D$Freq)] <- 0
plot((aggregate(ts(D$Freq), 1/7))) ## check that matches WHO epicurve,
                                   ## although that includes imputed
                                   ## onset dates
acf((aggregate(ts(D$Freq), 1/7))) ## ac very low

## rambaut data

rmd <- read.csv("../data/rambaut_cases.csv")

rmd %>% filter(country == "KSA") %>% mutate(odate = ymd(as.character(onset))) %>% select(odate) -> odate1

D1 <- data.frame(table(odate1))
D1$odate1 <- as.Date(D1$odate1)

D1 <- padr::pad(D1)
D1$Freq[is.na(D1$Freq)] <- 0


rts <- aggregate.ts(D1$Freq, 1 / 7)
rts_dates <- seq(ymd("2012-06-06"), by = "7 days", length = length(rts))



## http://www.who.int/emergencies/mers-cov/epi-17-november-2017.png
## epicurve digitization start 2012, week 12

c(0, 0, 0, 0,
  0, 0, 0, 0,
  0, 0, 0, 1,
  0, 0, 0, 0,
  0, 0, 0, 0,
  0, 0, 0, 0,
  0, 0, 0, 0,
  1, 2, 0, 1,
  0, 1, 0, 0,
  0, 0, 0, 0,
  0, 0, 0, 0,
  1, 0, 1, 0,
  0, 2, 0, 0,
  0, 0, 1, 0,
  6, 9, 6, 4,
  4, 4, 7, 5,
  5, 6, 5, 1,
  4, 3, 4, 1,
  3, 7, 6, 8,
  12, 6, 6, 3,
  2, 3, 2, 3,
  1, 2, 2, 2,
  2, 2, 0, 6,
  4, 4, 0, 0,
  2, 1, 3, 1,
  2, 7, 5, 7,
  19, 10, 51, 56,
  90, 92, 55, 26,
  31, 26, 7, 13,
  6, 5, 6, 3,
  0, 0, 0, 0,
  2, 0, 2, 1,
  5, 1, 5, 3,
  5, 6, 14, 12,
  7, 5, 4, 5,
  5, 1, 3, 1,
  4, 4, 3, 6,
  6, 11, 19, 25,
  20, 16, 18, 10,
  5, 7, 1, 3,
  1, 1, 12, 7,
  4, 12, 14, 6,
  2, 6, 3, 2,
  1, 6, 7, 14,
  23, 29, 43, 33,
  32, 7, 6, 4,
  3, 3, 9, 7,
  3, 1, 0, 0,
  2, 0, 1, 3,
  0, 0, 2, 2,
  2, 2, 6, 6,
  10, 7, 19, 14,
  10, 4, 2, 5,
  4, 1, 2, 1,
  1, 0, 1, 2,
  5, 26, 13, 5,
  3, 1, 2, 1,
  3, 2, 1, 1,
  2, 1, 4, 1,
  0, 1, 8, 4,
  4, 5, 5, 3,
  9, 11, 2, 4,
  9, 8, 7, 3,
  5, 4, 6, 4,
  5, 4, 5, 7,
  2, 5, 4, 3,
  2, 3, 5, 4,
  3, 6, 2, 15,
  28, 12, 4, 1,
  1, 0, 3, 5,
  4, 14, 3, 6,
  6, 1, 3, 2,
  2, 2, 1, 3,
  2) -> ec

Ywk <- c(paste0("2012-W", sprintf("%02d", seq(12, 53))),
         paste0("2013-W", sprintf("%02d", seq(1, 52))),
         paste0("2014-W", sprintf("%02d", seq(1, 52))),
         paste0("2015-W", sprintf("%02d", seq(1, 52))),
         paste0("2016-W", sprintf("%02d", seq(1, 52))),
         paste0("2017-W", sprintf("%02d", seq(1, 43))))
Ywk <- paste0(Ywk, "-1")

ecdf <- data.frame(date = ISOweek::ISOweek2date(Ywk), cases = ec)


## Compare with rambaut csv
png("epicurve-vs-rambaut.png")
plot(ecdf, xlab = "Week of Onset", ylab = "Cases", type = 'h',
     main = "Confirmed MERS-CoV cases in Saudi Arabia")
points(rts_dates, rts, col =2)
dev.off()

png("epicurve-sa.png")
plot(ecdf, xlab = "Week of Onset", ylab = "Cases", type = 'h',
     main = "Confirmed MERS-CoV cases in Saudi Arabia")
dev.off()

is_win1 <- ecdf$date < "2015-01-01"
ac1 <- acf(ecdf$cases[is_win1], plot = FALSE)
ac2 <- acf(ecdf$cases[!is_win1], plot = FALSE)

png("autocor-ec.png")
plot(ac1$acf[,1,1], xlab = "Lag", ylab = "Autocorrelation",
     xlim = c(1, 10), pch = 16)
legend("topright", col = c(1, 2), c("before 2015","2015 and later"),
       pch = 16, title = "Window")
points(ac2$acf[,1,1], col = 2, pch = 16)
dev.off()
