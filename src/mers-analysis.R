
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

matplot(dates, mdr, type = 's')
matplot(dates, apply(mdr, 2, cumsum), type = 's')

## scraps

ArRiyadCols <- grep("^Ar\\.Riyad[\\.0-9]{0,3}$", names(cmd))
ArRiyadCases <- rowSums(cmd[,ArRiyadCols])
