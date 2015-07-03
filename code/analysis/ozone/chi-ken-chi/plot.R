rm(list=ls())
library(evd)

# We should have 15 Apr - 15 Oct for most years
for (year in 1997:2014) {
  # get the days
  startdate <- as.Date(paste("04/15/", year, sep = ""), "%m/%d/%Y")
  enddate   <- as.Date(paste("10/15/", year, sep = ""), "%m/%d/%Y")
  dates     <- seq(from = startdate, to = enddate, by = 1) 
  
  nrows <- length(dates)
  chicago <- matrix(NA, nrows, 2)
  kenosha <- matrix(NA, nrows, 2)
  chicago[, 1] <- kenosha[, 1] <- dates
  
  chicago.file <- paste(year, "chicago.csv", sep = "")
  kenosha.file <- paste(year, "kenosha.csv", sep = "")
  chicago.this <- read.csv(file = chicago.file)
  kenosha.this <- read.csv(file = kenosha.file)
  
  for (i in 1:length(chicago.this$Daily.Max.8.hour.Ozone.Concentration)) {
    today <- as.Date(chicago.this$Date[i], "%m/%d/%Y")
    this.ozone <- chicago.this$Daily.Max.8.hour.Ozone.Concentration[i]
    if (today %in% dates) {  # only want to include data for dates in range
      idx <- which(dates == today)  # find which row of chicago to fill
      chicago[idx, 2] <- max(c(chicago[idx, 2], this.ozone), na.rm = TRUE)
    }
  }
  
  for (i in 1:length(chicago.this$Daily.Max.8.hour.Ozone.Concentration)) {
    today <- as.Date(chicago.this$Date[i], "%m/%d/%Y")
    this.ozone <- chicago.this$Daily.Max.8.hour.Ozone.Concentration[i]
    if (today %in% dates) {  # only want to include data for dates in range
      idx <- which(dates == today)  # find which row of chicago to fill
      chicago[idx, 2] <- max(c(chicago[idx, 2], this.ozone), na.rm = TRUE)
    }
  }
  
  for (i in 1:length(kenosha.this$Daily.Max.8.hour.Ozone.Concentration)) {
    today <- as.Date(kenosha.this$Date[i], "%m/%d/%Y")
    this.ozone <- kenosha.this$Daily.Max.8.hour.Ozone.Concentration[i]
    if (today %in% dates) {  # only want to include data for dates in range
      idx <- which(dates == today)  # find which row of chicago to fill
      kenosha[idx, 2] <- max(c(kenosha[idx, 2], this.ozone), na.rm = TRUE)
    }
  }
  
  if (year == 1997) {
    chicago.all <- chicago
    kenosha.all <- kenosha
  } else {
    chicago.all <- rbind(chicago.all, chicago)
    kenosha.all <- rbind(kenosha.all, kenosha)
  }
}

chiplot(data = cbind(chicago.all[, 2], kenosha.all[, 2]), 
        which = 1, ylim1 = c(0, 1))
