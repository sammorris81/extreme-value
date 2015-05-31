# a is for method 5
# b is for method 4
# c is for methods 1 - 3
#   c: fit.1 is method 1
#   c: fit.2 is method 2
#   c: fit.3 is method 3
setwd("/Volumes/sam-ext/spatial-skew-t/code/analysis/simstudy/results")

for (setting in 7) {
  for (g in 2:10) {
    oldfile <- paste(setting, "-a-", g, ".RData", sep="")
    load(oldfile)
    fit.2 <- fit.1
    
    for (i in 1:5) {
      # filename should be setting-method-dataset.RData
      filename <- paste(setting, "-5-", (g - 1) * 5 + i, ".RData", sep="")
      cat(filename, ": start ")
      fit.1 <- fit.2[[i]]
      save(fit.1, file=filename)
      cat("end \n")
    }
  }
}

for (setting in 1:7) {
  for (g in 1:10) {
    oldfile <- paste(setting, "-b-", g, ".RData", sep="")
    load(oldfile)
    fit.2 <- fit.1
    
    for (i in 1:5) {
      # filename should be setting-method-dataset.RData
      filename <- paste(setting, "-4-", (g - 1) * 5 + i, ".RData", sep="")
      cat(filename, ": start ")
      fit.1 <- fit.2[[i]]
      save(fit.1, file=filename)
      cat("end \n")
    }
  }
}

for (setting in 1:7) {
  for (g in 1:10) {
    oldfile <- paste(setting, "-c-", g, ".RData", sep="")
    load(oldfile)
    fit.4 <- fit.1
    
    for (i in 1:5) {
      # filename should be setting-method-dataset.RData
      filename <- paste(setting, "-3-", (g - 1) * 5 + i, ".RData", sep="")
      cat(filename, ": start ")
      fit.1 <- fit.3[[i]]
      save(fit.1, file=filename)
      cat("end \n")
      
      # filename should be setting-method-dataset.RData
      filename <- paste(setting, "-2-", (g - 1) * 5 + i, ".RData", sep="")
      cat(filename, ": start ")
      fit.1 <- fit.2[[i]]
      save(fit.1, file=filename)
      cat("end \n")
      
      # filename should be setting-method-dataset.RData
      filename <- paste(setting, "-1-", (g - 1) * 5 + i, ".RData", sep="")
      cat(filename, ": start ")
      fit.1 <- fit.4[[i]]
      save(fit.1, file=filename)
      cat("end \n")
    }
  }
}