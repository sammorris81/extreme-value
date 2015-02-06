# data is the data for the time series
# phi is the parameter for the time series
#   data_t ~ N(phi * data_{t-1}, sqrt(1 - phi^2))
# att.phi is the current number of attempts
# acc.phi is the current number of acceptances
# mh.phi is the candidate standard deviation
# day.mar is the margin in the data that represents time

# MAKE SURE THAT TIME IS LAST MARGIN OF DATA
updatePhiTS <- function(data, phi, day.mar, att, acc, mh) {
  att <- att + 1
  nt <- dim(data)[day.mar]
  if (day.mar == 2) {
    data.up1  <- data[, -1]
    data.lag1 <- data[, -nt]
  } else if (day.mar == 3) {
    data.up1  <- data[, , -1]
    data.lag1 <- data[, , -nt]
  }

  cur.mean <- phi * data.lag1
  cur.sd   <- sqrt(1 - phi^2)

  cur.phi.star <- transform$probit(x=phi, lower=-1, upper=1)
  can.phi.star <- rnorm(1, cur.phi.star, mh)
  can.phi      <- transform$inv.probit(x=can.phi.star, lower=-1, upper=1)
  can.mean     <- can.phi * data.lag1
  can.sd       <- sqrt(1 - can.phi^2)

  # likelihood of data impacted by phi does not include day 1
  R <- sum(dnorm(data.up1, can.mean, can.sd, log=T)) -
       sum(dnorm(data.up1, cur.mean, cur.sd, log=T)) +
       dnorm(can.phi.star, log=T) - dnorm(cur.phi.star, log=T)

  if (!is.na(R)) { if (log(runif(1)) < R) {
    acc <- acc + 1
    phi <- can.phi
  }}

  results <- list(phi=phi, att=att, acc=acc)
}