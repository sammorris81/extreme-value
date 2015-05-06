#############################################################:
###       THE MAIN MCMC FUNCTION TO FIT THE PS MODEL      ###:
#############################################################:

#############################################################:
#
# INPUTS:
#
#   y      = nyears x nsites matrix of data
#   x      = nyears x nsites matrix of covariates
#   s      = nsites x 2 matrix of spatial locations
#   knots  = nknots x 2 matrix of knots
#   thresh = the threshold
#   xp     = covariates at the prediction sites
#   sp     = locations of the prediction sites
#   mnB    = prior means (5 x 3 matrix)
#   sdB    = prior sds (5 x 3 matrix)
#   iters  = number of MCMC samples to generate
#   burn   = number of discard as burn-in
#   update = samples between graphical displays
#
# OUTPUTS:
#
#   samples = Posterior samples of the model parameters. The 2 columns correspond to the
#             intercept and covariate. The five rows correspond to
#
#             1) probability below the threshold
#             2) GPD scale
#             3) GPD shape
#             4) alpha
#             5) kernel bandwidth
#
###################################################################################
if (!exists("dPS_cpp_sca") | !exists("dPS_cpp_mat")) {
  library(Rcpp)
  Sys.setenv("PKG_CXXFLAGS"="-fopenmp")
  Sys.setenv("PKG_LIBS"="-fopenmp")
  sourceCpp(file = "./llps.cpp")
}

maxstable<-function(y, x, s, thresh, knots, sp = NULL, xp = NULL,
                    mnB = NULL, sdB = NULL, initB = NULL,
                    iters = 50000, burn = 10000, thin = 5,
                    threads = 1, update = 100, iterplot = FALSE){


    ny  <- dim(y)[1]
    ns  <- dim(y)[2]
    nk  <- dim(knots)[1]
    dw2 <- as.matrix(rdist(s, knots))^2

    #INITIAL VALUES:
    nb <- 3  # redoing to fit int, lat, long
    if (is.null(initB)) {
      beta <- matrix(0, nb, 5)
      beta[1, 1] <- 2
      beta[1, 2] <- 1
      beta[1, 3] <- 0.01
      beta[1, 4] <- 0
      beta[1, 5] <- 0
    }

    if (is.null(mnB)) {
      mnB = matrix(0, nb, 5)
    } else if (length(mnB) == 1) {
      mnB = matrix(mnB, nb, 5)
    } else if (nrow(mnB) != nb) {
      stop("mnB must have 3 rows")
    }

    if (is.null(sdB)) {
      sdB = matrix(1, nb, 5)
    } else if (length(sdB) == 1) {
      sdB = matrix(sdB, nb, 5)
    } else if (nrow(sdB) != nb) {
      stop("sdB must have 3 rows")
    }

    B <- array(0, c(ny, ns, 5))
    for (j in 1:5) {
      B[, , j] <- make.B(s, beta[, j], j)
    }

    npts <- 50
    Ubeta <- qbeta(seq(0, 1, length=npts + 1), 0.5, 0.5)
    MidPoints <- (Ubeta[-1] + Ubeta[-(npts + 1)]) / 2
    psi <- MidPoints * pi
    BinWidth <- Ubeta[-1] - Ubeta[-(npts + 1)]

    FAC <- fac2FAC(make.fac(dw2, B[1, 1, 5]))
    a   <- matrix(1, ny, nk)
    A   <- matrix(1, ny, ns)
    for (t in 1:ny) {
      A[t, ]    <- a2A(FAC, a[t, ], B[1, 1, 4])
    }
    curll <- matrix(0, ny, ns)
    curlp <- dPS_cpp_mat(a, B[1, 1, 4], psi, BinWidth, threads)
    for (t in 1:ny) {
      curll[t, ] <- loglike(y[t, ], A[t, ], B[t, , ], thresh)
#       for (k in 1:nk) {
#         curlp[t, k] <- dPS_cpp_sca(a[t, k], B[1, 1, 4], MidPoints, BinWidth, threads)
#       }
    }

    samples <- array(0, c(iters, nrow(beta), ncol(beta)))
    dimnames(samples)[[1]] <- paste("Iteration", 1:iters)
    dimnames(samples)[[2]] <- c("Intercept", "Latitude", "Longitude")
    dimnames(samples)[[3]] <- c("Prob extreme", "GPD scale", "GPD shape", "Alpha", "Bandwidth")


    cuts <- exp(c(-1, 0, 1, 2, 5, 10))
    MHa  <- rep(1, 100)
    atta <- acca <- 0 * MHa
    attb <- accb <- MHb <- matrix(c(0.2, 0.1, 0.1, 0.02, 0.02), nb, 5, byrow=T)

    yp <- NULL
    np <- 0
    if (!is.null(sp)) {
      np   <- nrow(sp)
      yp   <- array(0, c(iters, ny, np))
      dwp2 <- as.matrix(rdist(sp, knots))^2
    }

    for (i in 1:iters) {

     for (ttt in 1:thin) {
      ####################################################
      ##############      Random effects a    ############
      ####################################################

      olda <- a
      FAC.star <- FAC^(1 / B[1, 1, 4])
      for (t in 1:ny) {
        parts <- loglikeparts(y[t, ], A[t, ], B[t, , ], thresh)
        W     <- ifelse(parts$above, 1, 0)
        ccc   <- -A[t, ] * parts$expo + W * log(A[t, ])
        alpha <- B[1, 1, 4]

        for (k in 1:nk) {
         # l1     <- get.level(a[t, k], cuts)
         l1     <- get.level.sca(a[t, k], cuts)
         cana   <- exp(rnorm(1, log(a[t, k]), MHa[l1]))
         # l2     <- get.level(cana, cuts)
         l2     <- get.level.sca(a[t, k], cuts)
         WWW    <- FAC.star[, k]
         canA   <- A[t, ] + WWW * (cana - a[t, k])
         cc     <- -canA * parts$expo + W * log(canA)
         # canlp  <- dPS(cana, alpha)
         canlp  <- dPS_cpp_sca(cana, alpha, psi, BinWidth, threads)
         loga <- log(a[t, k])
         logcana <- log(cana)
         R <- sum(cc - ccc) + canlp - curlp[t, k] +
              dlognormal(loga, logcana, MHa[l2]) -
              dlognormal(logcana, loga, MHa[l1])
         if (!is.na(exp(R))) { if (runif(1) < exp(R)) {
             a[t, k]     <- cana
             A[t, ]      <- canA
             ccc         <- cc
             curlp[t, k] <- canlp
         } }
       }

       curll[t, ] <- loglike(y[t, ], A[t, ], B[t, , ], thresh)
      }

      ####################################################
      ##############    Model parameters, B   ############
      ####################################################

      for (j in 1:nb) { for (k in 1:5) { if ((k == 1) | (k == 2) | (j == 1)) {
       attb[j, k] <- attb[j, k] + 1
       canbeta <- beta
       canB    <- B
       canA    <- A
       canll   <- curll
       canlp   <- curlp
       canFAC  <- FAC

       canbeta[j, k] <- beta[j, k] + MHb[j, k] * rnorm(1)
       canB[, , k]   <- make.B(s, canbeta[, k], k)
       if (k > 3) {
         alpha  <- canB[1, 1, 4]
         gamma  <- canB[1, 1, 5]
         canFAC <- fac2FAC(make.fac(dw2, gamma))
         for (t in 1:ny) {
           canA[t, ] <- a2A(canFAC, a[t, ], alpha)
         }
       }

       for (t in 1:ny) {
         canll[t, ] <- loglike(y[t, ], canA[t, ], canB[t, , ], thresh)
       }

       if (k == 4) {
         canlp <- dPS_cpp_mat(a, canB[1, 1, 4], psi, BinWidth, threads)
       }
#        if (k == 4) { for (t in 1:ny) { for (l in 1:nk) {
#          canlp[t, l] <- dPS_cpp_sca(a[t, l], canB[1, 1, 4], MidPoints,
#                                     BinWidth, threads)
#        } } }

       R <- sum(canll - curll) +
            sum(canlp - curlp) +
            dnorm(canbeta[j, k], mnB[j, k], sdB[j, k], log=T) -
            dnorm(beta[j, k], mnB[j, k], sdB[j, k], log=T)

       if (!is.na(exp(R))) { if (log(runif(1)) < R) {
          accb[j, k] <- accb[j, k] + 1
          beta       <- canbeta
          B          <- canB
          A          <- canA
          curll      <- canll
          curlp      <- canlp
          FAC        <- canFAC
       } }
      } } }
     }#end thin

     samples[i, , ] <- beta

     #Predictions
     if (np > 0) {
      Bp <- array(0, c(ny, np, 5))
      for (j in 1:5) {
        Bp[, , j] <- make.B(sp, beta[, j], j)
      }

      for (t in 1:ny) {
        Ap         <- a2AP(dwp2, a[t, ], B[1, 1, 4], B[1, 1, 5])
        U          <- rGEV(np, 1, B[1, 1, 4], B[1, 1, 4])
        X          <- Ap * U
        tau        <- exp(-1 / X)
        yp[i, t, ] <- qGPD(tau, Bp[t, , ], thresh)
      }
     }


     level <- get.level(olda, cuts)
     for (j in 1:length(MHa)) {
         acca[j] <- acca[j] + sum(olda[level == j] != a[level == j])
         atta[j] <- atta[j] + sum(level == j)
         if ((i < burn / 2) & (atta[j] > 100)) {
           if (acca[j] / atta[j] < 0.3) { MHa[j] <- MHa[j] * 0.9 }
           if (acca[j] / atta[j] > 0.6) { MHa[j] <- MHa[j] * 1.1 }
           acca[j] <- atta[j] <- 0
         }
     }


      for (j in 1:nb) { for (k in 1:5) { if ((i < burn / 2) & (attb[j, k] > 50)) {
        if (accb[j, k] / attb[j, k] < 0.3) { MHb[j, k] <- MHb[j, k] * 0.9 }
        if (accb[j, k] / attb[j, k] > 0.6) { MHb[j, k] <- MHb[j, k] * 1.1 }
         accb[j, k] <- attb[j, k] <- 0
      } } }


     #DISPLAY CURRENT VALUE:
      if (iterplot) { if((i %% update) == 0) {
        par(mfrow=c(5, 2), mar=c(2, 2, 2, 2))
        plot(samples[1:i, 1, 1], main="Prob0 Int", type="l")
        abline(h=0)
        plot(samples[1:i, 2, 1], main="Prob0 Slope", type="l")
        abline(h=0)
        plot(samples[1:i, 1, 2], main="Scale Int", type="l")
        abline(h=0)
        plot(samples[1:i, 2, 2], main="Scale Slope", type="l")
        abline(h=0)
        plot(samples[1:i, 1, 3], main="Shape Int", type="l")
        abline(h=0)
        plot(samples[1:i, 2, 3], main="Shape Slope", type="l")
        abline(h=0)
        plot(samples[1:i, 1, 4], main="Alpha Int", type="l")
        abline(h=0)
        plot(samples[1:i, 2, 4], main="Alpha Slope", type="l")
        abline(h=0)
        plot(samples[1:i, 1, 5], main="BW Int", type="l")
        abline(h=0)
        plot(samples[1:i, 2, 5], main="BW Slope", type="l")
        abline(h=0)
      } }

      if ((i %% update) == 0) {
        cat("\t iter", i, "\n")
      }
    }

  results <- list(samples=samples, yp=yp)
  return(results)
}

# maxstableold<-function(y, x, s, thresh, knots, sp = NULL, xp = NULL,
#                     mnB = NULL, sdB = NULL, initB = NULL,
#                     iters = 50000, burn = 10000, thin = 5,
#                     threads = 1, update = 100, iterplot = FALSE){
#
#
#   ny  <- dim(y)[1]
#   ns  <- dim(y)[2]
#   nk  <- dim(knots)[1]
#   dw2 <- as.matrix(rdist(s, knots))^2
#
#   #INITIAL VALUES:
#   nb <- 3  # redoing to fit int, lat, long
#   if (is.null(initB)) {
#     beta <- matrix(0, nb, 5)
#     beta[1, 1] <- 2
#     beta[1, 2] <- 1
#     beta[1, 3] <- 0.01
#     beta[1, 4] <- 0
#     beta[1, 5] <- 0
#   }
#
#   if (is.null(mnB)) {
#     mnB = matrix(0, nb, 5)
#   } else if (length(mnB) == 1) {
#     mnB = matrix(mnB, nb, 5)
#   } else if (nrow(mnB) != nb) {
#     stop("mnB must have 3 rows")
#   }
#
#   if (is.null(sdB)) {
#     sdB = matrix(1, nb, 5)
#   } else if (length(sdB) == 1) {
#     sdB = matrix(sdB, nb, 5)
#   } else if (nrow(sdB) != nb) {
#     stop("sdB must have 3 rows")
#   }
#
#   B <- array(0, c(ny, ns, 5))
#   for (j in 1:5) {
#     B[, , j] <- make.B(s, beta[, j], j)
#   }
#
#   npts <- 100
#   Ubeta <- qbeta(seq(0, 1, length=npts + 1), 0.5, 0.5)
#   MidPoints <- (Ubeta[-1] + Ubeta[-(npts + 1)]) / 2
#   BinWidth <- Ubeta[-1] - Ubeta[-(npts + 1)]
#
#   FAC <- fac2FAC(make.fac(dw2, B[1, 1, 5]))
#   a   <- matrix(1, ny, nk)
#   A   <- matrix(1, ny, ns)
#   for (t in 1:ny) {
#     A[t, ]    <- a2A(FAC, a[t, ], B[1, 1, 4])
#   }
#   curll <- matrix(0, ny, ns)
#   curlp <- matrix(0, ny, nk)
#   for (t in 1:ny) {
#     curll[t, ] <- loglike(y[t, ], A[t, ], B[t, , ], thresh)
#     for (k in 1:nk) {
#       curlp[t, k] <- dPS(a[t, k], B[1, 1, 4])
#     }
#   }
#
#   samples <- array(0, c(iters, nrow(beta), ncol(beta)))
#   dimnames(samples)[[1]] <- paste("Iteration", 1:iters)
#   dimnames(samples)[[2]] <- c("Intercept", "Latitude", "Longitude")
#   dimnames(samples)[[3]] <- c("Prob extreme", "GPD scale", "GPD shape", "Alpha", "Bandwidth")
#
#
#   cuts <- exp(c(-1, 0, 1, 2, 5, 10))
#   MHa  <- rep(1, 100)
#   atta <- acca <- 0 * MHa
#   attb <- accb <- MHb <- matrix(c(0.2, 0.1, 0.1, 0.02, 0.02), nb, 5, byrow=T)
#
#   yp <- NULL
#   np <- 0
#   if (!is.null(sp)) {
#     np   <- nrow(sp)
#     yp   <- array(0, c(iters, ny, np))
#     dwp2 <- as.matrix(rdist(sp, knots))^2
#   }
#
#   for (i in 1:iters) {
#
#     for (ttt in 1:thin) {
#       ####################################################
#       ##############      Random effects a    ############
#       ####################################################
#
#       olda <- a
#       for (t in 1:ny) {
#         parts <- loglikeparts(y[t, ], A[t, ], B[t, , ], thresh)
#         W     <- ifelse(parts$above, 1, 0)
#         ccc   <- -A[t, ] * parts$expo + W * log(A[t, ])
#         alpha <- B[1, 1, 4]
#
#
#         for (k in 1:nk) {
#           l1     <- get.level(a[t, k], cuts)
#           # l1     <- get.level.sca(a[t, k], cuts)
#           cana   <- exp(rnorm(1, log(a[t, k]), MHa[l1]))
#           l2     <- get.level(cana, cuts)
#           # l2     <- get.level.sca(a[t, k], cuts)
#           WWW    <- FAC[, k]^(1 / alpha)
#           canA   <- A[t, ] + WWW * (cana - a[t, k])
#           cc     <- -canA * parts$expo + W * log(canA)
#           canlp  <- dPS(cana, alpha)
#           # canlp  <- dPS_cpp_sca(cana, alpha, MidPoints, BinWidth, threads)
#           R <- sum(cc - ccc) +
#             canlp - curlp[t, k] +
#             dlognormal(a[t, k], cana, MHa[l2]) -
#             dlognormal(cana, a[t, k], MHa[l1])
#           if (!is.na(exp(R))) { if (runif(1) < exp(R)) {
#             a[t, k]     <- cana
#             A[t, ]      <- canA
#             ccc         <- cc
#             curlp[t, k] <- canlp
#           } }
#         }
#
#         curll[t, ] <- loglike(y[t, ], A[t, ], B[t, , ], thresh)
#       }
#
#       ####################################################
#       ##############    Model parameters, B   ############
#       ####################################################
#
#       for (j in 1:nb) { for (k in 1:5) { if ((k == 1) | (k == 2) | (j == 1)) {
#         attb[j, k] <- attb[j, k] + 1
#         canbeta <- beta
#         canB    <- B
#         canA    <- A
#         canll   <- curll
#         canlp   <- curlp
#         canFAC  <- FAC
#
#         canbeta[j, k] <- beta[j, k] + MHb[j, k] * rnorm(1)
#         canB[, , k]   <- make.B(s, canbeta[, k], k)
#         if (k > 3) {
#           alpha  <- canB[1, 1, 4]
#           gamma  <- canB[1, 1, 5]
#           canFAC <- fac2FAC(make.fac(dw2, gamma))
#           for (t in 1:ny) {
#             canA[t, ] <- a2A(canFAC, a[t, ], alpha)
#           }
#         }
#
#         for (t in 1:ny) {
#           canll[t, ] <- loglike(y[t, ], canA[t, ], canB[t, , ], thresh)
#         }
#
# #         if (k == 4) {
# #           canlp <- dPS_cpp_mat(a, canB[1, 1, 4], MidPoints, BinWidth, threads)
# #         }
#         if (k == 4) { for (t in 1:ny) { for (l in 1:nk) {
#           canlp[t, l] <- dPS(a[t, l], canB[1, 1, 4])
#         } } }
#
#         R <- sum(canll - curll) +
#           sum(canlp - curlp) +
#           dnorm(canbeta[j, k], mnB[j, k], sdB[j, k], log=T) -
#           dnorm(beta[j, k], mnB[j, k], sdB[j, k], log=T)
#
#         if (!is.na(exp(R))) { if (log(runif(1)) < R) {
#           accb[j, k] <- accb[j, k] + 1
#           beta       <- canbeta
#           B          <- canB
#           A          <- canA
#           curll      <- canll
#           curlp      <- canlp
#           FAC        <- canFAC
#         } }
#       } } }
#     }#end thin
#
#     samples[i, , ] <- beta
#
#     #Predictions
#     if (np > 0) {
#       Bp <- array(0, c(ny, np, 5))
#       for (j in 1:5) {
#         Bp[, , j] <- make.B(sp, beta[, j], j)
#       }
#
#       for (t in 1:ny) {
#         Ap         <- a2AP(dwp2, a[t, ], B[1, 1, 4], B[1, 1, 5])
#         U          <- rGEV(np, 1, B[1, 1, 4], B[1, 1, 4])
#         X          <- Ap * U
#         tau        <- exp(-1 / X)
#         yp[i, t, ] <- qGPD(tau, Bp[t, , ], thresh)
#       }
#     }
#
#
#     level <- get.level(olda, cuts)
#     for (j in 1:length(MHa)) {
#       acca[j] <- acca[j] + sum(olda[level == j] != a[level == j])
#       atta[j] <- atta[j] + sum(level == j)
#       if ((i < burn / 2) & (atta[j] > 100)) {
#         if (acca[j] / atta[j] < 0.3) { MHa[j] <- MHa[j] * 0.9 }
#         if (acca[j] / atta[j] > 0.6) { MHa[j] <- MHa[j] * 1.1 }
#         acca[j] <- atta[j] <- 0
#       }
#     }
#
#
#     for (j in 1:nb) { for (k in 1:5) { if ((i < burn / 2) & (attb[j, k] > 50)) {
#       if (accb[j, k] / attb[j, k] < 0.3) { MHb[j, k] <- MHb[j, k] * 0.9 }
#       if (accb[j, k] / attb[j, k] > 0.6) { MHb[j, k] <- MHb[j, k] * 1.1 }
#       accb[j, k] <- attb[j, k] <- 0
#     } } }
#
#
#     #DISPLAY CURRENT VALUE:
#     if (iterplot) { if((i %% update) == 0) {
#       par(mfrow=c(5, 2), mar=c(2, 2, 2, 2))
#       plot(samples[1:i, 1, 1], main="Prob0 Int", type="l")
#       abline(h=0)
#       plot(samples[1:i, 2, 1], main="Prob0 Slope", type="l")
#       abline(h=0)
#       plot(samples[1:i, 1, 2], main="Scale Int", type="l")
#       abline(h=0)
#       plot(samples[1:i, 2, 2], main="Scale Slope", type="l")
#       abline(h=0)
#       plot(samples[1:i, 1, 3], main="Shape Int", type="l")
#       abline(h=0)
#       plot(samples[1:i, 2, 3], main="Shape Slope", type="l")
#       abline(h=0)
#       plot(samples[1:i, 1, 4], main="Alpha Int", type="l")
#       abline(h=0)
#       plot(samples[1:i, 2, 4], main="Alpha Slope", type="l")
#       abline(h=0)
#       plot(samples[1:i, 1, 5], main="BW Int", type="l")
#       abline(h=0)
#       plot(samples[1:i, 2, 5], main="BW Slope", type="l")
#       abline(h=0)
#     } }
#
#     if ((i %% update) == 0) {
#       cat("\t iter", i, "\n")
#     }
#   }
#
#   results <- list(samples=samples, yp=yp)
#   return(results)
# }


#############################################################:
###            OTHER FUNCTION USED IN THE MCMC            ###:
#############################################################:


expit <- function(x) {
  1 / (1 + exp(-x))
}

make.prob <- function(s, beta) {
  eta <- beta[1] + s[, 1] * beta[2] + s[, 2] * beta[3]
  eta[eta > 10] <- 10
  return(expit(eta))
}

make.scale <- function(s, beta) {
  eta <- beta[1] + s[, 1] * beta[2] + s[, 2] * beta[3]
  eta[eta > 10] <- 10
  return(exp(eta))
}

make.shape <- function(s, beta) {
  eta <- beta[1] + s[, 1] * beta[2] + s[, 2] * beta[3]
  return(eta)
}

make.alpha <- function(s, beta) {
  eta <- beta[1] + s[, 1] * beta[2] + s[, 2] * beta[3]
  eta[eta > 10] <- 10
  return(expit(eta))
}

make.gamma <- function(s, beta) {
  eta <- beta[1] + s[, 1] * beta[2] + s[, 2] * beta[3]
  eta[eta > 10] <- 10
  return(eta)
}

make.B <- function(s, beta, type) {
  B <- NA
  if (type == 1) {
    B <- make.prob(s, beta)
  }
  if (type == 2) {
    B <- make.scale(s, beta)
  }
  if (type == 3) {
    B <- make.shape(s, beta)
  }
  if (type == 4) {
    B <- make.alpha(s, beta)
  }
  if (type == 5) {
    B <- make.gamma(s, beta)
  }

  return(B)
}

loglikeparts <- function(y, A, B, thresh) {

  junk <- is.na(y)
  y <- ifelse(junk, thresh, y)

  if (is.matrix(B)) {
    prob <- B[, 1]
    sig <- B[, 2]
    xi <- B[, 3]
    alpha <- B[, 4]
  }

  if (!is.matrix(B)) {
    prob <- B[1]
    sig <- B[2]
    xi <- B[3]
    alpha <- B[4]
  }

  ai <- 1 / alpha
  above <- y > thresh
  above <- ifelse(junk, FALSE, above)
  L1 <- L3 <- expo <- log(1 / prob)^ai

  if (sum(above) > 0) {
    ttt <- xi * (y - thresh) / sig + 1
    ttt <- ifelse(above, ttt, 1)
    t <- (ttt)^(-1 / xi)
    l <- 1 - (1 - prob) * t
    l <- ifelse(above, l, 0.5)
    toohigh <- (xi < 0) & (y > thresh - sig / xi)
    t[toohigh] <- l[toohigh] <- ttt[toohigh] <- 0.1
    L <- log(1 / l)
    L2 <- L^ai
    L3 <- (ai - 1) * log(L) +
          log(1 - prob) +
          (1 + xi) * log(t) -
          log(alpha) -
          log(l) -
          log(sig)
  }
  expo <- ifelse(above, L2, L1)
  expo <- ifelse(junk, 0, expo)

  log <- ifelse(junk, 0, L3)

  results <- list(expo=expo, log=log, above=above)

  return(results)
}

loglike <- function(y, A, B, thresh) {

  junk <- is.na(y)
  y <- ifelse(junk, thresh, y)

  if (is.matrix(B)) {
    prob <- B[, 1]
    sig <- B[, 2]
    xi <- B[, 3]
    alpha <- B[, 4]
  }

  if (!is.matrix(B)) {
    prob <- B[1]
    sig <- B[2]
    xi <- B[3]
    alpha <- B[4]
  }

  ai <- 1 / alpha
  above <- y > thresh
  above <- ifelse(junk, FALSE, above)

  lll <- -A * (log(1 / prob)^ai)

  if (sum(above) > 0) {
    ttt <- xi * (y - thresh) / sig + 1
    ttt <- ifelse(above, ttt, 1)
    t <- (ttt)^(-1 / xi)
    l <- 1 - (1 - prob) * t
    l <- ifelse(above, l, 0.5)
    toohigh <- (xi < 0) & (y > thresh - sig / xi)
    t[toohigh] <- l[toohigh] <- ttt[toohigh] <- 0.1
    L <- log(1 / l)
    logpabove <- -A * (L^ai) +
               log(A) +
               (ai - 1) * log(L) +
               log(1 - prob) +
               (1 + xi) * log(t) -
               log(alpha) -
               log(l) -
               log(sig)
    lll[above] <- logpabove[above]
    lll[toohigh] <- -Inf
  }
  lll <- ifelse(junk, 0, lll)

  return(lll)
}

a2A <- function(FAC, a, alpha) {
  # theta is nxnF
  # s is nFxnt
  # alpha in (0,1)
  W <- FAC^(1 / alpha)
  if (length(a) == 1) {xxx <- W * a}
  if (length(a) > 1) {xxx <- W %*% a}
  return(xxx)
}

a2AP <- function(dw2,a,alpha,gamma){
    rho2 <- exp(gamma)^2
    fac  <- exp(-0.5 * dw2 / rho2)
    FAC  <- sweep(fac, 1, rowSums(fac), "/")
    W    <- FAC^(1 / alpha)
return((W %*% a)^alpha)}

fac2FAC <- function(x, single=F) {
  if (single) {x <- x / sum(x)}
  if (!single) {x <- sweep(x, 1, rowSums(x), "/")}
  return(x)
}

make.fac <- function(dw2, gamma) {
  rho2 <- exp(gamma)^2
  fac <- exp(-0.5 * dw2 / rho2)
  return(fac)
}

rGPD <- function(X, B, thresh) {

  if (is.matrix(B)) {
    prob <- B[, 1]
    sig <- B[, 2]
    xi <- B[, 3]
  }

  if (!is.matrix(B)) {
    prob <- B[1]
    sig <- B[2]
    xi <- B[3]
  }

  U <- exp(-1 / X)
  U2 <- 1 - (U - prob) / (1 - prob)
  Y <- thresh + ifelse(U < prob, 0, sig * (U2^(-xi) -1) / xi)

  return(Y)
}

qGPD <- function(tau, B, thresh){

  if (is.matrix(B)) {
    prob <- B[, 1]
    sig <- B[, 2]
    xi <- B[, 3]
  }

  if (!is.matrix(B)) {
    prob <- B[1]
    sig <- B[2]
    xi <- B[3]
  }

  tau2 <- 1 - (tau - prob) / (1 - prob)

  Y <- thresh + ifelse(tau < prob, 0, sig * (tau2^(-xi) - 1) / xi)
Y}


rGEV <- function(n,mu,sig,xi) {
  tau <- runif(n)
  x <- -1 / log(tau)
  x <- x^(xi) -1
  x <- mu + sig * x / xi

  return(x)
}


ld <- function(u, A, alpha) {
  psi <- pi * u
  c <- (sin(alpha * psi) / sin(psi))^(1 / (1 - alpha))
  c <- c * sin((1 - alpha) * psi) / sin(alpha * psi)
  logd <- log(alpha) - log(1 - alpha) - (1 / (1 - alpha)) * log(A) +
          log(c) - c * (1 / A^(alpha / (1 - alpha)))

  return(exp(logd))
}


ld2 <- function(u, logs, alpha, shift=0, log=T) {

  logs <- logs - shift / alpha
  s <- exp(logs)
  psi <- pi * u
  c <- (sin(alpha * psi) / sin(psi))^(1 / (1 - alpha))
  c <- c * sin((1 - alpha) * psi) / sin(alpha * psi)

  logd <- log(alpha) - log(1 - alpha) - (1 / (1 - alpha)) * logs +
          log(c) - c * (1 / s^(alpha / (1 - alpha))) +
          logs

  return(logd)
}


rPS <- function(n, alpha, MHA=1, iters=10, initA=NULL) {

  A <- matrix(0, iters, n)
  accA <- attA <- 0
  if (!is.null(initA)) { logs <- log(initA) }
  if (is.null(initA)) { logs <- rep(0, n) }
  lll <- rep(0, n)
  B <- runif(n)
  lll <- ld2(B, logs, alpha)
  for (i in 1:iters) {
    can <- rnorm(n, logs, MHA * (1 - alpha))
    ccc <- ld2(B, can, alpha)
    acc <- runif(n) < exp(ccc - lll)
    acc[is.na(acc)] <- F
    logs <- ifelse(acc, can, logs)
    lll <- ifelse(acc, ccc, lll)

    attA <- attA + length(acc)
    accA <- accA + sum(acc)
    if ((i < iters / 2) & (attA > 100)) {
      if (accA / attA < 0.4) {MHA <- MHA * 0.8}
      if (accA / attA > 0.5) {MHA <- MHA * 1.2}
      accA <- attA <- 0
    }

    canB <- rtnorm(B)
    ccc <- ld2(canB, logs, alpha)
    R <- ccc - lll +
         dtnorm(B, canB) -
         dtnorm(canB, B)
    acc <- runif(n) < exp(R)
    acc[is.na(acc)] <- F
    B <- ifelse(acc, canB, B)
    lll <- ifelse(acc, ccc, lll)
    A[i,] <- exp(logs)
  }

  return(A)
}


rtnorm <- function(mn, sd=0.2, fudge=0.001) {
  upper <- pnorm(1 - fudge, mn, sd)
  lower <- pnorm(fudge, mn, sd)
  U <- lower + (upper - lower) * runif(prod(dim(mn)))

  return(qnorm(U, mn, sd))
}

dtnorm <- function(y, mn, sd=0.2, fudge=0.001) {
  upper <- pnorm(1 - fudge, mn, sd)
  lower <- pnorm(fudge, mn, sd)
  l <- dnorm(y, mn, sd, log=T) - log(upper - lower)

  return(l)
}


ECkern <- function(h, alpha, gamma, Lmax=50) {
  dw2 <- rdist(c(0, h), seq(-Lmax, Lmax, 1))
  W <- fac2FAC(make.fac(dw2, gamma))^(1 / alpha)
  for (j in 1:length(h)) {
    h[j]<-sum((W[1, ] + W[j + 1, ])^alpha)
  }

  return(h)
}

npts <- 50
Ubeta <- qbeta(seq(0, 1, length=npts + 1), 0.5, 0.5)
MidPoints <- (Ubeta[-1] + Ubeta[-(npts + 1)]) / 2
BinWidth <- Ubeta[-1] - Ubeta[-(npts + 1)]

dPS <- function(A, alpha, npts=100) {
  l <- -Inf
  if (A > 0) {
    l <- log(sum(BinWidth * ld(MidPoints, A, alpha)))
  }

  return(l)
}


get.level <- function(A, cuts) {
  lev <- A * 0 + 1
  for (j in 1:length(cuts)) {
    lev <- ifelse(A > cuts[j], j + 1, lev)
  }
  return(lev)
}

get.level.sca <- function(A, cuts) {  # minor speedup in scalar case
  lev <- sum(A > cuts) + 1  # stability for when A > max(cuts)
  return(lev)
}

dlognormal <- function(logx, logmu, sig) {
  dnorm(logx, logmu, sig, log=T) - logx
}

