
library(gsignal)
source("./generic.R")

desh_method <- function(x, y, lims, n) {
  slopes <- c()
  lengths <- c()
  errs <- c()
  xL <- c()
  xR <- c()
  new_weights <- c()
  ints <- c()
  start_time <- Sys.time()
  print(lims)
  for (lhs in 1:lims) {
    for (rhs in 1:lims) {
      if (rhs-lhs > n) {
        mod <- lm(y[lhs:rhs] ~ x[lhs:rhs])
        coeffs <- coef(mod)
        m <- coeffs[2]
        c <- coeffs[1]
        fit_length <- get_fit_length(x[rhs], x[lhs], m)
        fit_err <- get_fit_err(x, y, m, c, lhs, rhs)
        
        slopes <- append(slopes, m)
        lengths <- append(lengths, fit_length)
        
        errs <- append(errs, fit_err)
        xL <- append(xL, x[lhs])
        xR <- append(xR, x[rhs])
        
        ints <- append(ints, c)
      }
    }
  }
  end_time <- Sys.time()
  print(end_time - start_time)
  retVal <- data.frame(slopes, lengths, errs, xL, xR, ints)
  return(retVal)
}

get_desh_weights <- function(results, p, q, r) {
  new_weights <- c()
  #print(results)
  for(row in c(1:nrow(results))) {
    if (results$xL[row] == 0) {
      results$xL[row] <- 0.1
    }
    weight <- (results$lengths[row] ^ p) / (results$errs[row] ^ q)
    new_weights <- append(new_weights, weight)
  }
  print(length(new_weights))
  return(new_weights)
}

get_fit_length <- function(rhs_val, lhs_val,m) {
  coeff <- sqrt(1 + (m^2))
  retVal <- abs(rhs_val - lhs_val)
  retVal <- coeff * retVal
  return(retVal)
}

get_fit_err <- function(x, y, m, c, lhs, rhs) {
  to_sum <- y - (m*x) - c
  to_sum <- to_sum ^ 2
  sum <- sum(to_sum[lhs:rhs])
  num <- sqrt(sum)
  denom <- rhs - lhs
  return(num / denom)
}

get_slopes_sides <- function(results, weights, plot=TRUE) {
  kern <- density(results$slopes, bw = "nrd", adjust = 0.25 , kernel = "gaussian" , weights = weights)
  best_slope <- kern$x[which.max(kern$y)]
  
  #lhs_weighted <- get_lhs_weights(results)
  lhs_kern <- density(results$xL, bw = "nrd", adjust = .1, kernel = "gaussian" ,weights = weights)
  #lhs_kern <- density(lhs_weighted$unique_lhs, bw = "nrd", adjust = 0.25, kernel = "gaussian" ,weights = lhs_weighted$weights)
  #plot(lhs_kern$x, lhs_kern$y, pch = 19)
  best_lhs <- lhs_kern$x[which.max(lhs_kern$y)]
  #best_lhs <- 0
  best_lhs
  
  #rhs_weighted <- get_rhs_weights(results)
  rhs_kern <- density(results$xR, bw = "nrd",adjust = 1, kernel = "gaussian" ,weights = weights)
  best_rhs <- rhs_kern$x[which.max(rhs_kern$y)]
  
  if(plot) {
    plot(kern$x, kern$y, pch = 19)
    plot(rhs_kern$x, rhs_kern$y, pch = 19)
    plot(lhs_kern$x, lhs_kern$y, pch = 19)
    
  }
  
  retVal <- data.frame(best_slope, best_lhs, best_rhs)
  return(retVal)
}


