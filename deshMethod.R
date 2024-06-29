
library(gsignal)
library("./generic.R")

desh_method <- function(x, y, lims, n, p, q) {
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
        
        weight <- (fit_length ^ p) / (fit_err ^ q)
        new_weights <- append(new_weights, weight)
        ints <- append(ints, c)
      }
    }
  }
  end_time <- Sys.time()
  print(end_time - start_time)
  retVal <- data.frame(slopes, lengths, errs, xL, xR, new_weights, ints)
  return(retVal)
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

get_slopes_sides <- function(results) {
  kern <- density(results$slopes, bw = "nrd", adjust = 0.25 , kernel = "gaussian" , weights = results$new_weights)
  plot(kern$x, kern$y, pch = 19)
  best_slope <- kern$x[which.max(kern$y)]
  best_slope
  
  lhs_weighted <- get_lhs_weights(results)
  lhs_kern <- density(lhs_weighted$unique_lhs, bw = "nrd", adjust = 0.25, kernel = "gaussian" ,weights = lhs_weighted$weights)
  plot(lhs_kern$x, lhs_kern$y, pch = 19)
  best_lhs <- lhs_kern$x[which.max(lhs_kern$y)]
  best_lhs
  
  #rhs_weighted <- get_rhs_weights(results)
  rhs_kern <- density(results$xR, bw = "nrd",adjust = 0.25, kernel = "gaussian" ,weights = results$new_weights)
  plot(rhs_kern$x, rhs_kern$y, pch = 19)
  best_rhs <- rhs_kern$x[which.max(rhs_kern$y)]
  best_rhs
  
  retVal <- data.frame(best_slope, best_lhs, best_rhs)
  return(retVal)
}


