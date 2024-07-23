
library(gsignal)
source("./generic.R")

# Performs the algorithm from the paper. 
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
  for(row in c(1:nrow(results))) {
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

# Uses KDEs for each slope and side using the weights, returning a data frame of the recommended values. 
get_slopes_sides <- function(results, weights, plot_title = "",plot_right=TRUE, plot_others=FALSE) {
  kern <- density(results$slopes, bw = "nrd", adjust = 0.25 , kernel = "gaussian" , weights = weights)
  best_slope <- kern$x[which.max(kern$y)]
  
  lhs_kern <- density(results$xL, bw = "nrd", adjust = .1, kernel = "gaussian" ,weights = weights)
  best_lhs <- lhs_kern$x[which.max(lhs_kern$y)]
  
  rhs_kern <- density(results$xR, bw = "nrd",adjust = 1.5, kernel = "gaussian" ,weights = weights)
  best_rhs <- rhs_kern$x[which.max(rhs_kern$y)]
  
  if(plot_others) {
    plot(kern$x, kern$y, pch = 16)
    plot(lhs_kern$x, lhs_kern$y, pch = 16)
  }
  
  if(plot_right) {
    plot(rhs_kern$x, rhs_kern$y, pch = 16)
    title(plot_title)
  }
  
  retVal <- data.frame(best_slope, best_lhs, best_rhs)
  return(retVal)
}


