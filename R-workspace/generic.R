# Generic interface for stochastic processes

library(ggplot2)
library(LaplacesDemon)

# Applies the window lambda to the acfs given. 
apply_window <- function(lambda, acf_all, M) {
  c_new <- c()
  n <- 0
  N <- 1000
  for (acf in acf_all$acf) {
    rat <- lambda(n, M)
    to_add <- rat * acf
    c_new <- append(c_new, to_add)
    n <- n+1
  }
  as_df <- data.frame(c(1:steps), c_new, acf_all$acf)
  colnames(as_df) <- c("Lag", "C_new", "C_orig")
  return(as_df)
}

# General layout of stochastic process. uses start_state and get_next_step to simulate iters number of steps of a stochastic process. 
#   Input: start_state(numeric), get_next_step(function of previous step), and iters(numeric)
#   Returns: path of stochastic process
asim <- function(start_state, get_next_step, iters = 1000) {
  curr_state <- start_state
  path <- c()
  for(i in c(1:iters)) {
    path <- append(path,curr_state)
    curr_state <- get_next_step(curr_state)
  }
  #acf_vals <- acf(path, plot=FALSE)$acf
  #sum_acf <- mean(acf_vals)
  #acfs <<- append(acfs, sum_acf)
  return(path)
}

# Estimates C(t) for the given visited states and time lag. 
#   Returns: C(t)
c_est <- function(visited_states, t) {
  n <- length(visited_states) 
  const <- n - abs(t)
  coeff <- 1 / const
  #sprintf("const: %f, coeff: %f", const, coeff)
  mean_states <- mean(visited_states)
  sum <- 0
  for (i in 0:(const-1)) {
    to_add1 <- visited_states[i+1] - mean_states
    to_add2 <- visited_states[i+abs(t)+1] - mean_states
    to_add <- to_add1 * to_add2
    #print(sprintf(" adding first: %f, second: %f at i=%f", to_add1, to_add2, i))
    sum <- sum + to_add
  }
  return(coeff * sum)
}

# Estimates C(t) for each value in the list. 
#   Returns: a dataframe with one column, titled acf
c_est_series <- function(visited_states) {
  retVal <- c()
  for(i in 0:length(visited_states)) {
    print(sprintf("I ======== %f", i))
    retVal <- append(retVal, c_est(visited_states, i))
  }
  ret_acfs <- data.frame(retVal)
  colnames(ret_acfs) <- c("acf")
  return(ret_acfs)
}


# Estimates the mean state of generated path with 100 simulations. 
#   Inputs: start_gen(function generating the start state based on stationary dist. ), get_next_step, num_steps, and states_visited(list to store the path) 
#   Returns: variance of the mean
est_sim <- function(start_gen, get_next_step, num_steps, states_visited) {
  results <- c()
  last_states <- c()
  for(i in c(1:1000)) {
    init_state <- start_gen()
    states_visited <<- asim(init_state, acfs, get_next_step, num_steps)
    last_states <- append(last_states, states_visited[length(states_visited)])
    res <- mean(states_visited)
    results <- append(results, res)
  }
  var_phi <- sd(results)
  ac <- acf(results, type = "covariance", plot=FALSE)
  return(var_phi)
}

# Plots the estimated c*tau_int(M) vs. M for comparison to determine autocorrelation time. 
#   Returns: the dataframe used to plot. 
find_m <- function(acf_all, c, M_upper_bound) { 
  poss_Ms <- c(1:M_upper_bound)
  print(poss_Ms)
  ct_ints <- c*t_int_series(poss_Ms, acf_all)
  #print(length(ct_ints))
  to_plot <- data.frame(poss_Ms, ct_ints)
  return(to_plot)
}

# Serves as the lambda function for windowing the acf. 
#   Inputs: t, the lag(numeric), M, the autocorrelation time(numeric), and type(string), default t for triangle
#   Returns: lambda(t)
lambda <- function(t, M, type="t") {
  t <- abs(t)
  if(type == "r") {
    if (abs(t) <= M) {
      retVal <- 1
      print("less than m")
      return(retVal)
    } else {
      print(sprintf("%f > %f", t, M))
      return(0)
    }
  }
  if (t <= M) {
    retVal <- 1 - (t/M)
    return(retVal)
  } else {
    return(0)
  }
}

# Finds tau_int for the given N and acfs based on formula from Sokal
t_int <- function(M, acf_all) {
  start_ind <- (1-length(acf_all$acf)) 
  end_ind <- length(acf_all$acf)-1
  sum <- 0
  print(sprintf("T _ INT FUNCTION n=%f, m=%f, start=%f, end=%f", length(acf_all), M, start_ind, end_ind))
  for (t in start_ind:end_ind) {
    p_hat <- acf_all$acf[abs(t)+1] / acf_all$acf[1] 
    new_term <- lambda(t, M,  type="t")*p_hat
    sum <- sum + new_term
    print(sprintf("to_add1: %f, 2: %f, adding %f, lambda=%f", acf_all$acf[abs(t)+1], acf_all$acf[1], new_term, lambda(t, M,  type="r")))
  }
  #retVal <- 0.5 *sum
  retVal <- sum
  return(retVal)
}

# Finds tau_int for each item in a list. 
t_int_series <- function(n_series, acf_all) {
  retVal <- c()
  for(n in c(1:length(n_series))) {
    to_add <- t_int(n_series[n],acf_all)
    retVal <- append(retVal, to_add)
  }
  return(retVal)
}

get_lhs_weights <- function(results) {
  weights <- c()
  unique_lhs <- unique(results$xL)
  for(lhs in unique_lhs) {
    selected <- results[results$xL == lhs,]
    weight <- max(selected$new_weight)
    weights <- append(weights, weight)
  }
  retVal <- data.frame(unique_lhs, weights)
  return(retVal)
}

get_rhs_weights <- function(results) {
  weights <- c()
  unique_rhs <- unique(results$xR)
  for(rhs in unique_rhs) {
    selected <- results[results$xR == rhs,]
    weight <- max(selected$new_weight)
    weights <- append(weights, weight)
  }
  retVal <- data.frame(unique_rhs, weights)
  return(retVal)
}


# acfs turn negative, so find first time and discard those data points
get_log <- function(acfs) {
  #first_negative <- which(acfs <= 0)[1]
  return(log(abs(acfs)))
}

generate_piecewise <- function(x) {
  y <- rep(0, length(x))
  for(i in 1:length(x)) {
    to_add <- 10 * exp(-0.5*x[i])
    to_add <- to_add * sin(5*x[i])
    if(x[i] < 1) {
      y[i] <- 0 + to_add
    } else if(x[i] > 9) {
      y[i] <- 40 + to_add
    } else {
      y[i] <- (5*x[i]) - 5 + to_add
    }
  }
  return(y)
}


get_table <- function(slopes) {
  unique_vals <- unique(slopes)
  retVal <- rbind(label=unique_vals, count=sapply(unique_vals,function(x)sum(slopes==x)))
  return(retVal)
}

get_weights_endpoints <- function (vals, results, col_select) {
  weights <- c()
  for(v in vals$val) {
    v_scenarios <- results[results[col_select] == v,]
    weight <- sum(v_scenarios$weights) * vals$Freq
    weights <- append(weights, weight)
  }
  return(data.frame(vals$val,weights))
}





