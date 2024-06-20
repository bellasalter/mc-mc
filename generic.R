# Generic interface for stochastic processes

library(ggplot2)
library(LaplacesDemon)

# General layout of stochastic process. uses start_state and get_next_step to simulate iters number of steps of a stochastic process. 
#   Input: start_state(numeric), get_next_step(function of previous step), and iters(numeric)
#   Returns: path of stochastic process
sim_steps <- function(start_state, get_next_step, iters = 1000) {
  curr_state <- start_state
  path <- c()
  for(i in c(1:iters)) {
    path <- append(path,curr_state)
    curr_state <- get_next_step(curr_state)
  }
  #acf_vals <- acf(path, plot=FALSE)$acf
  sum_acf <- mean(acf_vals)
  #acfs <<- append(acfs, sum_acf)
  return(path)
}

# Estimates the mean state of generated path with 100 simulations. 
#   Inputs: start_gen(function generating the start state based on stationary dist. ), get_next_step, num_steps, and states_visited(list to       store the path) 
#   Returns: variance of the mean
est_sim <- function(start_gen, get_next_step, num_steps, states_visited) {
  results <- c()
  last_states <- c()
  for(i in c(1:1000)) {
    init_state <- start_gen()
    states_visited <<- sim_steps(init_state, acfs, get_next_step, num_steps)
    last_states <- append(last_states, states_visited[length(states_visited)])
    res <- mean(states_visited)
    results <- append(results, res)
  }
  var_phi <- sd(results)
  ac <- acf(results, type = "covariance", plot=FALSE)
  return(var_phi)
}





