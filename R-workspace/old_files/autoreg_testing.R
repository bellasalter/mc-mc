library(ggplot2)
library(LaplacesDemon)
source("./generic.R")

# big function to do whole test, return t_exp suggested
test_m <- function(a, b, c, steps) {
  states_visited <- c()
  acfs <- c()
  exp_m <- 1 + 2*(a/(1-a))
  states_visited <- asim(get_init_state(), get_next_step, steps)
  my_acfs <- c_est_series(states_visited)
  
  M_upper_bound <- 6000
  triangle_df <- find_m(my_acfs, c, M_upper_bound)
  m_suggest <- triangle_df[triangle_df$poss_Ms >= triangle_df$ct_ints,]
  m_suggest <- m_suggest$poss_Ms[1]
  return(m_suggest)
}

get_next_step <- function(prev) {
  w <- rnorm(n=1, mean=0, sd=1)
  next_val <- a * prev + b * w
  return(next_val)
}

get_init_state <- function() {
  init_sd <- b^2 / (1-a^2)
  init_state <- rnorm(n=1, mean=0, sd=init_sd)
  return(init_state)
}

