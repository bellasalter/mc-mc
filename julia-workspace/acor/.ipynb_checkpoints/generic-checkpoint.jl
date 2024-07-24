

# does the simulation iters number of times, using the start state(double) and step_function. 
function asim(start_state, step_function, iters = 1000) 
    curr_state = start_state
    path = Float64[]
    for i in 1:iters 
        path = push!(path,curr_state)
        curr_state = step_function(curr_state)
    end
    return path
end

# Computes one value of the autocorrelation function of the visited states with time lag t
# C^(t) in Sokal paper
function c_est(visited_states, t) 
    n = length(visited_states)
    denom = n - abs(t)
    coeff = 1 / denom
    mean_states = mean(visited_states)
    sum = 0
    for i in 0:(denom-1)
        to_add1 = visited_states[i+1] - mean_states
        to_add2 = visited_states[i+abs(t)+1] - mean_states
        to_add = to_add1 * to_add2
        sum = sum + to_add
    end
    return coeff * sum
end

# Computes the entire autocorrelation function for every time lag, stored as a list
function c_est_series(visited_states)
    retVal = Float64[]
    for i in 0:(length(visited_states)-1)
        retVal = push!(retVal, c_est(visited_states, i))
    end
    return retVal
end

# Lambda function for triangular window. 
function lambda(t, M)
    if abs(t) <= M
        retVal = 2(1 - (abs(t) / M))
        #retVal = 1 #used for testing against rectangular window
        return retVal
    else 
        return 0
    end
end 

# Calculates the integrated autocorrelation time given M and the acfs. 
function calc_t_int(M,acfs) 
    #sum = 0
    #n = length(acfs)-1
    n = M
    lower_bound = 1 #1-n
    upper_bound = n-1
    sum = lambda(0, M)
    for t in lower_bound:upper_bound
        p_hat = acfs[t+1] / acfs[1]
        to_add = p_hat * lambda(t,M)
        sum = sum + (2*to_add)
    end
    return (0.5*sum)
end

# Calculates the integrated autocorrelation time corresponding to each M value. 
function calc_t_series(poss_Ms, acfs) 
    retVal = Float64[]
    for n in 1:length(poss_Ms) 
        to_add = calc_t_int(poss_Ms[n],acfs)
        retVal = push!(retVal, to_add)
    end
    return retVal
end

function find_greater(arr1, arr2) 
    for val in 1:length(arr1)
        if arr1[val] > arr2[val]
            return val
        end
    end
    return -1
end

function sokal(acfs, c) 
    poss_Ms = 1:(length(acfs)-1)
    t_series = calc_t_series(poss_Ms, acfs)
    #println(t_series)
    ct_ints = c*t_series
    retVal = find_greater(poss_Ms, ct_ints)
    return poss_Ms[retVal]
end
