# this file contains functions related to the scaling region algorithm located in Toward automated extraction and characterization of scaling regions in dynamical systems


# New method from paper!
function get_desh_weights(results, p, q)
    new_weights = Float64[]
    #println("weights function")
    for row in 1:nrow(results)
        #println(row)
        weight = (results.lengths[row] ^ p) / (results.errs[row] ^ q)
        push!(new_weights, weight)
    end
    return new_weights
end

# Calculates the fit length
function get_fit_length(rhs_val, lhs_val, m)
    coeff = sqrt(1 + (m^2))
    retVal = abs(rhs_val - lhs_val)
    retVal = coeff * retVal
    return retVal
end

# Calculates the fit error
function get_fit_err(x, y, m, c, lhs, rhs)
    to_repeat = [c]
    to_sum = y - (m * x) - repeat(to_repeat, length(y))
    to_sum = to_sum .^ 2
    sum_val = sum(to_sum[lhs:rhs])
    num = sqrt(sum_val)
    denom = rhs - lhs
    return num / denom
end

# Does the method from the paper
function desh_method(x, y, lims, n)
    slopes = Float64[]
    lengths = Float64[]
    errs = Float64[]
    xL = Float64[]
    xR = Float64[]
    ints = Float64[]
    start_time = now()
    for lhs in 1:lims
        for rhs in 1:lims
            if (rhs - lhs) > n
                data = DataFrame(X = x[lhs:rhs], Y = x[lhs:rhs])
                lin_model = GLM.lm(@formula(Y ~ X), data)
                
                m = GLM.coef(lin_model)[2] 
                c = GLM.coef(lin_model)[1]
                fit_length = get_fit_length(x[rhs], x[lhs], m)
                fit_err = get_fit_err(x, y, m, c, lhs, rhs)

                push!(slopes, m)
                push!(lengths, fit_length)
                push!(errs, fit_err)
                push!(xL, x[lhs])
                push!(xR, x[rhs])
            end
        end
    end
    
    end_time = now()
    println(end_time - start_time)
    retVal = DataFrame(slopes=slopes, lengths=lengths, errs=errs, xL=xL, xR=xR)
    return retVal
end

# Calculates KDE densities for the results(slopes, left hand sides, and right hand sides)
function get_slopes_sides(results, weights; plot_title="", plot_right=false, plot_others=false)
    kern = KernelDensity.kde(results.slopes, weights=weights)  # Adjust bandwidth if needed
    #println(kern)
    #best_slope = findmax(kern).index
    
    lhs_kern = KernelDensity.kde(results.xL, weights=weights)
    #best_lhs = lhs_kern.x[findmax(lhs_kern.y).index]
    
    rhs_kern = KernelDensity.kde(results.xR, weights=weights)
    #best_rhs = rhs_kern.x[findmax(rhs_kern.y).index]
    
    if plot_others
        plot(kern.x, kern.y, label="Slopes")
        plot!(lhs_kern.x, lhs_kern.y, label="LHS")
    end
    
    if plot_right
        plot(rhs_kern.x, rhs_kern.y, label="RHS")
        title!(plot_title)
    end
    
    #return DataFrame(best_slope=best_slope, best_lhs=best_lhs, best_rhs=best_rhs)
    return rhs_kern
end

# Gets the log of the absolute value of the acfs. 
function get_log(acfs::Vector{<:Real})
    return log.(abs.(acfs))
end