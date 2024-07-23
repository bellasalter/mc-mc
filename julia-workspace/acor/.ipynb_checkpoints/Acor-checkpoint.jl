module Acor

@assert VERSION >= v"1.1"

export acor

MYDIR = dirname(Base.source_path())
cd(()->run(`make`), MYDIR)

ACORPATH = MYDIR*"/acor"

function acor(X::Vector{Cdouble}, maxlag::Int)

    mean  = [0.0]
    sigma = [0.0]
    tau   = [0.0]

    status = ccall(("acor", ACORPATH),
                   Cint, # int
                   (
                       Ptr{Cdouble}, # double *mean
                       Ptr{Cdouble}, # double *sigma
                       Ptr{Cdouble}, # double *tau
                       Ptr{Cdouble}, # double X[]
                       Cint,         # int L
                       Cint,         # int maxlag
                   ),
                   pointer(mean),
                   pointer(sigma),
                   pointer(tau),
                   pointer(X),
                   length(X),
                   Cint(maxlag))

    (mean=mean[1],sigma=sigma[1],tau=tau[1],status=status)
end

end#module
