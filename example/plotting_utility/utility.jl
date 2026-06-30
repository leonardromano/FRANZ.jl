############################################################################################################################################
# Functions to assist with the data management in the example plots

"Return indices of snapshots at which time equals roughly certain values."
get_timepoints(time::Union{AbstractRange{<:Real}, Vector{<:Real}}, timepoints::Union{AbstractRange{<:Real}, Vector{<:Real}}) = [findmin(abs.(time .- t))[2] for t in timepoints]

"Return a vector of N logarithmically spaced values between 10^logx0 and 10^logx1."
LogRange(logx0::Real, logx1::Real, N::Integer) =10.0.^(LinRange(logx0, logx1, N))

"Check if a number is NaN, Inf or negative"
isbad(x::Real) = isnan(x) || isinf(x)

hasbad(a::Array{<:Number}) = sum(isbad.(a)) .> 0

"Calculate the gradient on the grid of midpoints (N-1) points"
function gradient_mid(x::Vector{Float64}, y::Vector{Float64})
    Ngrad  = length(x)-1 
    x_mid  = [0.5 * (x[i+1]+x[i]) for i in 1:Ngrad] 
    grad   = [(y[i+1] - y[i]) / (x[i+1]-x[i]) for i in 1:Ngrad]
    return x_mid, grad
end

"Calculate the second-order accurate gradient on the grid"
function gradient_second_order(x::Vector{Float64}, y::Vector{Float64})
    Nx     = length(x)
    grad   = zeros(Nx)
    for i in 2:Nx-1
        λ = ( (x[i+1] - x[i]) / (x[i] - x[i-1]) )^2
        grad[i] = (y[i+1] + (λ-1) * y[i] - λ * y[i-1]) / (x[i+1] + (λ-1) * x[i] - λ * x[i-1])
    end
    return x[2:end-1], grad[2:end-1]
end

"Call the corresponding gradient method"
function gradient(x::Vector{Float64}, y::Vector{Float64}; method::Symbol=:second_order)
    if method == :midpoints
        return gradient_mid(x,y)
    elseif method == :second_order
        return gradient_second_order(x,y)
    else
        println("ERROR: gradient method $method not yet implemented!")
        exit(86)
        return x, y
    end
end

############################################################################################################################################################
# Function wrapper for model arguments

"Wrappers to convert model arguments to functions of position and time, if they are not already functions."
input2func(x::Union{Real, AbstractArray}) = (args...) -> x
input2func(f::Function) = f

############################################################################################################################################################
# Allow multiplying & dividing functions with scalars and functions

Base.:*(f::Function, s::Real)     = (args...) -> s * f(args...)          # scalar on the right
Base.:*(s::Real, f::Function)     = (args...) -> s * f(args...)          # scalar on the left
Base.:*(f::Function, g::Function) = (args...) -> f(args...) * g(args...) # pointwise function multiplication

Base.:/(f::Function, s::Real)     = (args...) -> f(args...) / s         # scalar on the right
Base.:/(s::Real, f::Function)     = (args...) -> s / f(args...)         # scalar on the left
Base.:/(f::Function, g::Function) = (args...) -> f(args...) / g(args...)      # pointwise function division

Base.:max(f::Function, s::Real)     = (args...) -> max(f(args...), s)       # pointwise comparison
Base.:max(s::Real, f::Function)     = (args...) -> max(s, f(args...))       # pointwise comparison
Base.:max(f::Function, g::Function) = (args...) -> max(f(args...), g(args...))    # pointwise comparison

Base.:min(f::Function, s::Real)     = (args...) -> min(f(args...), s)       # pointwise comparison
Base.:min(s::Real, f::Function)     = (args...) -> min(s, f(args...))       # pointwise comparison
Base.:min(f::Function, g::Function) = (args...) -> min(f(args...), g(args...))    # pointwise comparison

############################################################################################################################################################
# Observables

"Compute numerical solution for a shearing bubble in a gravitational field with a vertical velocty gradient and derive shape tensor from ellipsoid approximation."
function compute_deformation(x_shell::Matrix{Vector{Float64}})
    # Get number of timebins
    Ntime = size(x_shell)[1]                           

    # compute time-dependent geometry
    x_c = Matrix{Float64}(undef, Ntime, 3)
    I_s = Array{Float64}(undef, Ntime, 3, 3)

    Threads.@threads for it in 1:Ntime
        # compute center of SNR
        center = mean(x_shell[it, :])

        # compute distance vectors, radii and angles from center
        dx = [xi - center for xi in x_shell[it, :]]
        R  = [sqrt(δx'δx) for δx in dx]
        dV = [r^3 / 3 for r in R]

        # compute geometric moments assuming dΩ = const.
        x_c[it, :] = center + 0.75 * sum(dV .* dx) / sum(dV)

        # compute shape tensor relative to x_c
        dx = [xi - x_c[it, :] for xi in x_shell[it, :]]

        # recompute Volume (elements) relative to center of volume
        dV = [(δx'δx)^1.5 / 3 for δx in dx]

        I_s[it, :, :] = 0.6 * sum(dV .* [δx'δx * I - δx * δx' for δx in dx]) / sum(dV)
    end

    return x_c, I_s
end

"return ellipsoid radii and principal directions of 3x3 matrix"
function get_ellipsoid(M::Matrix{<:Real})
    if hasbad(M)
        return fill(NaN, 3), fill(fill(NaN, 3), 3)
    else
        eig  = eigen(M)
        R    = sqrt.(2.5 * tr(M) .- 5 * eig.values)
        vecs = [eig.vectors[:, i] for i in 1:3]
        return R, vecs
    end
end

"Check if l and be are in the forward direction (positive-x) and if not mirror them."
function forward_frame(l::Real, b::Real)
    if -90 < l < 90
        return l, b
    elseif l < -90
        return l + 180.0, -b
    else
        return l - 180.0, -b
    end
end