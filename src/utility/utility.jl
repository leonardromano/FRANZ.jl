############################################################################################################################################################
# Utility functions useful for numerical integration, setting up arrays, array manipulations etc.

############################################################################################################################################################
# Unit vectors on sphere

"Return the radial unit vector on a sphere corresponding to a given direction"
get_e_r(θ::Real, φ::Real) = [sin(θ) * cos(φ), sin(θ) * sin(φ), cos(θ)]

"Return the radial unit vector on a sphere corresponding to a given direction"
get_e_θ(θ::Real, φ::Real) = [cos(θ) * cos(φ), cos(θ) * sin(φ), -sin(θ)]

"Return the radial unit vector on a sphere corresponding to a given direction"
get_e_φ(θ::Real, φ::Real) = [-sin(φ), cos(φ), 0]

"Return the unit vectors on a sphere corresponding to a given direction."
get_direction(θ::Real, φ::Real) = get_e_r(θ, φ), get_e_φ(θ, φ), get_e_θ(θ, φ)

############################################################################################################################################################
# Function wrapper for model arguments

"Wrappers to convert model arguments to functions of position and time, if they are not already functions."
input2func(x::Union{Real, AbstractArray}) = (args...) -> x
input2func(x::Real, x0::Real) = (t::Real) -> t >= 0.0 ? x : x0
input2func(f::Function) = f
input2func(f::Function, x0::Real) = (t::Real) -> t >= 0.0 ? f(t) : x0

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
# Gradient of a function of time

"Compute the gradient of a function using evaluation at two nearby points."
d_dt(f::Function; dt::Real=SMALL) = (t::Real) -> (f(t + dt/2) - f(t - dt/2)) / dt