############################################################################################################################################################
# Functions for starting calculation

"""
    numerical_solution(t_Myr; kwargs...)

Compute the numerical evolution of a blastwave shock surface and return its
position, velocity, surface normal, swept-up mass, and finishing status at the requested
output times.

The shock surface can be sampled in two ways:

# Angular grid sampling

Specify one or more values of `cosθ` and `ϕ`:

```julia
numerical_solution(
    t_Myr;
    cosθ = LinRange(-1, 1, 181),
    ϕ    = LinRange(0, 2π, 361)
)
```

Each pair `(cosθ, ϕ)` defines one radial direction. Scalar values may be used
to evolve a single direction.

# Healpix sampling

Specify a Healpix map:

```julia
numerical_solution(
    t_Myr;
    map = HealpixMap{Float64, NestedOrder}(4)
)
```

Each Healpix pixel defines one radial direction.

# Arguments

- `t_Myr` : Output times in Myr.

# Keyword Arguments

- `model=Model()` : Blastwave model.
- `environment=Environment()` : External environment model.
- `cosθ` : Sampling in cosine of polar angle.
- `ϕ` : Sampling in azimuthal angle (radians).
- `map` : Healpix map defining the angular sampling.
- `dlogt=1e-2` : Internal logarithmic timestep.
- `δα=SMALL` : Surface regularization parameter.

# Returns

Returns

```julia
t, x, v, n, M, f
```

where

- `t[it, idir]` : Time in Myr for the corresponding direction.
- `x[it, idir]` : Cartesian position vector in pc.
- `v[it, idir]` : Cartesian velocity vector in km/s (in frame moving along external velocity field).
- `n[it, idir]` : Surface-normal vector in pc^2.
- `M[it, idir]` : Swept-up mass in Msol.
- `f[idir]` : Status flag for the corresponding direction.

The returned times may differ from the requested values if the integration
terminates early (and importantly between different directions!).

# Status Flags

- `-1` : NaN value encountered.
- `0`  : Integration reached the final output time.
- `1`  : Singular surface normal encountered.
- `2`  : User-defined stopping condition triggered.

# Examples

```julia
# Default spherical Sedov blastwave in uniform medium
t, x, v, n, M, f = numerical_solution(t_Myr)

# Single direction
t, x, v, n, M, f = numerical_solution(
    t_Myr;
    cosθ = 1.0,
    ϕ = 0.0
)

# Healpix sampling
t, x, v, n, M, f = numerical_solution(
    t_Myr;
    map = HealpixMap{Float64, NestedOrder}(4)
)
```
"""
function numerical_solution(t_Myr::Vector{<:Real},                          # Snapshot times
                            directions::Vector{<:Tuple{<:Real, <:Real}};    # Sampling of shock surface                                                      
                            model::Model=Model(),                           # model parameters
                            environment::Environment=Environment(),         # model for environment                                                                        
                            dlogt::Real=1e-2, δα::Real=SMALL)               # Time resolution
    # setup problem
    x0, v0, t0, r0, M0, E_0, env, physics, Sedov = setup_problem(model, t_Myr[1], environment=environment)

    # Initial external velocity vector
    v0_ext  = env.v_ext(x0, 0.0)
    
    # Number of time bins
    Ntime = length(t_Myr)

    # output arrays
    x = Array{Vector{Float64}}(undef, Ntime, length(directions)) # Cartesian position in pc
    v = Array{Vector{Float64}}(undef, Ntime, length(directions)) # Cartesian velocity in km/s
    n = Array{Vector{Float64}}(undef, Ntime, length(directions)) # Surface normal in pc^2
    M = Array{Float64}(undef, Ntime, length(directions))         # Mass in Msol
    t = Array{Float64}(undef, Ntime, length(directions))         # Time in Myr
    f = Array{Int}(undef, length(directions))                    # Flag to classify outcomes

    Threads.@threads for idir in eachindex(directions)
        # get unit vectors
        θ, φ = directions[idir]
        e_r, e_φ, e_θ = get_direction(θ, φ)

        # Get initial state
        x_0 = x0 + v0_ext * t0 + r0 * e_r
        v_0 = v0 * e_r

        # get tangential vectors
        ∂_φ = r0 * e_φ  # rescaled with 1 / sin(\theta)
        ∂_θ = r0 * e_θ  # actually ∂_(cosθ), rescaled with -sin(\theta)
        # These rescalings are introduced to remove the singularity at the poles.
        d∂_φ = v0 / r0 * ∂_φ
        d∂_θ = v0 / r0 * ∂_θ

        # Initial state (x, y, z, δv_x, δv_y, δv_z, M, E, Eloss, tangential vectors, tangential vector rates of change)
        u_ini = [x_0..., v_0..., M0, E_0..., ∂_θ..., ∂_φ..., d∂_θ..., d∂_φ...]

        # solve numerical time evolution
        t_out, pos_time, vel_time, dA_time, M_time, f[idir] = numerical_full(t_Myr, t0=t0, u_ini=u_ini, physics=deepcopy(physics), env=env, dlogt=dlogt, δα=δα, Sedov=Sedov)
        
        # assign vectors
        for it in eachindex(t_Myr)
            x[it, idir] = pos_time[it, :]
            v[it, idir] = vel_time[it, :]
            n[it, idir] = dA_time[it, :]
        end

        # assign scalars
        M[:, idir] = M_time
        t[:, idir] = t_out
    end

    # convert velocity to km/s
    v *= pc_Myr / km_s

    return t, x, v, n, M, f
end

"Method to compute numerical solution indexed by (range of) values for cosθ and ϕ."
function numerical_solution(t_Myr::Vector{<:Real};                                                       # Snapshot times
                            model::Model=Model(),                                                        # model parameters
                            environment::Environment=Environment(),                                      # model for environment                                             
                            cosθ::Union{Nothing, AbstractRange, Vector{<:Real}, Real}=nothing,                 
                            ϕ::Union{Nothing, AbstractRange, Vector{<:Real}, Real}=nothing,              # Sampling of shock surface
                            map::Union{Nothing, AbstractHealpixMap}=nothing,                             # Healpix map: shorthand notation for fullsky coverage
                            dlogt::Real=1e-2, δα::Real=SMALL)

    # some checks
    map_exists   = !(isnothing(map))
    cosθ_exists  = !(isnothing(cosθ))
    ϕ_exists     = !(isnothing(ϕ))
    angles_exist = cosθ_exists && ϕ_exists

    # check that input is valid
    if map_exists && (cosθ_exists || ϕ_exists)
        error("Specify either map or cosθ and ϕ. Both at the same time (or any combinations) do not work.")
    elseif !(map_exists) && !(angles_exist)
        error("Need to specify either map or cosθ and ϕ.")
    end

    # precompute initial direction coordinates
    directions = map_exists ? [pix2ang(map, ihp) for ihp in eachindex(map)] : [(acos(cθ), φ) for cθ in cosθ for φ in ϕ]

    return numerical_solution(t_Myr, directions, model=model, environment=environment, dlogt=dlogt, δα=δα)
end