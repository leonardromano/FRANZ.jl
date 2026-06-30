##########################################################################################
# Numerical solutions to thin shell model

##########################################################################################
# Differential equation

"Compute the acceleration for a given state."
function get_acceleration(x::Vector{<:Real}, v::Vector{<:Real}, t::Real, x_ref::Vector{<:Real}, δM::Real, δE::Real, dA::Vector{<:Real}, dt::Real, Sedov::Bool, 
                          rates::Rates, physics::Physics, env::Environment)
    # unit normal vector
    n  = normalize(dA)

    # external velocity field
    v_ext = env.v_ext(x, t)
    dx_dt = v + v_ext

    # local density
    ρ = env.density(x, t)

    # velocity squared
    v2 = v'v

    # get volume element
    δV = dot(x-x_ref, dA) / 3

    # Mass growth rate
    dM_dt = rates.Mdot + ρ * max(dA'v, 0) # Inside of Bubble assumed to be empty (dM_dt = 0 for inward moving shell)

    # drag force (due to shell inertia)
    F_drag = - dM_dt * v

    # radial force (wind injection, pressure gradients, etc.)
    F_rad = rates.pdot * n

    # compute changes due to pressure gradient / cooling
    if Sedov
        # Pressure gradient force
        Δρ = env.density(x + 0.5 * dt * dx_dt, t+0.5dt) - env.density(x - 0.5 * dt * dx_dt, t-0.5dt)
        r  = norm(x - x_ref)
        Δr = norm(x - x_ref + (dx_dt - env.v_ext(x_ref, t)) * dt) - r
        k_ρ   = - r / ρ * Δρ / Δr * max(sign(v'n), 0) #  only if expansion is outward (!)

        F_rad += LP1 * (LP_E * δE / (3 * δV) + 3 * ρ * v2 + (k_ρ - LP_k) * δM * v2 / (3 * δV)) * dA
    end

    # add force contribution from CRs to radial force
    F_rad += max((physics.CR_injection(t) / δV / 3 - env.P_CR(x, t)), 0) * dA # 1/3 here comes from P_CR = (γ_CR - 1) * ϵ_CR with γ_CR = 4/3

    # acceleration due to change in external Velocity (shear)
    dv_ext_dt = env.∇v_ext(x, t) * dx_dt + env.dv_ext_dt(x, t)

    # combine to get accelerations (radial expansion + external acceleration + external acceleration due to change in external velocity)
    return F_rad / δM + F_drag / δM + env.g_ext(x, t) - dv_ext_dt

end

"Equation of motion for shell being deformed by galactic shear"
function du_dt(u::Vector{<:Real}, t::Real, dt::Real; rates::Rates=Rates(), physics::Physics=Physics(), env::Environment=Environment(), 
               Sedov::Bool=false, δα::Real=SMALL)::Vector{<:Real}
    # get position vector
    x   = u[1:3]
    x_c = env.x_c(t)

    # get surface normal and tangential vectors
    ∂_θ = u[10:12]
    ∂_φ = u[13:15]
    dA = cross(∂_θ, ∂_φ)

    # compute velocity vector
    vel     = u[4:6]
    dx_dt   = vel + env.v_ext(x, t)
    v2      = vel'vel
    v       = sqrt(v2)

    # get Mass & density
    M = u[7]
    ρ = env.density(x, t)
    
    # Energy gain and loss rates
    dEgain_dt = Sedov ? rates.Edot : 0.0
    dEloss_dt = Sedov ? physics.cooling(ρ, T_shock(v)) * ρ * M : 0.0
    dE_dt     = dEgain_dt - dEloss_dt + M * dx_dt'env.g_ext(x, t)

    # net range of energy due to CRs
    dE_CR_dt = rates.Edot_CR

    # change in Mass
    Mdot_sw = ρ * max(dA'vel, 0) # Inside of Bubble assumed to be empty (dM_dt = 0 for inward moving shell)
    dM_dt   = Mdot_sw + rates.Mdot

    # acceleration
    acc = get_acceleration(x, vel, t, x_c, M, u[8], dA, dt, Sedov, rates, physics, env)

    # get tangential vector rates of change
    # local velocity gradients
    d∂_θ = u[16:18]
    d∂_φ = u[19:21]

    # external velocity gradient 
    ∇v_ext = env.∇v_ext(x, t)
    d∂_θ_dt = d∂_θ + ∇v_ext * ∂_θ
    d∂_φ_dt = d∂_φ + ∇v_ext * ∂_φ

    # get rate of change of velocity gradients
    δdA = δα * norm(dA) * normalize(∂_θ)
    da_dθ = ( get_acceleration(x + δα * ∂_θ, vel + δα * d∂_θ, t, x_c, M, u[8], dA + δdA, dt, Sedov, rates, physics, env) 
            - get_acceleration(x - δα * ∂_θ, vel - δα * d∂_θ, t, x_c, M, u[8], dA - δdA, dt, Sedov, rates, physics, env)) / (2*δα)

    δdA = δα * norm(dA) * normalize(∂_φ)
    da_dφ = ( get_acceleration(x + δα * ∂_φ, vel + δα * d∂_φ, t, x_c, M, u[8], dA + δdA, dt, Sedov, rates, physics, env) 
            - get_acceleration(x - δα * ∂_φ, vel - δα * d∂_φ, t, x_c, M, u[8], dA - δdA, dt, Sedov, rates, physics, env)) / (2*δα)

    return [dx_dt..., acc..., dM_dt, dE_dt, dEloss_dt, d∂_θ_dt..., d∂_φ_dt..., da_dθ..., da_dφ...] 
end

##########################################################################################
# RK4-Integrator

"Integrate the equation of motion for one timestep using the Runge-Kutta-Algorithm"
function integrate_EoM(u0::Vector{<:Real}, t0::Real, t1::Real; physics::Physics=Physics(), env::Environment=Environment(), Sedov::Bool=false, δα::Real=SMALL)
    # get differentials
    dt = t1 - t0

    # estimate instantaneous rates from yield functions
    rates = Rates()
    rates.Edot    = (physics.energy_injection(t1) - physics.energy_injection(t0)) / dt
    rates.Mdot    = (physics.mass_injection(t1)   - physics.mass_injection(t0)) / dt
    rates.Edot_CR = (physics.CR_injection(t1)     - physics.CR_injection(t0)) / dt
    rates.pdot    = sqrt(2 * rates.Edot * rates.Mdot)
    rates.pdot   *= Sedov ? 1.0 : physics.boost_factor # wind momentum injection with boost factor (rapidly cooling wind; Lancaster et al. 2024)

    # compute differentials at different stages of the RK4 algorithm
    k1 = du_dt(u0,              t0,         dt, rates=rates, physics=physics, env=env, Sedov=Sedov, δα=δα)
    k2 = du_dt(u0 + 0.5dt * k1, t0 + 0.5dt, dt, rates=rates, physics=physics, env=env, Sedov=Sedov, δα=δα)
    k3 = du_dt(u0 + 0.5dt * k2, t0 + 0.5dt, dt, rates=rates, physics=physics, env=env, Sedov=Sedov, δα=δα)
    k4 = du_dt(u0 + dt * k3,    t0 + dt,    dt, rates=rates, physics=physics, env=env, Sedov=Sedov, δα=δα)

    # combine to get total differential
    du = dt * (k1 + 2k2 + 2k3 + k4) / 6

    # catch singular tangent vectors
    if dot(u0[10:12], u0[10:12] + du[10:12]) < 0
        du[10:12] = -u0[10:12]
    elseif dot(u0[13:15], u0[13:15] + du[13:15]) < 0
        du[13:15] = -u0[13:15]
    end

    return u0 + du
end

"Compute timestep based on a number of timestep criteria."
function get_timestep(t0::Real, u0::Vector{<:Real}, dlogt::Real; physics::Physics=Physics(), env::Environment=Environment(), Sedov::Bool=false)
    # most basic timestep criterium
    dt = t0 * dlogt

    # compute velocity vector
    vel  = u0[4:6]
    v2   = vel'vel
    v    = sqrt(v2)

    # get density and radius squared
    ρ  = env.density(u0[1:3], t0)

    # get mass changing rate
    dM_dt_in = (physics.mass_injection(t0 + dt) - physics.mass_injection(t0)) / dt
    dM_dt = ρ * max(vel'cross(u0[10:12], u0[13:15]), 0) + dM_dt_in

    # Properly resolve velocity changes due to mass changes
    dt = min(dt, 0.1 * u0[7] / dM_dt)

    # resolve cooling during Sedov-phase
    if Sedov
        dt = min(dt, 0.1 * u0[8] / abs(physics.cooling(ρ, T_shock(v))) / ρ / u0[7])
    end

    # resolve changes in external velocity
    dt = min(dt, 0.1 / norm(env.∇v_ext(u0[1:3], t0)))

    # resolve changes in external gravitational field
    dt = get_gravitational_timestep(dt, u0[1:3], u0[4:6] + env.v_ext(u0[1:3], t0), t0, env.g_ext)

    return dt
end

"iteratively reduce the timestep to avoid large variations in the gravitational field BUT only if it significantly affects the dynamics."
function get_gravitational_timestep(dt::Real, x::Vector{<:Real}, dx_dt::Vector{<:Real}, t::Real, g_ext::Function)
    Δt  = dt
    g_x = g_ext(x, t)
    TOL_g = 0.01 * norm(g_x)
    TOL_v = 1e-4 * norm(dx_dt)

    while norm(g_ext(x + dx_dt * Δt + 0.5 * g_x * Δt^2, t + Δt) - g_x) > TOL_g && TOL_g * Δt > TOL_v
        Δt *= 0.5
    end

    return Δt
end

##########################################################################################
# Solvers

"Integrate EoM to get numerical trajectory for given model."
function numerical_full(t_Myr::Vector{<:Real}; t0::Real=0.0, u_ini::Vector{<:Real}=rand(21), # state variables and initial time
                        physics::Physics=Physics(),                                          # physics model
                        env::Environment=Environment(),                                      # environment model
                        dlogt::Real=1e-2, δα::Real=SMALL, Sedov::Bool=true)                  # Numerical criteria

    # get number of snapshots
    N = length(t_Myr)
    
    # get output objects
    t_out = Vector{Float64}(undef, N)
    x_out = Matrix{Float64}(undef, N, 3)
    v_out = Matrix{Float64}(undef, N, 3)
    n_out = Matrix{Float64}(undef, N, 3)
    M_out = Vector{Float64}(undef, N)

    # Initial state (x, y, z, v_x, v_y, v_z, M, E, Eloss, tangential vectors, tangential vector rates of change)
    u0 = copy(u_ini)
    u1 = copy(u_ini)

    # initial time
    t0 = t0
    t1 = t0

    # current snapshot
    isnapshot = 1

    # completion flag:
    # -1 -> NaN value encountered
    #  0 -> default (computation reached final time)
    #  1 -> singular surface normal
    #  2 -> out-of-bounds (user specified stopping condition)
    flag = 0

    # integrate forward in time
    while isnapshot <= N
        # check if we need to write a snapshot
        if t1 >= t_Myr[isnapshot]

            if t0 == t1
                x_out[isnapshot, :] = u1[1:3]
                v_out[isnapshot, :] = u1[4:6]
                n_out[isnapshot, :] = cross(u1[10:12], u1[13:15])
                M_out[isnapshot]    = u1[7]
            else
                # linear interpolation
                x_0 = u0[1:3]
                x_1 = u1[1:3]
                v_0 = u0[4:6]
                v_1 = u1[4:6]
                n_0 = cross(u0[10:12], u0[13:15])
                n_1 = cross(u1[10:12], u1[13:15])

                x_out[isnapshot, :] = (x_0 * (t1 - t_Myr[isnapshot]) + x_1 * (t_Myr[isnapshot] - t0)) / (t1 - t0)
                v_out[isnapshot, :] = (v_0 * (t1 - t_Myr[isnapshot]) + v_1 * (t_Myr[isnapshot] - t0)) / (t1 - t0)
                n_out[isnapshot, :] = (n_0 * (t1 - t_Myr[isnapshot]) + n_1 * (t_Myr[isnapshot] - t0)) / (t1 - t0)
                M_out[isnapshot]    = (u0[7] * (t1 - t_Myr[isnapshot]) + u1[7] * (t_Myr[isnapshot] - t0)) / (t1 - t0)
            end

            t_out[isnapshot] = t_Myr[isnapshot]

            # update snapshot index
            isnapshot += 1

            if isnapshot > N
                break
            end
        end

        # check if multiple streamlines are converging or volume is zero
        if norm(u1[10:12]) == 0.0 || norm(u1[13:15]) == 0.0
            # we don't have a good way to obtain a solution beyond this point

            # capture the last accurate state
            t_out[isnapshot]    = t0
            x_out[isnapshot, :] = u0[1:3]
            v_out[isnapshot, :] = u0[4:6]
            n_out[isnapshot, :] = cross(u0[10:12], u0[13:15])
            M_out[isnapshot]    = u1[7]

            # set all later values to NaN
            t_out[(isnapshot+1):N]    .= NaN
            x_out[(isnapshot+1):N, :] .= NaN  
            v_out[(isnapshot+1):N, :] .= NaN 
            n_out[(isnapshot+1):N, :] .= NaN
            M_out[(isnapshot+1):N]    .= NaN
            flag = 1

            break
        end

        # check if out-of-bounds condition is met
        if env.out_of_bounds(u1, t1)
            # reached user-defined stopping condition

            # capture the last accurate state
            t_out[isnapshot]    = t0
            x_out[isnapshot, :] = u0[1:3]
            v_out[isnapshot, :] = u0[4:6]
            n_out[isnapshot, :] = cross(u0[10:12], u0[13:15])
            M_out[isnapshot]    = u1[7]

            # set all later values to NaN
            t_out[(isnapshot+1):N]    .= NaN
            x_out[(isnapshot+1):N, :] .= NaN  
            v_out[(isnapshot+1):N, :] .= NaN 
            n_out[(isnapshot+1):N, :] .= NaN
            M_out[(isnapshot+1):N]    .= NaN
            flag = 2

            break
        end

        if sum(isnan.(u1)) > 0
            println("Encountered NaN in model: $u_ini")
            println("t0 = $t0, u0 = $u0")
            println("t1 = $t1, u1 = $u1")

            # capture the last non-NaN state
            t_out[isnapshot]    = t0
            x_out[isnapshot, :] = u0[1:3]
            v_out[isnapshot, :] = u0[4:6]
            n_out[isnapshot, :] = cross(u0[10:12], u0[13:15])
            M_out[isnapshot]    = u1[7]

            # set all later values to NaN
            t_out[(isnapshot+1):N]    .= NaN
            x_out[(isnapshot+1):N, :] .= NaN  
            v_out[(isnapshot+1):N, :] .= NaN 
            n_out[(isnapshot+1):N, :] .= NaN
            M_out[(isnapshot+1):N]    .= NaN
            flag = -1

            break
        end

        # update old state
        t0 = t1
        u0 = u1

        # assign increase in t
        dt = get_timestep(t0, u0, dlogt, physics=physics, env=env, Sedov=Sedov)

        # check that we don't skip any snapshots
        if isnapshot < N
            dt = min(dt, 0.9 * (t_Myr[isnapshot+1] - t0))
        end

        # update new state
        t1 += dt
        u1  = integrate_EoM(u0, t0, t1, physics=physics, env=env, Sedov=Sedov, δα=δα)

        # Sedov phase ends once cooling losses start to become dominant
        if Sedov && (u1[9] >= 0.1 * physics.energy_injection(t1))
            Sedov = false
        end
    end

    return t_out, x_out, v_out, n_out, M_out, flag
end