############################################################################################################################################################
# Problem Setup

"Setup initial conditons as well as rates and functions describing problem based on given model."
function setup_problem(model::Model, t_first_snapshot::Real; environment::Environment=Environment())
    # define conversion factors to code units
    n_scale = model.n_0 / Dunit # Scale density in Msol / pc^3
    P_scale = n_scale * vunit^2 # Scale pressure code units                   

    # make input environment callable
    density       = input2func(environment.density)
    v_ext         = input2func(environment.v_ext)
    ∇v_ext        = input2func(environment.∇v_ext)
    dv_ext_dt     = input2func(environment.dv_ext_dt)
    x_c           = input2func(environment.x_c)
    g_ext         = input2func(environment.g_ext)
    P_CR          = input2func(environment.P_CR)
    out_of_bounds = input2func(environment.out_of_bounds)

    # copy external environment in code units
    env = Environment()
    env.density       = (x::Vector{<:Real}, t::Real) -> (n_scale * density(x, t)::Real)
    env.v_ext         = (x::Vector{<:Real}, t::Real) -> (vunit   * v_ext(x, t)::Vector{<:Real})
    env.∇v_ext        = (x::Vector{<:Real}, t::Real) -> (vunit   * ∇v_ext(x, t)::Matrix{<:Real})
    env.dv_ext_dt     = (x::Vector{<:Real}, t::Real) -> (vunit   * dv_ext_dt(x, t)::Vector{<:Real})
    env.x_c           = (t::Real)           -> (x_c(t)::Vector{<:Real}) 
    env.g_ext         = (x::Vector{<:Real}, t::Real) -> (vunit^2 * g_ext(x, t)::Vector{<:Real})
    env.P_CR          = (x::Vector{<:Real}, t::Real) -> (P_scale * model.f_CR * P_CR(x, t)::Real)
    env.out_of_bounds = (u::Vector{<:Real}, t::Real) -> (out_of_bounds(u, t)::Bool)

    # normalize cooling function to code units
    Λ_22(n::Real, T::Real)   = model.Λ_22 isa Function ? model.Λ_22(n, T) : model.Λ_22 * Λ_default(n, T)
    Λ_cool = 1e-22 * χ * Λ_22 * Dunit * Msol / μ / Eunit * Myr

    # define various quantities in code units based on model
    E_51(t::Real) = input2func(model.E_51, 0.0)(t)             # all-sky explosion energy in 10^51 erg
    M_ej(t::Real) = input2func(model.M_ej, 0.0)(t)             # all-sky ejecta mass in Msol
    E_in    = 1e51 * (1-model.f_CR) * E_51 / 4π / Eunit        # Explosion Energy in code units
    E_CR_in = 1e51 * model.f_CR * E_51 / 4π / Eunit            # Explosion Energy in Cosmic Rays (CRs)
    M_in    = M_ej / 4π                                        # Ejecta mass in code units
    v_ej    = sqrt ∘ (2 * d_dt(E_in) / max(d_dt(M_in), SMALL)) # Initial ejecta speed in code units  

    # create physics object
    physics = Physics()
    physics.energy_injection = E_in
    physics.mass_injection   = M_in
    physics.CR_injection     = E_CR_in
    physics.cooling          = Λ_cool
    physics.boost_factor     = model.α_p

    # Initial Conditions
    # Coordinates
    x0 = env.x_c(0.0)

    # Get timescale for end of ejecta dominated phase
    n_0 = env.density(x0, 0.0)
    t_ED = 1.9e-4 * (M_ej(SMALL))^(5/6) / max(SMALL, sqrt(E_51(SMALL))) * (n_0 * Dunit)^(-1/3)

    # Initial time radial speed and radius
    t0 = max(min(t_first_snapshot, 0.1 * t_ED), SMALL)
    v0 = sqrt(2 * (E_in(SMALL) - E_in(-SMALL)) / max(M_in(SMALL) - M_in(-SMALL), SMALL))
    r0 = v0 * t0

    # Initial swept up mass and Energy
    M0      = physics.mass_injection(t0) + env.density(x0, 0.0) * r0^3 / 3
    Egain_0 = physics.energy_injection(t0)
    Eloss_0 = physics.cooling(n_0, T_shock(v0)) * t0 * n_0 * M0
    E_CR_0  = physics.CR_injection(t0)
    E_0     = [Egain_0 - Eloss_0, Eloss_0]

    # check if we compute blastwave expansion
    Sedov = E_51(t0) > 0.0

    return x0, v0, t0, r0, M0, E_0, env, physics, Sedov
end