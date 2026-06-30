"""
Mutable Struct: Store the Physics describing the physical model
"""
Base.@kwdef mutable struct Physics
    energy_injection::Function = (t::Real) -> 1.0 # Energy yield as a function of time in Msol pc^2 / Myr^2
    mass_injection::Function   = (t::Real) -> 1.0 # Mass yield as a function of time in Msol
    CR_injection::Function     = (t::Real) -> 1.0 # CR Energy yield as a function of time in Msol pc^2 / Myr^2
    cooling::Function          = Λ_default        # Cooling function in Msol pc^2 / Myr^3
    boost_factor::Float64      = 1.0              # Momentum boost factor (see Lancaster et al. 2024) 
end

"""
Mutable Struct: Store the instantaneous injection rates.
"""
Base.@kwdef mutable struct Rates
    Edot::Real    = 0.0      # Instantaneous energy injection rate in Msol pc^2 / Myr^3
    Mdot::Real    = 0.0      # Instantaneous mass injection rate in Msol / Myr
    Edot_CR::Real = 0.0      # Instantaneous CR-energy injection rate in Msol pc^2 / Myr^3
    pdot::Real    = 0.0      # Instantaneous Momentum injection rate in Msol * pc / Myr^2
end

"""
    Environment

Container describing the external physical environment in which the blastwave
evolves.

All fields may be provided either as constants or as functions of
position `x` and time `t`. During setup, constants (scalars, vectors, and matrices) are automatically
lifted to (constant) functions.

# Fields

- `density(x, t)`  
  Ambient density profile in units of Model.n_0 -> scalar valued function.

- `v_ext(x, t)`  
  Background velocity field in km/s -> vector valued function.

- `∇v_ext(x, t)`  
  Velocity gradient tensor in km/s/pc -> matrix valued function.

- `dv_ext_dt(x, t)`  
  Time derivative of background velocity in km/s/Myr -> vector valued function.

- `x_c(t)`  
  Center position of the explosion in pc-> vector valued function.

- `g_ext(x, t)`  
  External gravitational field in (km/s)^2/pc -> vector valued function.

- `P_CR(x, t)`  
  Cosmic-ray pressure background in units of Model.n_0 * (km/s)^2-> scalar valued function.

- `out_of_bounds(u, t)`  
  Boolean condition used to stop integration early -> boolean valued function.
  Here u = (x, y, z, v_x, v_y, v_z, M, E, Eloss, tangential vectors (∂_θ, ∂_φ), tangential vector rates of change (d∂_θ_dt, d∂_φ_dt)) 
  is the instantaneous state vector of the shock surface.

# Notes

- Constants (e.g. `1.0`, `zeros(3)`) are interpreted as spatially and temporally uniform fields.
- Functions are used directly without modification.
- Internally, all inputs are normalized to code units (Msol, pc, Myr) and converted to functions if they are not already.

# Examples

```julia
# Uniform, stationary medium of density 1 amu/cm^3
env = Environment()

# Stratified density
env = Environment(density = (x,t) -> cosh(x[3] / z_s)^-2)

# Moving background flow
env = Environment(v_ext = [0.0, 0.0, 10.0])
```
"""
Base.@kwdef mutable struct Environment
    density       = 1.0          # Density Profile as a function of position and time
    v_ext         = zeros(3)     # Background velocity field as a function of position and time
    ∇v_ext        = zeros(3,3)   # Background velocity gradient field as a function of position and time
    dv_ext_dt     = zeros(3)     # Change of background velocity field as a function of position and time
    x_c           = zeros(3)                        # Center coordinates in pc as a function of time
    g_ext         = zeros(3)     # Background gravitational field as a function of position and time
    P_CR          = 0.0          # Background CR Pressure as a function of position and time
    out_of_bounds = false        # Function to check specific stopping conditions for a given model as a function of time
end

"""
    Model

Container describing the physical parameters of the blastwave model.

Parameters control energy injection, ejecta properties, cosmic-ray feedback,
cooling, and background density normalization.

Some parameters (E_51, M_ej & Λ_22) may be constants or functions of time depending on the chosen
model complexity.

# Fields

- `E_51`  
  Explosion energy as a function of time (or constant) in units of 10^51 erg.

- `M_ej`  
  Ejecta mass as a function of time (or constant) in solar masses.

- `f_CR`  
  Fraction of explosion energy converted into cosmic rays.

- `α_p`  
  Momentum boost factor (see Lancaster et al. 2024).

- `Λ_22`  
  Cooling function normalized to its value at 1e6 K in 10^-22 cm^3 erg/s.
  Providing a non-function value will be interpreted as the default powerlaw cooling function with slope -0.7.

- `n_0`  
  scale density in amu cm⁻³.

# Notes

- `E_51` & `M_ej` can represent:
  - constant values (e.g. single explosion),
  - time-dependent feedback history (e.g. from a star formation history),
  - burst-like (SN) or continuous wind models (via function form).
- Non-function inputs are automatically promoted to constant functions internally.
- In the current implementation, cooling becomes dynamically unimportant below 10^6 K (-> shell formation).

# Examples

```julia
# Single SN with 10^51 erg explosion energy and 1 Msol ejecta mass in a medium with density 1 amu/cm^3 and no cooling
model = Model()

# Time-dependent continuous SN injection with 10^51 erg explosion energy and 1 Msol ejecta mass in a medium with density 1 amu/cm^3
model = Model(E_51 = t -> 1.0 + t / 1.0, M_ej = t -> 1.0 + t / 1.0)

# Burst-like star formation history
model = Model(E_51 = t -> 1.0 + floor(Int, t/1.0), M_ej = t -> 1.0 + floor(Int, t/1.0))

# Single SN with CR injection
model = Model(f_CR = 0.2, α_p = 2.0)
"""
Base.@kwdef mutable struct Model
    # Explosion parameters
    E_51 = 1.0 # Explosion Energy as a function of time in 1e51 erg
    M_ej = 1.0 # Ejecta Mass as a function of time in Msol
    f_CR::Real = 0.0 # Fraction of Explosion Energy converted into non-thermal (Cosmic Ray) Energy
    α_p::Real  = 1.0 # Momentum boost factor (see Lancaster et al. 2024)

    # additional physics parameters
    Λ_22 = 0.0 # Cooling function normalized to its value at 1e6 K in cm^3 erg/s (default is powerlaw cooling function with slope -0.7)

    # environment parameters
    n_0::Real = 1.0  # Midplane density in amu / cm^3
end
