# Physical constants
const Msol = 1.989e33
const mH   = 1.67262192e-24            # proton (hydrogen mass)
const X_H  = 0.76                      # primordial hydrodgen fraction
const Y_He = (1 - X_H) / X_H / 4       # primordial abundance of Helium
const μ    = (1 / X_H) * mH            # primordial mean atomic weight
const kB   = 1.380649e-16              # Boltzmann constant
const G    = 6.6740800e-08              # Gravitational constant in cgs

# gas properties
const γ = 5/3                       # adiabatic index
const χ = (γ+1) / (γ-1)             # strong shock compression ratio
const ζ = 2 / (γ+1)                 # strong shock pressure factor
const τ = ζ / χ                     # strong shock temperature factor

# common physical units in cgs
const pc     = 3.086e18
const yr     = 3.155e7
const Myr    = 1e6 * yr
const km_s   = 1e5

# code units (M = Msol, L = pc, t = Myr)
const pc_Myr = pc / Myr                # pc / Myr
const Dunit  = Msol / pc^3 / μ         # Msol/pc^3 in amu/cm^3
const Eunit  = Msol * pc_Myr^2         # Msol (pc / Myr)^2 in cgs
const vunit  = km_s / pc_Myr

# conversion of shock velocity to shock temperature
T_shock(v_shock::Real) = τ * μ * (v_shock * pc_Myr)^2 / kB

# numerical constants
const SMALL = sqrt(eps(1.0))
const TINY = nextfloat(0.0)
const BIG = prevfloat(Inf)

# Sedov Taylor Pressure constants (Laumbach & Probstein 1969)
const LP1  = 0.5 * (γ-1) / (2γ-1)
const LP_E = 1.5 * (γ+1)^2
const LP_k = 3 + 4γ / (γ+1)

# Proportionality constants for self-similar solutions
const ξ_ST = (25/8π * LP1 * LP_E / (1 - 2 * LP1 * (3-LP_k/3)))^(1/5)
const ξ_W  = (25/4π * LP1 * LP_E / (7 - 9 * LP1 * (3-LP_k/3)))^(1/5)
const ξ_MCS = (3/π)^(1/4)
const ξ_MDW = (3/2π)^(1/4)