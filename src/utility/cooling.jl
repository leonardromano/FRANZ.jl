############################################################################################################################################################
# Default cooling function

"Default cooling function (powerlaw with slope -0.7 normalized to 1e-22 erg cm^3 / s at 1e6 K)."
Λ_default(n::Real, T::Real) = T > 1e4 ? (T / 1e6)^-0.7 : 0.0

############################################################################################################################################################
# Define here your fancy cooling function (normalized to 1e-22 erg cm^3 / s at 1e6 K) and assign it to the model parameters when calling `numerical_solution`