__precompile__(true)
module FRANZ
    using LinearAlgebra
    using Healpix

    export
        # model struct
        Model,
        Environment,
        # numerical solution
        numerical_full,
        # observables
        numerical_solution

    # constants, types & utility functions
    include("utility/utility.jl")
    include("utility/constants.jl")
    include("utility/structs.jl")
    include("utility/cooling.jl")

    # thin shell problem
    include("numerical/problem.jl")
    include("numerical/start_calculation.jl")
    include("numerical/solver.jl")
end