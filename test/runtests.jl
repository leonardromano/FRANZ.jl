using Test
using FRANZ

@testset "FRANZ loading" begin
    @test isdefined(FRANZ, :FRANZ)
end

@testset "Sedov Test" begin
    # Test for a Sedov Blastwave in a uniform medium with 
    # Density 1 amu/cc
    # 1e51 erg
    # 1 Msol ejecta
    # After 1 Myr
    t, x, v, n, M, f = numerical_solution([1.0], cosθ=0.0, ϕ=0.0)
    R    = x[1][1]
    R_ST = 81.8523973189277     # analytical solution
    res  = abs(R - R_ST) / R_ST

    # numerical precision test
    @test res < 1e-3
end
