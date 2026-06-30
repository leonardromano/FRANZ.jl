using Test
using FRANZ

@testset "FRANZ loading" begin
    @test isdefined(FRANZ, :FRANZ)
end
