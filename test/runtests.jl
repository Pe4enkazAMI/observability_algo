using Test

include("../Non-Linear.jl")

@testset "Identifiability of the whole system" begin

    @variables t x1(t) x2(t)
    xdot = [x2, exp(x2)]
    @test is_NL_Observable(xdot, [x1, x2], vec([x1]))

end
