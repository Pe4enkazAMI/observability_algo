using Test

using Logging
logger = Logging.SimpleLogger(stdout, Logging.Debug)
global_logger(logger)

include("../Non-Linear.jl")

@testset "Identifiability of the whole system" begin

    @variables t x1(t) x2(t)
    xdot = [x2, exp(x2)]
    y = x1
    @test is_NL_Observable(xdot, [y], [x1, x2], nothing, true)

    @variables t x(t) y(t) a b c d
    xdot = [(a-b*y)*x, (-c+d*x)*y, a, b, c, d]
    y = x + y
    @test is_NL_Observable(xdot, [y], [x], nothing, true)

    @variables x1 x2
    xdot = [x1^2 + 2*x1*x2, x2^2]
    y = x1 + x2    
    @test !is_NL_Observable(xdot, [y], [x1, x2], nothing, true)

    @variables x1 x2 x3
    xdot = [x1^2 + x2^2 + 2x1 + 2x2 - x1*x3, x3, x1]
    y = x1 + x3*x2
    @test is_NL_Observable(xdot, [y], [x1, x2, x3], nothing, true)

    @variables x1 x2
    xdot = [0, exp(2 * x1) + 2 * exp(x1) * x2 + x2^2]
    y = exp(x1) + x2
    # the actual output should be `false` but for the current version of the algorithm
    # the correct should be something like throwing an exception
    @test is_NL_Observable(xdot, [y], [x1, x2], nothing, true)
end
