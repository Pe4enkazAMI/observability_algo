using Arblib
using LinearAlgebra
using Symbolics
using Nemo

A = ArbMatrix(Arblib.Random.randn(5), prec = 30)
ATA = A * transpose(A)

B = ArbMatrix([1 24])
lu(B)

@variables x y z
xdo = [exp(x) y]

params = [x, y]
arb_change = ArbMatrix(Arblib.Random.randn(2), prec = 30)
valued_jacobs = Symbolics.value.(substitute(xdo, Dict(params[i] => arb_change[i] for i in 1:lastindex(params))))


params[1, 1]