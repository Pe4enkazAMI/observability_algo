using Arblib
using LinearAlgebra
using Symbolics
using Nemo

new = [1 2 0 1
       3 4 0 3]

ns = transpose(ArbMatrix(new, prec = 30))

l, u, p = lu(ns, check = false)

p
ns
