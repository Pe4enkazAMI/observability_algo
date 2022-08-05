using Arblib
using LinearAlgebra
using Symbolics
using Nemo
include("../Non-Linear.jl")
@variables t x(t) y(t) a b c d
xdot = [(a-b*y)*x, (-c+d*x)*y, a, b, c, d]
z = x + y 
params = [x, y, a, b, c, d]
rand_arr = rand(Int, size(params)[1])
new = is_NL_Observable(xdot, [z], [x, y, a, b, c, d], nothing, false)
new = Symbolics.value.(substitute(new, Dict(params[i] => rand_arr[i] for i in 1:length(params))))

new

matrix = [2 3 4 5 6
          3 4 5 6 7
          4 5 6 7 8
          5 6 7 8 9]

function find_linear_indep(matrix::Any)
    """ Input:
            matrix: matrix of size MxN
        Output:
            ind_array: indices of linear independent set of columns
    """
    m, n = size(matrix)
    mSpace = Nemo.MatrixSpace(ZZ, m, n)
    matr = Nemo.rref(mSpace(matrix))[2]
    ind_array = []
    for i in 1:min(m, n)
        if matr[i, i] != 0
            push!(ind_array, i)
        end
    end
    return ind_array
end

find_linear_indep(matrix)
 



