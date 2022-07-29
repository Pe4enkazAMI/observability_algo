using Arblib
using LinearAlgebra
using Symbolics
using Nemo
include("../Non-Linear.jl")
@variables t x(t) y(t) a b c d
xdot = [(a-b*y)*x, (-c+d*x)*y, a, b, c, d]
z = x + y 
params = [x, y, a, b, c, d]
rand_arr = rand(Float32, size(params)[1])
new = is_NL_Observable(xdot, [z], [x, y, a, b, c, d], nothing, false)
new = Symbolics.value.(substitute(new, Dict(params[i] => rand_arr[i] for i in 1:length(params))))


matrix = [2 3 4 5 6
          3 4 5 6 7
          4 5 6 7 8
          5 6 7 8 9]

L, U, P = lu(matrix, check = false)
U

(L*U)
matrix[P, :]



function find_linear_indep(matrix::Any, tol::Real = 1e-6)
    """ Input:
            matrix: matrix of size MxN
        Output:
            ind_array: indices of linear independent set of columns or row (IN A[P, :]!!!!), i haven't decided yet :)
    """
    L, U, P = lu(matrix)
    n, m = size(U);
    ind_array = []
    for i in 1:min(n, m)
        if U[i, i] > tol
            push!(ind_array, i)
        end
    end
    return ind_array 
end

find_linear_indep(matrix)
