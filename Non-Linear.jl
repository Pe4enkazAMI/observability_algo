using Symbolics
using LinearAlgebra
using Nemo

@variables t x1(t) x2(t) a(t) c(t) d(t)

# Example 1
DX = [exp(x2), x2] ##last 3 eq's are for params "a, c, d"
y1 = x1
y2 = x2



######################## Function for computing Jacobian rank
function get_rank(Matr::Any, threshold = 1e-12)
     """
          Input:
               Matr: matrix -_-
               threshold: -_-
          Output:
               length(S) == rank(Matr)
     """
     U, S, VT = svd(Matr)
     S = filter(n -> n > threshold, S)
     return length(S)
end 
########################


######################## Function for computing ans for exact parameter
function get_ans_specific(Matr::Any, specific::Any, vars::Any)
     """
     Input:
          Matr: Jacobian matrix
          specific: specific vector interesting for us for some reason
          vars: variables of our system
     Output:
          1 if interesting vector is observable.
          0 if not.
     """

    random_array = rand(Float32, length(vars))

    _Only_Here_ = copy(Matr)

    _Only_Here_ = Symbolics.value.(substitute(_Only_Here_, Dict(vars[i] => random_array[i] for i in 1:length(vars))))

    specific = Symbolics.value.(substitute(specific, Dict(vars[i] => random_array[i] for i in 1:length(vars))))

    rkMat = get_rank(_Only_Here_)

    push!(_Only_Here_, specific)

    rkExMat = get_rank(_Only_Here_)
    if rkExMat == rkMat
         return 1
    else
         return 0
    end
end
########################



function is_NL_Observable(sys::Any, viewable::Any, params::Any, specific::Any = nothing)
"""
     Input:
          sys: vector of functions that depicts some system. // STRICTLY in order: normal equation -> dummy parameter func.
          viewable: actually it should be called "output", it is the vector which we observe during something.
          params: vector of parameters obviously.
          specific: not essential parameter, helps to find out the obs. information about some specific vector.
     Output:
          1 if system is observable.
          0 if not.
"""
    rand_arr = rand(Float32, length(params))
     n = length(viewable)
     for i in (n+1):(length(sys)-1 + (n*length(params)))
         push!(viewable, dot(Symbolics.gradient(viewable[i - n], params, simplify = true), sys))
     end

    ans = Symbolics.jacobian(viewable, params, simplify=true)
    shape = size(ans)

    if specific != nothing

        sp_ans = get_ans_specific(ans, specific, params)

        return simplify(det(ans)), sp_ans
    end

    ans_rank = get_rank(Symbolics.value.(substitute(ans, Dict(params[i] => rand_arr[i] for i in 1:length(params)))))

    if ans_rank == min(shape[1], shape[2])
     return 1
    else 
     return 0
    end

end
#smth

println(expand(is_NL_Observable(DX, vec([y1, y2]), [x1, x2])))
