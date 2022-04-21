using Symbolics
using LinearAlgebra
using Nemo

@variables t x1(t) x2(t) x3(t)

# Example 1
DX = [exp(x3 + x2), exp(x2 + x3), exp(x2 + x3)] ##last 3 eq's are for params "a, c, d"
y1 = x1
y2 = x2



######################## Function for computing Jacobian rank #UPD
function get_rank(Matr::Any, threshold = 1e-5)
    """
         Input:
              Matr: matrix -_-
              threshold: -_-
         Output:
              length(S) == rank(Matr)
    """
    U, S, VT = svd(Matr)
    S = filter((n) -> n > threshold, S)
    ans = count((i->(i > 0)),S)
    return ans
end
########################


######################## Function for computing ans for exact parameter
function get_ans_specific!(Jacobian::Matrix{Num}, Variables::Vector{Num}, Specific::Vector{Num})
     random_array = rand(Float32, length(Variables)) ## made a random array for substit.

     Spec = copy(Specific)

     before_add = copy(Jacobian)

     before_add = Symbolics.value.(substitute.(before_add, (Dict((Variables[i] => random_array[i]) for i in 1:length(Variables)), )))

     rank1 = get_rank(before_add)
     reshape(Spec, size(Spec))
     before_add = vcat(before_add, Transpose(Spec))

     before_add = Symbolics.value.(substitute.(before_add, (Dict((Variables[i] => random_array[i]) for i in 1:length(Variables)), )))

     rank2 = get_rank(before_add)

     return rank1 == rank2 ? true : false

 end
########################



function is_NL_Observable(sys::Any, output::Any, params::Vector{Num}, specific::Vector{Num})
     """
          Input:
               sys: vector of functions that depicts some system. // STRICTLY in order: normal equation -> dummy parameter func.
               viewable: actually it should be called "output", it is the vector which we observe during something.
               params: vector of parameters obviously.
               specific: not essential parameter, helps to find out the obs. information about some specific vector.
          Output:
               True if system is observable.
               False if not.
     """
     rand_arr = rand(Float32, size(params)[1])

     n = length(output) ## [y1, y2, y3 ... yn] for each yi we find L(vars) - 1 deriv;

     for i in (n+1):(n*length(params))
         push!(output, dot(Symbolics.gradient(output[i - n], params, simplify = true), sys))
     end

     ans = Symbolics.jacobian(output, params, simplify=true)
     shape = size(ans)

     if specific != nothing

          sp_ans = get_ans_specific!(ans, params, specific)

          return sp_ans

     end

     ans_rank = get_rank(Symbolics.value.(substitute(ans, Dict(params[i] => rand_arr[i] for i in 1:length(params)))))

     println(ans)

     if ans_rank == min(shape[1], shape[2])

               return true

     else

               return false

     end
end


println(expand(is_NL_Observable(DX, [y1, y2], [x1, x2, x3], [x1, 0, 0])))
