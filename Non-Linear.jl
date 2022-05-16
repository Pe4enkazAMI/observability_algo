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
#############################


function build_evaluation!(dict, expression)
     """  
          Input: 
               dict: map where we store our changes
               expression: expression to parse
          Output:
               It's a void function, we expect dict to have all changes after calling our func
     """
     if typeof(expression) <: Int
         return
     end
     println(expression)
     if (typeof(expression) == Symbol) || (operation(expression) in [cos, sin, exp])
         if !(expression in keys(dict))
             dict[expression] = BigInt(rand(1:1000))
         end
         return
     end
     for expr in arguments(expression)
         build_evaluation!(dict, expr)
     end
 end

#############################

function evaluate(expression, dict)
     """
          Input:
               expression: expression to parse
               dict: map of changes
          Output:
               new_expression a copy of an old one, but with substituted values
     """
    if typeof(expression) <: Int
        return expression
    end
    if (typeof(expression) == Symbol) || (operation(expression) in [cos, sin, exp])
        return dict[expression]
    end
    new_tree = [evaluate(new_expression, dict) for new_expression in arguments(expression)]
    return Expr(:call, operation(expression), new_tree...)
 end


########################

function get_upper_rank(jacobian)
     """
          Input:
               jacobian: actually it is a matrix on which we want to know the upper bound for rank
          Output:
               rank: upper bound for rank of jacobian
     """
     var_change = Dict{Any, Any}()
     new_arr = []
     for i in jacobian
          expr = Symbolics.toexpr(i)
          build_evaluation!(var_change, expr)
          push!(new_arr, eval(evaluate(expr, var_change)))
     end
     new_arr = reshape(Vector{Int}(new_arr), size(jacobian))
     return rank(new_arr)
end

######################## Function for computing ans for exact parameter

function get_ans_specific!(Jacobian::Matrix{Num}, Variables::Vector{Num}, Specific::Vector{Num})
     """
          Input:
               Jacobian: Jacobian to procces
               Variables: Variables that we use in Jacobian
               Specific: Specific for which we need to answer: is it identifyable or not
          Output:
               Bool: True if specific is identifyable or False otherwise
     """
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



function is_NL_Observable(sys::Any, output::Any, params::Vector{Num}, specific::Vector{Num}, guarantee::Bool)
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

     if guarantee
          upper_rank = get_upper_rank(jacobian)
          if (upper_rank == ans_rank) && (ans_rank == min(shape[1], shape[2]))
               return true
          elseif (upper_rank != ans_rank) && (ans_rank == min(shape[1], shape[2]))
               return true
          else 
               return false
          end
     end

     if ans_rank == min(shape[1], shape[2])
          return true
     else
          return false
     end
end

println(expand(is_NL_Observable(DX, [y1, y2], [x1, x2, x3], [x1, 0, 0], true)))