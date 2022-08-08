using Symbolics
using LinearAlgebra
using Nemo
using Arblib

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

function find_linear_indep(matrix::Any)
     """ Input:
             matrix: matrix of size MxN
         Output:
             ind_array: indices of linear independent set of columns
     """
     m, n = size(matrix)
 
     if m == 1
         return [1];
     end
 
     mSpace = Nemo.MatrixSpace(QQ, m, n)
     matr = Nemo.rref(mSpace(matrix))[2]
     ind_array = []
     for i in 1:min(m, n)
         if matr[i, i] != 0
             push!(ind_array, i)
         end
     end
     return ind_array
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
     return rank(new_arr), new_arr
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
     random_array = rand(Int, length(Variables)) ## made a random array for substit.
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

########################

function get_lower_rank(matrix)
     matrix = transpose(ArbMatrix(matrix))
     L, U, P = lu(matrix, check = false)
     cols = []
     for i in 1:min(size(U)[1], size(U)[2])
          if U[i, i] != zero(U[i, i])
               push!(cols, i)
          end
     end

     for i in 1:lastindex(cols)
          cols[i] = P[cols[i]]
     end

     return cols
end

########################
function guarantee_func(ans, rank, valued_matr, shape)
     upper_rank, integer_matr = get_upper_rank(ans)

     upper_rank_indep_col = find_linear_indep(integer_matr)
     lower_rank_indep_col = get_lower_rank(valued_matr)

     if (length(upper_rank_indep_col) == length(lower_rank_indep_col) == min(shape[1], shape[2]))
          return true
     else
          return false
     end
end

########################

function is_NL_Observable(sys::Any, output::Any, params::Vector{Num}, specific::Any, guarantee::Bool)
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

     if !(specific === nothing)
          sp_ans = get_ans_specific!(ans, params, specific)
          return sp_ans
     end

     ans_valued = Symbolics.value.(substitute(ans, Dict(params[i] => rand_arr[i] for i in 1:length(params))))
     ans_rank = get_rank(ans_valued)

     if guarantee
          return guarantee_func(ans, ans_rank, ans_valued, shape)
     end

     if ans_rank == min(shape[1], shape[2])
          return true
     else
          return false
     end
end

