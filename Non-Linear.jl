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
function find_linear_indep(matrix::Any)
     """ Input:
            matrix: matrix of size MxN
        Output:
            ind_array: indices of linear independent set of columns
    """
    m, n = size(matrix)
    ind_col = []
    MatSpace = Nemo.MatrixSpace(QQ, m, n)
    rref_mat = Nemo.rref(MatSpace(matrix))[2]
    for i in 1:m
        for j in 1:n
            if rref_mat[i, j] != 0
                push!(ind_col, j)
                break
            end
         end
    end
    return ind_col
end

########################

function get_lower_rank(indep_cols, jacobian, params)
     arb_change = ArbMatrix(Arblib.Random.randn(length(params)), prec = 30)
     n, m = size(jacobian)
     valued_jacob = ArbMatrix(n, m)
     for i in 1:n
          for j in 1:m
               valued_jacob[i, j] = Symbolics.value(substitute(jacobian[i,j], Dict(params[k] => arb_change[k] for k in 1:lastindex(params))))
          end
     end
     @debug "Interval matrix after specialization $valued_jacob with precision $(valued_jacob.prec)"

     ATA = transpose(valued_jacob[:,indep_cols]) * valued_jacob[:,indep_cols]
     @debug "the ATA matrix is $ATA"
     try 
         jac_lu = lu(ATA) 
         @debug "Resulting LU is $jac_lu"
     catch e
         @debug "Caught $e in LU"
         return false
     end
     return true
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
     @debug "Resulting evaluation $var_change"
     new_arr = reshape(Vector{Int}(new_arr), size(jacobian))
     return rank(new_arr), find_linear_indep(new_arr)
end
########################

function guarantee_func(jacobian, params)
     rank, indep_cols = get_upper_rank(jacobian)
     @debug "Upper bound for the rank is $rank, conjectured independent columns are $indep_cols"
     return get_lower_rank(indep_cols, jacobian, params)
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

     if guarantee
          return guarantee_func(ans, params) == min(shape[1], shape[2])
     end

     rand_arr = rand(Float32, size(params)[1])
     ans_valued = Symbolics.value.(substitute(ans, Dict(params[i] => rand_arr[i] for i in 1:length(params))))
     ans_rank = get_rank(ans_valued)

     if ans_rank == min(shape[1], shape[2])
          return true
     else
          return false
     end
end


@variables x1 x2
xdot = [0, exp(2*x2) + exp(x1)*x2]
y = exp(x1) + x2
is_NL_Observable(xdot, [y], [x1, x2], nothing, true)
