using Symbolics
using LinearAlgebra


@variables t, x, y



num = (exp(t) + exp(x^2) + cos(y))*sin(x^2) - 10exp(t) + t
expression = Symbolics.toexpr(num)
args = arguments(expression)
ds = Dict{Any, Any}()



function build_evaluation!(dict, expression)
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



function evaluate(expression, dict)
    if typeof(expression) <: Int
        return expression
    end
    if (typeof(expression) == Symbol) || (operation(expression) in [cos, sin, exp])
        return dict[expression]
    end
    new_tree = [evaluate(new_expression, dict) for new_expression in arguments(expression)]
    return Expr(:call, operation(expression), new_tree...)
end
  

build_evaluation!(ds, expression)
evaluate(expression, ds)

ds

num = [(exp(t) + exp(x^2) + cos(y))*sin(x^2) - 10exp(t), exp(t), sin(y)]


expression = Symbolics.toexpr(num)


args = arguments(expression)


ds = Dict{Any, Any}()






jacob = Symbolics.jacobian(num, [t, x, y], simplify = true)
function get_upper_rank(jacobian)
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

