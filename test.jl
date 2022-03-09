using Symbolics
@variables x, y

S = [exp(x) exp(x^2) sin(x)  cos(x)]
s_expr = eval(build_function(S, x)[1])
s_expr(1)


function TO_DO(Matrix::Any, params::Any)
    nums = eval(build_function(Matrix, params)[1])
    gen = rand(Float32, (1, length(params)))
    return nums(gen)    
end
par = [x, y]
TO_DO(S, par)