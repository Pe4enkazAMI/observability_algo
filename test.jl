using Symbolics
@variables t, x(t), y(t)

function get_ans_exact(Matrix::Any, exact::Any, vars::Any)

    random_array = rand(Float32, length(vars))

    _Only_Here_ = copy(Matrix)

    _Only_Here_ = Symbolics.value.(substitute(_Only_Here_, Dict(vars[i] => random_array[i] for i in 1:length(vars))))

    exact = Symbolics.value.(substitute(exact, Dict(vars[i] => random_array[i] for i in 1:length(vars))))
    
    rkMat = get_rank(_Only_Here_)

    push!(_Only_Here_, exact)

    rkExMat = get_rank(_Only_Here_)
    if rkExMat == rkMat
         return 1
    else
         return 0
    end
end
