using Symbolics
using LinearAlgebra
using Nemo

@variables t x1(t) x2(t) a(t) c(t) d(t)

# Example 1
DX = [exp(x2), x2] ##last 3 eq's are for params "a, c, d"
y1 = x1
y2 = x2



######################## Function for computing Jacobian rank
function get_rank(Matr::Any)
     U, S, VT = svd(Matr)
     count = 0
     for i in S
          if i == 0
               break
          end
          count += 1
     end
     return count
end 
########################


######################## Function for computing ans for exact parameter
function get_ans_exact(Matr::Any, exact::Any, vars::Any)

    random_array = rand(Float32, length(vars))

    _Only_Here_ = copy(Matr)

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
########################



function is_NL_Observable(sys::Any, viewable::Any, params::Any, specific::Any = nothing)
    rand_arr = rand(Float32, length(params))
     n = length(viewable)
     for i in (n+1):(length(sys)-1 + (n*length(params)))
         push!(viewable, dot(Symbolics.gradient(viewable[i - n], params, simplify = true), sys))
     end

    ans = Symbolics.jacobian(viewable, params, simplify=true)
    shape = sizeof(shape(ans))

    if specific != nothing

        exact_ans = get_ans_exact(ans, specific, params)

        return simplify(det(ans)), exact_ans
    end

    ans = get_rank(Symbolics.value.(substitute(ans, Dict(params[i] => rand_arr[i] for i in 1:length(params))))

    return ans == min(shape)
end


println(expand(is_NL_Observable(DX, vec([y1, y2]), [x1, x2])))
