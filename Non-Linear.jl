using Symbolics
using LinearAlgebra
using Nemo

@variables t x1(t) x2(t) a(t) c(t) d(t)

# Example 1
DX = [a*x1 + x1*x2, c*x2 + d*x1*x2, 0, 0, 0] ##last 3 eq's are for params "a, c, d"
y = x1



######################## Function for computing Jacobian rank
function get_rank(Matrix::Any)
     U, S, VT = svd(Matrix)
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
function get_ans_exact(Matrix::Any, exact::Any, vars::Any)
     _Only_Here_ = copy(Matrix)

     before_add = eval(build_function(_Only_Here_, vars)[1])

     rkMat = get_rank(before_add(rand(Int, (1, length(vars)))))

     push!(_Only_Here_, exact)

     after_add = eval(build_function(_Only_Here_, vars)[1])

     rkExMat = get_rank(after_add(rand(Int, (1, length(vars)))))
     if rkExMat == rkMat
          return 1
     else
          return 0
     end
end
########################




#=function test(sys::Any, viewable::Any, params::Any, exact::Any = nothing)
     n = length(viewable)
     for i in (n+1):(length(sys)-1 + n)
         push!(viewable, dot(Symbolics.gradient(viewable[i - n], params, simplify = true), sys))
     end
 
     ans = Symbolics.jacobian(Y, params, simplify=true)
     if exact != nothing
 
         exact_ans = get_ans_exact(ans, exact, params)

         return simplify(det(ans)), exact_ans
     end
     return ans
 end

fff = test(DX, y, [x1, x2, a, c, d])

K = eval(build_function(fff, [x1, x2, a, c, d])[1])
get_rank(K([1, 2, 3, 4, 5])) =#






function is_NL_Observable(sys::Any, viewable::Any, params::Any, exact::Any = nothing)
     n = length(viewable)
     for i in (n+1):(length(sys)-1 + n)
         push!(viewable, dot(Symbolics.gradient(viewable[i - n], params, simplify = true), sys))
     end

    ans = Symbolics.jacobian(Y, params, simplify=true)

    if exact != nothing

        exact_ans = get_ans_exact(ans, exact, params)

        return simplify(det(ans)), exact_ans
    end

    return simplify(det(ans))
end

println(expand(is_NL_Observable(DX, y, [x1, x2, a, c, d])))
