using Symbolics
using LinearAlgebra
using Nemo

@variables t x1(t) x2(t) a(t) c(t) d(t)
D = Differential(t)
function my_dot_product(X::Any, Y::Any) ## why doesnt built_in work???
    if (size(X) != size(Y))
        return 0;
    end
    ans = 0
    for i in 1:size(Y)[1]
        ans += X[i] * Y[i]
    end
    return ans
end
# Example 1
DX = [a*x1 + x1*x2, c*x2 + d*x1*x2, 0, 0, 0] ##last 3 eq's are for params "a, c, d"
y = x1 
#println(Symbolics.gradient(a*exp(x2), [x1, x2, a], simplify = true))
#= Example 2
DX(t) = [x1(t)^2 + 2x1(t)*x2(t), x2(t)^2]
y(t) = x1(t) + x2(t) 
=#

function is_NL_Observable(sys::Any, viewable::Any, D::Differential, params::Any)

    Y = [viewable]

    for i in 1:size(sys)[1]-1

        Y = vcat(Y, [my_dot_product(Symbolics.gradient(Y[i], params, simplify = true), sys)])

    end

    ans = Symbolics.jacobian(Y, params, simplify=true)

    return simplify(det(ans))

end

print(is_NL_Observable(DX, y, D, [x1, x2, a, c, d]))