using Optim
#using Plots
G = 1.9/602.3^2
H = 1.2/602.3^2
k = 12.36/602.3^5
k⁺ = 3.2*k
using QuadGK
μ = 100
function f!(x)
    function J(χ,s)
        function Ω(χ,s,p)
            Δ   = -2*(H-1/4* k⁺*χ)*s
            M   = -4*(G-1/2*k*χ)*χ+1/4* k⁺*s^2
            E   =  sqrt(p^2 + M^2)
            ω₈⁺ = sqrt((E+μ)^2 + 4*Δ^2)
            ω₈⁻ = sqrt((E-μ)^2 + 4*Δ^2)
            ω₁⁺ = sqrt((E+μ)^2 + Δ^2)
            ω₁⁻ = sqrt((E-μ)^2 + Δ^2)
            return p^2*((ω₈⁺+ ω₈⁻)+ 8*(ω₁⁺+ ω₁⁻))
        end
        return QuadGK.quadgk(p->Ω(χ,s,p),0,602.3)[1]
end
return -J(x[1],x[2])/(2*3.14159^2) +6*G*x[1]^2+3*H*x[2]^2-4*k*x[1]^3-1.5*k⁺*x[2]^2*x[1]
end

res = optimize(f!,[-7.5e6,-3.3], method=ConjugateGradient())


print(res,"\n")
print(Δ(Optim.minimizer(res)[1],Optim.minimizer(res)[2]),"\n")
print(M(Optim.minimizer(res)[1],Optim.minimizer(res)[2]),"\n")