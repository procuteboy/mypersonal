using Optim
G = 1.926/602.3^2
H = 1.74/602.3^2

k = 12.36/602.3^5
k⁺ = 5.0*k
μ = 100.0
B = 100

function f(x::Vector)
    function J(χ,s1,s2)
        function ΩN(χ,s1,s2,p)
            Δ1   = -2*(H-1/4* k⁺*χ)*s1
            Δ2   = -2*(H-1/4* k⁺*χ)*s2
            Δa   = 0.5*( Δ1 +  sqrt(Δ1^2 + 8*Δ2^2))
            Δb   = 0.5*( Δ1 -  sqrt(Δ1^2 + 8*Δ2^2))
            s    = sqrt((s1^2 + 2*s2^2)/3)
            M   = -4*(G-1/2*k*χ)*χ+1/12* k⁺*(s1^2+2*s2^2)
            E   =  sqrt(p^2 + M^2)
            ϵ0  = sqrt((E+μ)^2 + Δ1^2)+sqrt((E-μ)^2 + Δ1^2)
            ϵ1  = sqrt((E+μ)^2 + Δa^2)+sqrt((E-μ)^2 + Δa^2)
            ϵ2  = sqrt((E+μ)^2 + Δb^2)+sqrt((E-μ)^2 + Δb^2)
            return exp(-p^2/602.3^2)*p^2*(6*ϵ0+2*ϵ1+2*ϵ2)/(4*pi^2)
        end
        function ΩC(χ,s1,s2,p)
            Δ1   = -2.0*(H-1/4* k⁺*χ)*s1
            Δ2   = -2.0*(H-1/4* k⁺*χ)*s2
            Δa   = 0.5*( Δ1 +  sqrt(Δ1^2 + 8*Δ2^2))
            Δb   = 0.5*( Δ1 -  sqrt(Δ1^2 + 8*Δ2^2))
            s    = sqrt((s1^2 + 2*s2^2)/3)
            M   = -4*(G-1/2*k*χ)*χ+1/12* k⁺*(s1^2+2*s2^2)
            E   =  sqrt(p^2 + M^2)
            ϵc  = sqrt((E+μ)^2+Δ2^2)+sqrt((E-μ)^2 + Δ2^2)
            y  =  exp(-p^2/602.3^2)*(8*ϵc)/(4*pi^2)*0.5
            for a = 1:10000
                Epa = sqrt(p^2+2*a*B+M^2)
                pa = sqrt(p^2+2*a*B)
                h  = exp(-pa^2/602.3^2)
                ϵca = sqrt((Epa+μ)^2+Δ2^2)+sqrt((Epa-μ)^2 + Δ2^2)
                y += h*(8*ϵca)/(4*pi^2)
            end
            return y
        end
        return quadgk(p->ΩN(χ,s1,s2,p),0,10000.0)[1]+ B*quadgk(p->ΩC(χ,s1,s2,p),0,10000.0)[1]
end
return -J(x[1],x[2],x[3])+ 6*G*x[1]^2+H*(x[2]^2+2*x[3]^2)-4*k*x[1]^3-0.5*k⁺*(x[2]^2+2*x[3]^2)*x[1]
end

Δ1(χ,s1)   = -2.0*(H-1/4* k⁺*χ)*s1
Δ2(χ,s2)   = -2.0*(H-1/4* k⁺*χ)*s2
M(χ,s1,s2)   = -4*(G-1/2*k*χ)*χ+1/12* k⁺*(s1^2+2*s2^2)
res = Optim.optimize(f,[-7.5e6,-0.0,-0.0], method = ConjugateGradient())
print(Δ1(Optim.minimizer(res)[1],Optim.minimizer(res)[2]),"\n")
print(Δ2(Optim.minimizer(res)[1],Optim.minimizer(res)[3]),"\n")
print(M(Optim.minimizer(res)[1],Optim.minimizer(res)[2],Optim.minimizer(res)[3]),"\n")

