using Optim
G = 1.926/602.3^2
H = 1.74/602.3^2

k = 12.36/602.3^5
k⁺ = 5.0*k
# a = zeros(15)
# b = zeros(15)
# c = zeros(15)
# d = zeros(15)
# for i = 1:15
#     d[i] =i
# end
B = 1.0e4
μ = 550.0
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
            return 1/(1+(p^2/602.3^2)^5)*p^2*(6*ϵ0+2*ϵ1+2*ϵ2)/(4*pi^2)
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
            y  =  1/(1+(p^2/602.3^2)^5)*(8*ϵc)/(4*pi^2)*0.5
            for a = 1:1000
                Epa = sqrt(p^2+2*a*B+M^2)
                pa = sqrt(p^2+2*a*B)
                h  = 1/(1+(pa^2/602.3^2)^5)
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
res = Optim.optimize(f,[-7.5,-3.3e5,-3.3e5], method = ConjugateGradient())
χ = Optim.minimizer(res)[1]
s1 = Optim.minimizer(res)[2]
s2 = Optim.minimizer(res)[3]
a= Δ1(χ,s1)
b= M(χ,s1,s2)
c= Δ2(χ,s2)
