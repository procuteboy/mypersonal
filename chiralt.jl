using Optim
G = 1.926/602.3^2
H = 1.74/602.3^2

k = 12.36/602.3^5
k⁺ = 5.0*k
μ = 10.0
B = 30000.0
T = 100
function f(x::Vector)
    function J(χ)
        function ΩN(χ,p)
            M   = -4*(G-1/2*k*χ)*χ
            E   =  sqrt(p^2 + M^2)
            ϵ0  = sqrt((E+μ)^2 )+sqrt((E-μ)^2 )
            ϵ00 = E+μ
            ϵ01 = E-μ
            ϵ1  = sqrt((E+μ)^2 )+sqrt((E-μ)^2 )
            ϵ10 = E+μ
            ϵ11 = E-μ
            ϵ2  = sqrt((E+μ)^2 )+sqrt((E-μ)^2 )
            ϵ20 = E+μ
            ϵ21 = E-μ
            temperterm = 12*T*(log(1+exp(-ϵ00/T))+log(1+exp(-ϵ01/T)))+2*T*(log(1+exp(-ϵ10/T))+log(1+exp(-ϵ01/T))+log(1+exp(-ϵ20/T))+log(1+exp(-ϵ21/T)))
            return 1/(1+(p^2/602.3^2)^5)*p^2*(6*ϵ0+2*ϵ1+2*ϵ2+temperterm)/(4*pi^2)
        end
        function ΩC(χ,p)
            M   = -4*(G-1/2*k*χ)*χ
            E   =  sqrt(p^2 + M^2)
            ϵc  = sqrt((E+μ)^2)+sqrt((E-μ)^2)
            ϵc0 = E + μ
            ϵc1 = E - μ
            temperterm0 =16*T* (log(1+exp(-ϵc0/T))+log(1+exp(-ϵc1/T)))
            y  =  1/(1+(p^2/602.3^2)^5)*(8*ϵc+temperterm0)/(4*pi^2)*0.5
            for a = 1:1000
                Epa = sqrt(p^2+2*a*B+M^2)
                pa = sqrt(p^2+2*a*B)
                h  = 1/(1+(pa^2/602.3^2)^5)
                ϵca = sqrt((Epa+μ)^2)+sqrt((Epa-μ)^2)
                ϵca0 =Epa+μ
                ϵca1 = Epa-μ
                temperterm1 = 16*T* (log(1+exp(-ϵca0/T))+log(1+exp(-ϵca1/T)))
                y += h*(8*ϵca+temperterm1)/(4*pi^2)
            end
            return y
        end
        return quadgk(p->ΩN(χ,p),0,10000.0)[1]+ B*quadgk(p->ΩC(χ,p),0,10000.0)[1]
end
return -J(x[1])+ 6*G*x[1]^2 -4*k*x[1]^3
end
res = Optim.optimize(f,[-7.5e6], method = ConjugateGradient())
M(χ)   = -4*(G-1/2*k*χ)*χ
χ = Optim.minimizer(res)[1]
print(M(χ))
