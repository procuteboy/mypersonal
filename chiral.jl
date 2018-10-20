using Optim
using QuadGK
G = 1.926/602.3^2
H = 1.74/602.3^2

k = 12.36/602.3^5
k⁺ = 4.2*k
μ = 0.0
#B= 10^5
a = zeros(50)
b = zeros(50)
for i = 1:50
    B =  0.04*i*90000.0
function f(x::Vector)
    function J(χ)
        function ΩN(χ,p)
            M   = -4*(G-1/2*k*χ)*χ
            E   =  sqrt(p^2 + M^2)
            ϵ0  = sqrt((E+μ)^2 )+sqrt((E-μ)^2 )
            ϵ1  = sqrt((E+μ)^2 )+sqrt((E-μ)^2 )
            ϵ2  = sqrt((E+μ)^2 )+sqrt((E-μ)^2 )
            return 1/(1+(p^2/602.3^2)^15)*p^2*(6*ϵ0+2*ϵ1+2*ϵ2)/(4*pi^2)
        end
        function ΩC(χ,p)
            M   = -4*(G-1/2*k*χ)*χ
            E   =  sqrt(p^2 + M^2)
            ϵc  = sqrt((E+μ)^2)+sqrt((E-μ)^2)
            y  =  1/(1+(p^2/602.3^2)^15)*(8*ϵc)/(4*pi^2)*0.5
            for a = 1:1000
                Epa = sqrt(p^2+2*a*B+M^2)
                pa = sqrt(p^2+2*a*B)
                h  = 1/(1+(pa^2/602.3^2)^15)
                ϵca = sqrt((Epa+μ)^2)+sqrt((Epa-μ)^2)
                y += h*(8*ϵca)/(4*pi^2)
            end
            return y
        end
        return QuadGK.quadgk(p->ΩN(χ,p),0,10000.0)[1]+ B*QuadGK.quadgk(p->ΩC(χ,p),0,10000.0)[1]
end
return -J(x[1])+ 6*G*x[1]^2 -4*k*x[1]^3
end
res = Optim.optimize(f,[-7.5e6],BFGS())
M(χ)   = -4*(G-1/2*k*χ)*χ
χ = Optim.minimizer(res)[1]
#print(M(χ))
a[i] = M(χ)
b[i] = χ
end

c = zeros(50)
for i = 1:50
    c[i] = 0.04*i*90000.0
end

