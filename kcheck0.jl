using NLsolve
G = 1.926/602.3^2
H = 1.74/602.3^2
χ=-1.78846e7
k = 12.36/602.3^5
B = 10^3
M= 479.57
μ = 370.0
function f!(x::Vector,fvec)
    function J(x)
        function Ω(p)
            E   =  sqrt(p^2 + M^2)
            return 1/(1+(p^2/602.3^2)^5)*p^2/pi^2*(1/(E+μ) + 1/(E-μ))
        end
        return quadgk(p->Ω(p),0,11602.3)[1]-1/(4*H-x*χ*12.36/602.3^5) 
    end
    fvec[1]=J(x[1])
end
r =nlsolve(f!,[3.1])
