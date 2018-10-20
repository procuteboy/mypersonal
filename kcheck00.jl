using NLsolve
using QuadGK
G = 1.926/602.3^2
H = 1.74/602.3^2
χ=-1.3893235264e7
k = 12.36/602.3^5
B = 10^3
M= 355.57
μ = 270.0
function f!(x::Vector,fvec)
    function J(x)
        function Ω(p)
            E   =  sqrt(p^2 + M^2)
            return p^2/pi^2*(1/(E+μ) + 1/(E-μ))
        end
        return QuadGK.quadgk(p->Ω(p),0,602.3)[1]-1/(4*H-x*χ*12.36/602.3^5) 
    end
    fvec[1]=J(x[1])
end
r =nlsolve(f!,[6.1])
