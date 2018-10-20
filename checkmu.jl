using NLsolve

χ = -1.8563798621160723e7
G = 1.926/602.3^2
H = 1.74/602.3^2

k = 12.36/602.3^5
k⁺ = 0.0*k

function f!(x::Vector,fvec)
    function J(x)
        function Ω(x,p)
            M   = -4.0*(G-1/2*k*χ)*χ
            E   =  sqrt(p^2 + M^2)
            return 1/(1+(p^2/602.3^2)^5)*p^2/pi^2*(1/(E+x) + 1/(E-x))
        end
        return quadgk(p->Ω(x,p),0,11602.3)[1]-1/(4*H-k⁺*χ) 
    end
    fvec[1]=J(x[1])
end
nlsolve(f!,[260.1])
