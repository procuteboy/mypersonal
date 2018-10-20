using NLsolve
d = zeros(100)
for i = 1:100
    χ = b[i]
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
r =nlsolve(f!,[260.1])
d[i] = r.zero[1]
end

using LaTeXStrings
plot(c,d,lab="μ")
yaxis!(L"$\mu^{MCFL}_c [MeV]$",(300,500))
xaxis!(L"eB$ [MeV^2]$",(0,3.0e5),0:5.0e4:3.0e5)
