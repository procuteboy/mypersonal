using NLsolve
G = 1.9/602.3^2
H = 1.74/602.3^2
k = 12.36/602.3^5
k⁺ = 5.0*k
d =zeros(100)
e = zeros(100)
for i = 1:100
    e[i] = 50+0.5*i
end
χ = 0.0
B = 3000.0
for i = 1:100
    T = 50+0.5*i
function f!(x::Vector,fvec)
    function J(x)
        function Ω(x,p)
            M   = -4.0*(G-1/2*k*χ)*χ
            E   =  sqrt(p^2 + M^2)
            f1  = 	1/(exp((E+x)/T)+1)
            f2  = 	1/(exp((E-x)/T)+1)
            return 1/(1+(p^2/602.3^2)^5)*p^2/pi^2*((1-2*f1)/(E+x) + (1-2*f2)/(E-x))
        end
        return quadgk(p->Ω(x,p),0,11602.3)[1]-1/(4*H-k⁺*χ) 
    end
    fvec[1]=J(x[1])
end
r=nlsolve(f!,[260.1])
d[i] = r.zero[1]
end
using Plots
scatter(d,e)
