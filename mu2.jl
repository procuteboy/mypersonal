data = readdlm("x-eb(1)(Lor5).DAT",' ')
using QuadGK
using Plots
pyplot()
a = data[:,2]
c = data[:,1]
b = data[:,3]
using NLsolve

G = 1.926/602.3^2
H = 1.74/602.3^2

k = 12.36/602.3^5
k⁺ = 4.2*k
e = zeros(45)
for i = 1:45
    χ = b[i]

B =  c[i]*300*300
function f!(F,x)
    function J(x)
        function Ω(x,p)
            M   = -4.0*(G-1/2*k*χ)*χ 
            E   =  sqrt(p^2 + M^2)
            return 1/(1+(p^2/602.3^2)^5)*p^2/pi^2*(1/(E+x) + 1/(E-x))
        end
        return QuadGK.quadgk(p->Ω(x,p),0,11602.3)[1]-1/(4*H-k⁺*χ) 
    end
    F[1]=J(x[1])
end
r = nlsolve(f!,[260.1])
e[i] = r.zero[1]
end
plot(c,e,lab="μ")
