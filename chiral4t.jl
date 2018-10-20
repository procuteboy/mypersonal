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
        function Ω1(x,p)
            M   = -4.0*(G-1/2*k*χ)*χ
            E   =  sqrt(p^2 + M^2)
            f1  = 	1/(exp((E+x)/T)+1)
            f2  = 	1/(exp((E-x)/T)+1)
            return 1/(1+(p^2/602.3^2)^5)*p^2/(2*pi^2)*((1-2*f1)/(E+x) + (1-2*f2)/(E-x))
        end
        function Ω2(x,p)
            M   = -4.0*(G-1/2*k*χ)*χ
            E   =  sqrt(p^2 + M^2)
            f1  = 	1/(exp((E+x)/T)+1)
            f2  = 	1/(exp((E-x)/T)+1)
            y  =  1/(1+(p^2/602.3^2)^5)*((1-2*f1)/(E+x) + (1-2*f2)/(E-x))/(4*pi^2)
            for a = 1:1000
                Epa = sqrt(p^2+2*a*B+M^2)
                pa = sqrt(p^2+2*a*B)
                f3  = 	1/(exp((Epa+x)/T)+1)
                f4  = 	1/(exp((Epa-x)/T)+1)
                y += 1/(1+(pa^2/602.3^2)^5)*((1-2*f3)/(Epa+x) + (1-2*f4)/(Epa-x))/(2*pi^2)
            end
            return y
        end
        return quadgk(p->Ω1(x,p),0,11602.3)[1] + B*quadgk(p->Ω2(x,p),0,11602.3)[1] - 1/(4*H-k⁺*χ) 
    end
    fvec[1]=J(x[1])
end
r=nlsolve(f!,[260.1])
d[i] = r.zero[1]
end
using Plots
scatter(d,e)
