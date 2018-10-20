using NLsolve

G = 1.926/602.3^2
H = 1.74/602.3^2

k = 12.36/602.3^5
k⁺ = 4.2*k
e = zeros(50)
for i = 1:50
    χ = b[i]
B =10000.0
B +=  i*10000.0
function f!(x::Vector,fvec)
    function J(x)
        function Ω1(x,p)
            M   = -4.0*(G-1/2*k*χ)*χ
            E   =  sqrt(p^2 + M^2)
            return 1/(1+(p^2/602.3^2)^5)*p^2/(2*pi^2)*(1/(E+x) + 1/(E-x))
        end
        function Ω2(x,p)
            M   = -4.0*(G-1/2*k*χ)*χ
            E   =  sqrt(p^2 + M^2)
            y  =  1/(1+(p^2/602.3^2)^5)*(1/(E+x) + 1/(E-x))/(4*pi^2)
            for a = 1:10000
                Epa = sqrt(p^2+2*a*B+M^2)
                pa = sqrt(p^2+2*a*B) 
                y += 1/(1+(pa^2/602.3^2)^5)*(1/(Epa+x) + 1/(Epa-x))/(2*pi^2)
            end
            return y
        end
        return quadgk(p->Ω1(x,p),0,11602.3)[1] + B*quadgk(p->Ω2(x,p),0,11602.3)[1] - 1/(4*H-k⁺*χ) 
    end
    fvec[1]=J(x[1])
end
r = nlsolve(f!,[260.1])
e[i] = r.zero[1]
end
plot(c,e,lab="μ")
yaxis!(L"$\mu^{SMCFL}_c [MeV]$",(300,500))
xaxis!(L"eB$ [MeV^2]$",(0,3.0e5),0:5.0e4:3.0e5)