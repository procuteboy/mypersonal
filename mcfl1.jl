using Optim

H = 1.74/602.3^2
μ = 400.0
a=zeros(20)
b=zeros(20)
B =3000.0
function f(x::Vector)
    function J(s1,s2)
        function Ω(s1,s2,p)
            Δ1   = -2*H*s1
            Δ2   = -2*H*s2
            Δa   = 0.5*(   Δ1 +  sqrt(Δ1^2 + 8*Δ2^2))
            Δb   = 0.5*( - Δ1 + sqrt(Δ1^2 + 8*Δ2^2))
            E   = p
            ϵ0  = sqrt((E+μ)^2 + Δ1^2)+sqrt((E-μ)^2 + Δ1^2)
            ϵ1  = sqrt((E+μ)^2 + Δa^2)+sqrt((E-μ)^2 + Δa^2)
            ϵ2  = sqrt((E+μ)^2 + Δb^2)+sqrt((E-μ)^2 + Δb^2)
            return p^2*(6*ϵ0+2*ϵ1+2*ϵ2)/(4*pi^2)
        end
        return quadgk(p->Ω(s1,s2,p),0,602.3)[1]
    end
function sp1(s1,s2,p)
        Δ1   = -2.0*H*s1
        Δ2   = -2.0*H*s2
        F(z) = sqrt(((sqrt(z)+μ))^2+Δ2^2)+sqrt(((sqrt(z)-μ))^2+Δ2^2)-2*sqrt(z+Δ2^2)
        xf = Δ2^2/(2*B)
        return -2*quadgk(y->F(p^2+2*y*B),0.0,10000.0)[1]
end
function sp2(s1,s2,p)
    Δ1   = -2.0*H*s1
    Δ2   = -2.0*H*s2
    F(z) = sqrt(((sqrt(z)+μ))^2+Δ2^2)+sqrt(((sqrt(z)-μ))^2+Δ2^2)-2*sqrt(z+Δ2^2)
        inter =  F(p^2)
        for a = 1:1000 
            inter += 2*F(p^2+2*a*B)
        end
        return inter
end
p1(s1,s2) = (quadgk(p->sp1(s1,s2,p),0,10000.0)[1] + quadgk(p->sp2(s1,s2,p),0,10000.0)[1])*B/(2*pi^2)
function p2(s1,s2)
        Δ1   = -2.0*H*s1
        Δ2   = -2.0*H*s2
        F(z) = sqrt(((sqrt(z)+μ))^2+Δ2^2)+sqrt(((sqrt(z)-μ))^2+Δ2^2)-2*sqrt(z+Δ2^2)
        xf = Δ2^2/(2*B)
        return (zeta(-1,xf) -0.5*(xf^2-xf)*log(xf)+xf^2/4)*B^2/(pi^2)
end
function sp3(s1,s2,p)
    Δ1   = -2.0*H*s1
    Δ2   = -2.0*H*s2
    return (sqrt(((sqrt(p^2)+μ))^2+Δ2^2)+sqrt(((sqrt(p^2)-μ))^2+Δ2^2))*p^2/(pi^2)
end
p3(s1,s2) = quadgk(p->sp3(s1,s2,p),0,602.3)[1]
p(s1,s2) = 4*(p1(s1,s2)+p2(s1,s2)+p3(s1,s2))

return -J(x[1],x[2])+ H*(x[1]^2+2*x[2]^2) - p(x[1],x[2])
end

Δ1(s1)   = -2*H*s1
Δ2(s2)   = -2*H*s2

res = Optim.optimize(f,[-1.3e6,-1.3e6], method = ConjugateGradient())
s1 = Optim.minimizer(res)[1]
s2 = Optim.minimizer(res)[2]
a[1]=Δ1(s1)
b[1]=Δ2(s2)
B=15000.0
res = Optim.optimize(f,[-1.3e6,-1.3e6], method = ConjugateGradient())
s1 = Optim.minimizer(res)[1]
s2 = Optim.minimizer(res)[2]
a[2]=Δ1(s1)
b[2]=Δ2(s2)
B=20000.0
res = Optim.optimize(f,[-1.3e6,-1.3e6], method = ConjugateGradient())
s1 = Optim.minimizer(res)[1]
s2 = Optim.minimizer(res)[2]
a[3]=Δ1(s1)
b[3]=Δ2(s2)
B=25000.0
res = Optim.optimize(f,[-1.3e6,-1.3e6], method = ConjugateGradient())
s1 = Optim.minimizer(res)[1]
s2 = Optim.minimizer(res)[2]
a[4]=Δ1(s1)
b[4]=Δ2(s2)
B=30000.0
res = Optim.optimize(f,[-1.3e6,-1.3e6], method = ConjugateGradient())
s1 = Optim.minimizer(res)[1]
s2 = Optim.minimizer(res)[2]
a[5]=Δ1(s1)
b[5]=Δ2(s2)
B=35000.0
res = Optim.optimize(f,[-1.3e6,-1.3e6], method = ConjugateGradient())
s1 = Optim.minimizer(res)[1]
s2 = Optim.minimizer(res)[2]
a[6]=Δ1(s1)
b[6]=Δ2(s2)
B=40000.0
res = Optim.optimize(f,[-1.3e6,-1.3e6], method = ConjugateGradient())
s1 = Optim.minimizer(res)[1]
s2 = Optim.minimizer(res)[2]
a[7]=Δ1(s1)
b[7]=Δ2(s2)
B=45000.0
res = Optim.optimize(f,[-1.3e6,-1.3e6], method = ConjugateGradient())
s1 = Optim.minimizer(res)[1]
s2 = Optim.minimizer(res)[2]
a[8]=Δ1(s1)
b[8]=Δ2(s2)
B=50000.0
res = Optim.optimize(f,[-1.3e6,-1.3e6], method = ConjugateGradient())
s1 = Optim.minimizer(res)[1]
s2 = Optim.minimizer(res)[2]
a[9]=Δ1(s1)
b[9]=Δ2(s2)
B=55000.0
res = Optim.optimize(f,[-1.3e6,-1.3e6], method = ConjugateGradient())
s1 = Optim.minimizer(res)[1]
s2 = Optim.minimizer(res)[2]
a[10]=Δ1(s1)
b[10]=Δ2(s2)
B=60000.0
res = Optim.optimize(f,[-1.3e6,-1.3e6], method = ConjugateGradient())
s1 = Optim.minimizer(res)[1]
s2 = Optim.minimizer(res)[2]
a[11]=Δ1(s1)
b[11]=Δ2(s2)
B=65000.0
res = Optim.optimize(f,[-1.3e6,-1.3e6], method = ConjugateGradient())
s1 = Optim.minimizer(res)[1]
s2 = Optim.minimizer(res)[2]
a[12]=Δ1(s1)
b[12]=Δ2(s2)
B=70000.0
res = Optim.optimize(f,[-1.3e6,-1.3e6], method = ConjugateGradient())
s1 = Optim.minimizer(res)[1]
s2 = Optim.minimizer(res)[2]
a[13]=Δ1(s1)
b[13]=Δ2(s2)
B=75000.0
res = Optim.optimize(f,[-1.3e6,-1.3e6], method = ConjugateGradient())
s1 = Optim.minimizer(res)[1]
s2 = Optim.minimizer(res)[2]
a[14]=Δ1(s1)
b[14]=Δ2(s2)
B=80000.0
res = Optim.optimize(f,[-1.3e6,-1.3e6], method = ConjugateGradient())
s1 = Optim.minimizer(res)[1]
s2 = Optim.minimizer(res)[2]
a[15]=Δ1(s1)
b[15]=Δ2(s2)
B=85000.0
res = Optim.optimize(f,[-1.3e6,-1.3e6], method = ConjugateGradient())
s1 = Optim.minimizer(res)[1]
s2 = Optim.minimizer(res)[2]
a[16]=Δ1(s1)
b[16]=Δ2(s2)
B=95000.0
res = Optim.optimize(f,[-1.3e6,-1.3e6], method = ConjugateGradient())
s1 = Optim.minimizer(res)[1]
s2 = Optim.minimizer(res)[2]
a[17]=Δ1(s1)
b[17]=Δ2(s2)
B=100000.0
res = Optim.optimize(f,[-1.3e6,-1.3e6], method = ConjugateGradient())
s1 = Optim.minimizer(res)[1]
s2 = Optim.minimizer(res)[2]
a[18]=Δ1(s1)
b[18]=Δ2(s2)
B=105000.0
res = Optim.optimize(f,[-1.3e6,-1.3e6], method = ConjugateGradient())
s1 = Optim.minimizer(res)[1]
s2 = Optim.minimizer(res)[2]
a[19]=Δ1(s1)
b[19]=Δ2(s2)
B=110000.0
res = Optim.optimize(f,[-1.3e6,-1.3e6], method = ConjugateGradient())
s1 = Optim.minimizer(res)[1]
s2 = Optim.minimizer(res)[2]
a[20]=Δ1(s1)
b[20]=Δ2(s2)