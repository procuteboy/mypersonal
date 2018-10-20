
H = 1.74/602.3^2
μ = 400.0
B = 1000.0
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
p(s1,s2) = 2*(p1(s1,s2)+p2(s1,s2)+p3(s1,s2))
p(-1.3e6,-1.3e6)


H = 1.74/602.3^2
μ = 400.0
B = 1000.0

 function ΩC(s1,s2,p)
            Δ1   = -2.0*H*s1
            Δ2   = -2.0*H*s2

            E   =  p
            ϵc  = sqrt((E+μ)^2+Δ2^2)+sqrt((E-μ)^2 + Δ2^2)
            y  = exp(-p^2/602.3^2)*(8*ϵc)/(8*pi^2)
            for a = 1:1000
                pa = sqrt(p^2+2*a*B)
                Epa = sqrt(p^2 + 2*a*B)
#               h  = 1/(1+exp(20*(Epa/602.3-1)))
                h  = 1/(1+(pa^2/602.3^2)^5)
                ϵca = sqrt((Epa+μ)^2+Δ2^2)+sqrt((Epa-μ)^2 + Δ2^2)
                y += h*(8*ϵca)/(4*pi^2)
            end
            return y
        end
Q(s1,s2) = B*quadgk(p->ΩC(s1,s2,p),0.0,10000.0)[1]


Q(-1.3e6,-1.3e6)
H = 1.74/602.3^2
μ = 400.0
B = 1000.0

        function ΩN(s1,s2,p)
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
Z(s1,s2)=quadgk(p->ΩN(s1,s2,p),0,602.3)[1]






