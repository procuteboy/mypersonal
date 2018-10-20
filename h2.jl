using PyPlot, LaTeXStrings, PyCall
PyPlot.rc("font", family="serif", serif="Times")
rcParams = PyDict(PyPlot.matplotlib["rcParams"])
rcParams["text.usetex"] = true
rcParams["xtick.labelsize"] = 8
rcParams["ytick.labelsize"] = 8
rcParams["axes.labelsize"] = 8
figsize(3.487,3.487/1.618)
using QuadGK

data = readdlm("x-eb(1)(lor5).DAT",' ')
M =470
mu =315
a = data[:,2]
c = data[:,1]
function f1(x,B)
	function ω(x,p)
		E =sqrt(p^2+x^2)
	    y=1/(1+(p^2/602.3^2)^5)*(1/(E+mu) + 1/(E-mu))/(pi^2)
	    for a =1:1000
	    	Epa=sqrt(x^2+p^2+2*a*B)
	    	pa = sqrt(p^2+2*a*B)
            y += 1/(1+(pa^2/602.3^2)^5)*(1/(Epa+x) + 1/(Epa-x))/(2*pi^2)
	end
	return y
end
	return B*QuadGK.quadgk(p->ω(x,p),0,10602.3)[1]
end

function f2(x)
	        function Ω1(x,p)
            E   =  sqrt(p^2 + x^2)
            return 1/(1+(p^2/602.3^2)^5)*p^2/(pi^2)*(1/(E+mu) + 1/(E-mu))
        end
	return QuadGK.quadgk(p->Ω1(x,p),0,10602.3)[1]
end
function h2(x)
   return 1/(2*f2(x))
end
function h1(x,B)
   return 1/(f2(x)+f1(x,B) )
end
b = zeros(15)
d = zeros(15)
e = zeros(15)
f = zeros(15)
for i =1:15
    b[i] = h1(a[i],c[i+30]*300*300)/h2(M)
end
for i =1:15
    e[i] = h2(M)/h2(M)
end
for i =1:15
    f[i] = h1(M,c[i+30]*300*300)/h2(M)
end
for i =1:15
    d[i] = h2(a[i])/h2(M)
end
c =zeros(15)
for i =1:15
    c[i] = 0.1*(i+30)*300^2/140^2
end
subplot()
plot(c,b, color="red", linewidth=2.0,)
plot(c,f, color="red", linewidth=2.0, linestyle="--")
ylabel(L"H'_{I}/H'_{0}")
xlabel(L"eB/m_\pi^2")
subplots_adjust(left=.15, bottom=.16, right=.99, top=.97)
savefig("h.pdf")