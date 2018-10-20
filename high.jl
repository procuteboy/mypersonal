data = readdlm("x-eb(1)(lor5).DAT",' ')
using QuadGK
using PyPlot
a = data[:,2]
c = data[:,1]
function f2(x)
	return QuadGK.quadgk(p->sqrt(x^2+p^2)*p^2/((x^2+p^2)-350^2),0,602.3)[1]
end
function h2(x)
    pi^2/(8*f2(x))
end
function f1(x)
	return QuadGK.quadgk(p->sqrt(x^2+p^2)/((x^2+p^2)-350^2),0,602.3)[1]
end
function h1(x,y)
   return 1/(1/(2*h2(x))+y/(2*pi^2)*f1(x))
end
b = zeros(15)
d = zeros(15)
e = zeros(15)
f = zeros(15)
for i =1:15
    b[i] = h1(a[i],c[i+30]*300*300)/h2(400)
end
for i =1:15
    e[i] = h2(400)/h2(400)
end
for i =1:15
    f[i] = h1(400,c[i+30]*300*300)/h2(400)
end
for i =1:15
    d[i] = h2(a[i])/h2(400)
end
x_vals = Array{Vector}(2)
y_vals = Array{Vector}(2)
y_vals[1]=b
y_vals[2]=f
x_vals[1]=d
x_vals[2]=e
c =zeros(15)
for i =1:15
    c[i] = 0.1*(i+30)
end
p1=plot(c,y_vals,frame=:box)
