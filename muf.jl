using Plots
pyplot()
data = readdlm("mu1.data",' ')

x_vals = Array{Vector}(3)
x_vals[1] = data[:,2]
c = data[:,1]
data = readdlm("mu2.data",' ')
x_vals[2] = data[:,2]
data = readdlm("mu3.data",' ')
x_vals[3] = data[:,2]
plot(c,x_vals)