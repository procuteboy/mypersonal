
data = readdlm("mu1.data",' ')
mu1=data[:,2]
data = readdlm("mu2.data",' ')
mu2=data[:,2][1:18]
eB2 = data[:,1][1:18]*300^2/140^2
data = readdlm("mu3.data",' ')
mu3=data[:,2]
eB = data[:,1]*300^2/140^2
eb = linspace(0.0, 0.23, 23)*10^6/140^2
mu4= [315,315,315.2,315.6,316.2,316.7,317.6,318.6,319,320,322.0,324.0,
326.1,328.4,329.2,329.6,329.78,330.0,330.6,331.8,334.2,335.6,337.1]

using PyPlot, LaTeXStrings, PyCall
plot( eB, mu1, color="red", linewidth=2.0, label=L"\mu_c^I")
plot(eB2, mu2, color="green", linewidth=2.0, linestyle="--", label=L"\mu_c^{II}")
plot(eB, mu3, color="blue", linewidth=2.0, linestyle="dotted", label=L"\mu_X")
plot(eb, mu4, color="black", linewidth=2.0, linestyle="dashdot", label=L"\mu_c")
ylabel(L"\mu[MeV]")
xlabel(L"eB/m_\pi^2")
ylim(300,400)
xlim(0.6,20)

subplots_adjust(left=.15, bottom=.16, right=.99, top=.97)
savefig("third.pdf")