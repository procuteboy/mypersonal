
import matplotlib.pyplot as plt
plt.rcParams["font.family"] = "Times New Roman"
SMALL_SIZE = 15
MEDIUM_SIZE = 15
BIGGER_SIZE = 18

          # controls default text sizes
plt.rc('axes', titlesize=16)     # fontsize of the axes title
plt.rc('axes', labelsize=18)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=18)    # fontsize of the tick labels
plt.rc('ytick', labelsize=18)    # fontsize of the tick labels
plt.rc('legend', fontsize=12)    # legend fontsize
#plt.rc('figure', titlesize=16)  # fontsize of the figure title
import numpy as np
from math import log# Simple data to display in various forms
x = np.linspace(1, 20, 50)
y0 = np.log(x)
y1 = 4*np.log(x)
y2 = 4*np.log(x)- 0.01*x*x
plt.figure(figsize=(8,6), dpi=100)
plt.subplot(211)
plt.plot(x, y0,'k')
plt.ylabel(r'$\mathcal{T}_1$')
plt.title('(a)', fontsize=20)
plt.subplot(212)
plt.plot(x, y1, 'k',label =r'$\lambda =0$')
plt.plot(x, y2, 'k--',label = r'$\lambda =0.5$')
plt.title('(b)', fontsize=20)
plt.subplot(212)
plt.xlabel("L/R")
plt.ylabel(r'$\mathcal{T}_2$')
plt.tight_layout()
plt.legend()
plt.show()
