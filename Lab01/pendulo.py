import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import odeint

m=0.2
l=0.5
k=0.1
g=9.8

def pendulum(x,t):
    return [x[1], - g/l * x[0] - k/m * x[1] ]

P0 = [0.3,0]
a=0
b=10
trange=np.linspace(a,b,200)
sol=odeint(pendulum,P0,trange)
y1 = sol[:,1]    
y = sol[:,0]
plt.plot(y, y1)
plt.show()