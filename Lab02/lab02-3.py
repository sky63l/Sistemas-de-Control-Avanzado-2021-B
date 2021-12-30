import matplotlib.pyplot as plt
import numpy as np
from numpy.core.fromnumeric import size
from scipy.integrate import odeint

to=0
tf=60
t=np.linspace(to,tf,200)

a = 2.99
k = 1

def nonlinear(x,t):
    u = (-2 -a * np.cos(x[0])) * (x[0]+a * np.sin(x[0]) + x[1]) - x[0] - k * (x[1] +2 * x[0]+ a* np.sin(x[0]))
    return [x[0] + a * np.sin(x[0]) + x[1], u ]

x0 = [-20,30]

x_sol=odeint(nonlinear,x0,t)
x_2 = x_sol[:,1]    
x_1 = x_sol[:,0]

u_matrix= np.ones(size(t))

for i in range(0,size(t)):
    u = (-2 -a * np.cos(x_1[i])) * (x_1[i] +a * np.sin(x_1[i]) + x_2[i]) - x_1[i] - k * (x_2[i] +2 * x_1[i]+ a* np.sin(x_1[i]))
    u_matrix[i]=u_matrix[i]*u



plt.plot( t, u_matrix, label = "u")
plt.plot( t, x_1, label = "x_1" )
plt.plot( t, x_2, label = "x_2" )
plt.legend()
plt.title("Respuesta al Sistema Controlado")
plt.ylabel("x")
plt.xlabel("t")
plt.grid()
plt.show()