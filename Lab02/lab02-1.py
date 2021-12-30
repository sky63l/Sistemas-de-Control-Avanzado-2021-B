import matplotlib.pyplot as plt
import numpy as np
from numpy.core.fromnumeric import size
from scipy.integrate import odeint

to=0
tf=60
t=np.linspace(to,tf,200)

a = 1
k = 1000

def nonlinear(x,t):
    u = (-2 -a * np.cos(x[0])) * (x[0]+a * np.sin(x[0]) + x[1]) - x[0] - k * (x[1] +2 * x[0]+ a* np.sin(x[0]))
    return [x[0] + a * np.sin(x[0]) + x[1], u ]

x0 = [-20,30]

x_sol=odeint(nonlinear,x0,t)
x_dot = x_sol[:,1]    
x = x_sol[:,0]

u_matrix= np.ones(size(t))

for i in range(0,size(t)):
    u = (-2 -a * np.cos(x[i])) * (x[i] +a * np.sin(x[i]) + x_dot[i]) - x[i] - 1000 * (x_dot[i] +2 * x[i]+ a* np.sin(x[i]))
    u_matrix[i]=u_matrix[i]*u
print(u_matrix)

#plt.plot( t, u_matrix, label = "u")
plt.title("Respuesta del Sistema Controlado")
plt.plot( t, x, label = "x_1" )
plt.plot( t, x_dot, label = "x_2" )
plt.legend()
plt.grid()
plt.xlabel("t")
plt.ylabel(" x ")
plt.show()