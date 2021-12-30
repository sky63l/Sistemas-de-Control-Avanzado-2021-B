
import numpy as np
from scipy import linalg
import matplotlib.pyplot as plt


class Control():

    @staticmethod
    def modos_deslizantes(dyn,X,d_set,_):
        
        Lambda = 0.3
        k = 2
        
        q_deseado = np.array([[d_set[0]],[d_set[1]]])
        q_deseado_dot = np.array([[d_set[2]],[d_set[3]]])
        q_deseado_ddot = np.array([[d_set[4]],[d_set[5]]])

        q = X[0:2,0]
        q_dot = X[2:4,0]
        
        q_tilda = q-q_deseado
        q_tilda_dot = q_dot-q_deseado_dot

        S = q_tilda_dot + Lambda*q_tilda
        f = -dyn.M_inv @ (dyn.V@q_dot + dyn.G)
        u_hat = q_deseado_ddot - Lambda*q_tilda_dot - f
        print(dyn.G)
        #
        def sat(S):
            d = np.array([[0.1],[0.1]]) 
            if S[0,:] > d[0,:] and S[1,:] > d[1,:] == True:
                return np.sign(S)
            else:
                return S/d
        #

        #u = dyn.M @ (u_hat - k*np.sign(S))
        u = dyn.M @ (u_hat - k*sat(S))
        
        return u,None

class Cinematica(object):
   
    m_1 = 1.0               #masa del primer frame
    m_2 = 1.0               #masa del segundo frame
    L_c1 = 1.0              #longitud del primer frame al centro de masa del primer brazo
    L_c2 = 1.0              #longitud del segundo frame al centro de masa del primer brazo
    L_1 = 2*L_c1            #longitud del primer frame
    L_2 = 2*L_c2            #longitud del segundo frame
    I_1 = (m_1*L_1**2)/3    #momento de inercia del primer frame con respecto al origen
    I_2 = (m_2*L_2**2)/3    #momento de inercia del segundo frame con respecto al primer frame
    G = 9.8                 #constante de gravedad

    def theta_nominal(self):
        
        theta_1 = self.m_1*self.L_c1**2 + self.m_2*(self.L_1**2+self.L_c2**2) + self.I_1 + self.I_2
        theta_2 = self.L_1*self.L_c2
        theta_3 = self.m_2*self.L_c2**2 + self.I_2
        theta_4 = self.m_1*self.L_c1 + self.m_2*self.L_1
        theta_5 = self.m_2*self.L_c2
        return np.array([[theta_1],[theta_2],[theta_3],[theta_4],[theta_5]])

class Dinamica(Cinematica):
    
    # M =  matriz de inercia
    # M_inv = inversa de la matriz de inercia
    # V : matriz de par de Coriolis centrífugo
    # G : matriz sobre pares gravitacionales
    
    def __init__(self,X,dt):
        self.X = X
        self.dt =dt
        self.M = self.matriz_M()
        self.M_inv = self.matriz_M_inv()
        self.V = self.matriz_V()
        self.G = self.matriz_G()

    def matriz_M(self):
        
        q_2 = self.X[1,0]
        d_11 = self.m_1*self.L_c1**2 + self.m_2*(self.L_1**2+self.L_c2**2+2*self.L_1*self.L_c2*np.cos(q_2))+self.I_1+self.I_2
        d_12 = self.m_2*(self.L_c2**2+self.L_1*self.L_c2*np.cos(q_2))+self.I_2
        d_22 = self.m_2*self.L_c2**2+self.I_2
        return np.array([[d_11,d_12],[d_12,d_22]])

    def matriz_M_inv(self):
       
        q_2 = self.X[1,0]
        d_11 = self.m_1*self.L_c1**2 + self.m_2*(self.L_1**2+self.L_c2**2+2*self.L_1*self.L_c2*np.cos(q_2))+self.I_1+self.I_2
        d_12 = self.m_2*(self.L_c2**2+self.L_1*self.L_c2*np.cos(q_2))+self.I_2
        d_22 = self.m_2*self.L_c2**2+self.I_2
        det = d_11*d_22-d_12**2
        inv_11 = d_22/det
        inv_12 = -d_12/det
        inv_22 = d_11/det
        return np.array([[inv_11,inv_12],[inv_12,inv_22]])

    def matriz_V(self):
        
        q_2 = self.X[1,0]
        q_dot_1 = self.X[2,0]
        q_dot_2 = self.X[3,0]
        h = -self.m_2*self.L_1*self.L_c2*np.sin(q_2)
        c_11 = h*q_dot_2
        c_12 = h*q_dot_2 + h*q_dot_1
        c_21 = -h*q_dot_1
        return np.array([[c_11,c_12],[c_21,0]])

    def matriz_G(self):
       
        q_1 = self.X[0,0]
        q_2 = self.X[1,0]
        g_1 = (self.m_1*self.L_c1+self.m_2*self.L_1)*self.G*np.cos(q_1) + self.m_2*self.L_c2*self.G
        g_2 = self.m_2*self.L_c2*self.G*np.cos(q_1+q_2)
        return np.array([[g_1],[g_2]])
    
    def dyn_update(self,tau):
        
        dt = self.dt
        q_ddot = self.M_inv @ (tau - self.V@np.array([[self.X[2,0]],[self.X[3,0]]]) - self.G)
        return self.X + np.mat([[self.X[2,0]*dt + (q_ddot[0,0]*dt**2)/2],\
            [self.X[3,0]*dt + (q_ddot[1,0]*dt**2)/2],[q_ddot[0,0]*dt],[q_ddot[1,0]*dt]])


if __name__ == '__main__':

    f = 1e3
    t_max = 30
    t_lista = np.linspace(0,t_max,int(t_max*f))
    t = t_lista
    dt = 1/f

    X = np.mat([[np.pi/6],[1],[0],[0]])

    q1_lista = [X[0,0]]
    q2_lista = [X[1,0]]

    q1_deseado = [1 for _ in t]
    q2_deseado = [0.5*np.sin(t_i) for t_i in t]
    q1_deseado_dot = [0 for _ in t]
    q2_deseado_dot = [0.5*np.cos(t_i) for t_i in t]
    q1_deseado_ddot = [0 for _ in t]
    q2_deseado_ddot = [-0.5*np.sin(t_i) for t_i in t]
    q_deseado_set = np.array([d_set for d_set in zip(q1_deseado,q2_deseado,q1_deseado_dot,q2_deseado_dot,q1_deseado_ddot,q2_deseado_ddot)])

    def normalizar_angulo(theta):
        return (((theta+np.pi) % (2*np.pi)) - np.pi)  

    q1_deseado_set = [normalizar_angulo(q1_d) for q1_d in q_deseado_set[:,0]]
    q2_deseado_set = [normalizar_angulo(q2_d) for q2_d in q_deseado_set[:,1]]

    q1_e_lista = [q_deseado_set[:,0][0]-X[0,0]]
    q2_e_lista = [q_deseado_set[:,1][0]-X[1,0]]

    control = Control()
    Cinematica = Cinematica()
    theta = Cinematica.theta_nominal()
    #print(theta)

    def plotear():


        plt.figure(1)
        plt.subplot(221)
        line1,=plt.plot(t,q1_lista,'r-')
        line1_d,=plt.plot(t,q1_deseado,'g--')
        plt.legend(handles=(line1,line1_d),labels=("q1","q1 deseado"))
        plt.title("Control de q1 por modos deslizantes")
        plt.xlabel("tiempo[s]")
        plt.ylabel("theta[rad]")
        plt.grid()

        plt.subplot(223)
        line2,=plt.plot(t,q2_lista,'b-')
        line2_d,=plt.plot(t,q2_deseado,'c--')
        plt.legend(handles=(line2,line2_d),labels=("q2","q2 deseado"))
        plt.title("Control de q2 por modos deslizantes")
        plt.xlabel("tiempo[s]")
        plt.ylabel("theta[rad]")
        plt.grid()

        plt.subplot(222)
        plt.plot( t, data[:,1], 'r-', label = "tau_0")
        plt.legend()
        plt.title("Ley de control para q1 con saturación")
        plt.ylabel("theta[rad]")
        plt.xlabel("t[s]")
        plt.grid()

        plt.subplot(224)
        plt.plot( t, data[:,2], label = "tau_1" )
        plt.legend()
        plt.title("Ley de control para q2 con saturación")
        plt.ylabel("theta[rad]")
        plt.xlabel("t[s]")
        plt.grid()
        plt.show()


    tau0_matrix= np.ones(len(t_lista))   
    tau1_matrix= np.ones(len(t_lista))   


    for i in range(len(t_lista)-1):
        dyn = Dinamica(X,dt)
        #print(dyn)
        
        tau,theta = control.modos_deslizantes(dyn,X,q_deseado_set[i+1],theta)
        
        #print(tau)
        tau0_matrix[i]=tau0_matrix[i]*tau[0,:]
        tau1_matrix[i]=tau1_matrix[i]*tau[1,:]
        X = dyn.dyn_update(tau)
        
        q1_lista.append(normalizar_angulo(X[0,0]))
        q2_lista.append(normalizar_angulo(X[1,0]))

        q1_e_lista.append(normalizar_angulo(q_deseado_set[:,0][i+1]-X[0,0]))
        q2_e_lista.append(normalizar_angulo(q_deseado_set[:,1][i+1]-X[1,0]))
  
    #print(tau0_matrix)
    data = np.column_stack([t, tau0_matrix, tau1_matrix])
    np.savetxt("data.txt", data )
    plotear()
