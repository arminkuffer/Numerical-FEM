import numpy as np
import matplotlib.pyplot as plt

def Vorwärtseuler(dt):
    phi = np.zeros(len(dt)) # bei zwei endet das Zeitintervall
    phi[0] = 0
    n = 0
    while n<=(len(dt)-2):
        phi[n+1] = phi[n]+dt[1]*fvorwärts(phi[n],dt[n])
        n+=1
    return phi

def Rückwärtseuler(dt):
    phi = np.zeros(len(dt)) # bei zwei endet das Zeitintervall
    phi[0] = 0
    n = 0
    while n<=(len(dt)-2):
        phi[n+1] = (phi[n]+dt[1]*frückwärts(dt[n+1]))/(6*dt[1]+1)
        n+=1
    return phi

def Trapez(dt):
    phi = np.zeros(len(dt)) # bei zwei endet das Zeitintervall
    phi[0] = 0
    n = 0
    while n<=(len(dt)-2):
        phi[n+1] = (phi[n]+dt[1]*ftrapez(dt[n+1])+1/2*dt[1]*fvorwärts(phi[n],t[n]))/(4*dt[1]+1)
        n+=1
    return phi

def fvorwärts(phi,t):
    return (t**2)*np.exp(-5*t)-6*phi

def frückwärts(t):
    return (t**2)*np.exp(-5*t)

def ftrapez(t):
    return (1/2)*(t**2)*np.exp(-5*t)

dt = 0.2
t = np.arange(0,2+dt,dt)
x_ex = np.linspace(0,2,100)
y_ex = np.exp(-5*x_ex)*((x_ex**2)-2*x_ex+2)-2*np.exp(-6*x_ex)
"""plt.plot(x_ex,y_ex)
plt.plot(t,Vorwärtseuler(t))
plt.plot(t,Rückwärtseuler(t), color="red")
plt.plot(t,Trapez(t),color="green")
plt.show()"""

#FKT IX:
def OST(theta,timestep,M,B,C,sol):
    LhsRhs = np.zeros((2,),dtype=object)
    LhsRhs[0] = M-theta*timestep*B[0]
    LhsRhs[1] = np.dot((M + (1 - theta) * timestep * B[1]), sol) + timestep * (theta * C[0] + (1 - theta) * C[1])
    return LhsRhs


#FKT X:
def AB2(timestep,M,B,C,sol):
    LhsRhs = np.zeros((2,),dtype=object)
    LhsRhs[0] = M
    LhsRhs[1] = np.dot(M,sol[0]) + timestep* ((3/2) * (np.dot(B[0],sol[0]) + C[0]) - (1/2) * (np.dot(B[1],sol[1]) -C[1]))
    return LhsRhs

#FKT XI:
def AM3(timestep,M,B,C,sol):
    LhsRhs = np.zeros((2,),dtype=object)    
    LhsRhs[0] = M-timestep*(1/12)*5*B[0]
    LhsRhs[1] = np.dot(M,sol[0])+timestep*(1/12)*(5*C[0]+8*np.dot(B[1],sol[0])+8*C[1]-np.dot(B[2],sol[1])-C[2])
    return LhsRhs

#FKT XII:
def BDF2(timestep,M,B,C,sol):
    LhsRhs = np.zeros((2,),dtype=object)
    LhsRhs[0] = 1.5*M-timestep*B
    LhsRhs[1] = 2*np.dot(M,sol[0])-(1/2)*np.dot(M,sol[1])+timestep*C
    return LhsRhs



#Anwendung auf DGL
def DGLOST(phi0,timestep,theta):
    phi = np.zeros(int(2/timestep+1))
    phi[0] = phi0
    for n in range(1,len(phi)):
        phi[n] = (OST(theta,timestep,[1.0],[-6.0,-6,0],[frückwärts((n)*timestep),frückwärts((n-1)*timestep)],[phi[n-1]])[1])/(OST(theta,timestep,[1.0],[-6.0,-6,0],[frückwärts((n)*timestep),frückwärts((n-1)*timestep)],[phi[n-1]])[0])
    return phi

def DGLAB2(phi0,timestep):
    phi = np.zeros(int(2/timestep+1))
    phi[0] = phi0
    for n in range(1,len(phi)):
        if n == 1: 
            phi[n] = (OST(0.5,timestep,[1.0],[-6.0,-6.0],[frückwärts((n-1)*timestep),frückwärts((n-2)*timestep)],[phi[n-1]])[1])/(OST(0.5,timestep,[1.0],[-6.0,-6,0],[frückwärts((n)*timestep),frückwärts((n-1)*timestep)],[phi[n]])[0])
        else: 
            phi[n] = (AB2(timestep,[1.0],[-6.0,-6.0],[frückwärts((n-1)*timestep),frückwärts((n-2)*timestep)],[phi[n-1],phi[n-2]])[1])/(AB2(timestep,[1.0],[-6.0,-6.0],[frückwärts((n-1)*timestep),frückwärts((n-2)*timestep)],[phi[n-1],phi[n-2]])[0])
    return phi

def DGLAM3(phi0,timestep):
    phi = np.zeros(int(2/timestep+1))
    phi[0] = phi0
    for n in range(1,len(phi)):
        if n == 1: 
            phi[n] = (OST(0.5,timestep,[1.0],[-6.0,-6.0],[frückwärts((n)*timestep),frückwärts((n-1)*timestep)],[phi[n-1]])[1])/(OST(0.5,timestep,[1.0],[-6.0,0],[frückwärts((n)*timestep),frückwärts((n-1)*timestep)],[phi[n-1]])[0])
        else: 
            phi[n] = (AM3(timestep,[1.0],[-6.0,-6.0,-6.0],[frückwärts((n)*timestep),frückwärts((n-1)*timestep),frückwärts((n-2)*timestep)],[phi[n-1],phi[n-2]])[1])/(AM3(timestep,[1.0],[-6.0,-6.0,-6.0],[frückwärts((n)*timestep),frückwärts((n-1)*timestep),frückwärts((n-2)*timestep)],[phi[n-1],phi[n-2]])[0])
    return phi

def DGLBDF2(phi0,timestep):
    phi = np.zeros(int(2/timestep+1))
    phi[0] = phi0
    for n in range(1,len(phi)):
        if n == 1: 
            phi[n] = (OST(0.5,timestep,[1.0],[-6.0,-6.0],[frückwärts((n)*timestep),frückwärts((n-1)*timestep)],[phi[n-1]])[1])/(OST(0.5,timestep,[1.0],[-6.0,0],[frückwärts((n)*timestep),frückwärts((n-1)*timestep)],[phi[n-1]])[0])
        else: 
            phi[n] = (BDF2(timestep,[1.0],[-6.0],[frückwärts((n)*timestep)],[phi[n-1],phi[n-2]])[1])/(BDF2(timestep,[1.0],[-6.0],[frückwärts((n)*timestep)],[phi[n-1],phi[n-2]])[0])
    return phi

"""plt.plot(t,Vorwärtseuler(t),color="orange")
plt.plot(t,DGLOST(0,dt,0.5),color="red")
plt.plot(t,DGLAB2(0,dt),color="blue")
plt.plot(t,DGLAM3(0,dt),color="green")
plt.plot(t,DGLBDF2(0,dt),color="yellow")
plt.plot(x_ex,y_ex,color="black")
plt.show()"""
#print(OST(0.5, 0.2, np.array([1.1]), np.array([1.4, 1.5]), np.array([1.7, 1.8]), np.array([2.0])))
# print(AB2(0.2, np.array([1.1]), np.array([1.5, 1.6]), np.array([1.8, 1.9]), np.array([2.0, 2.1])))
"""print(AM3(0.2, np.array([1.1]), np.array([1.4, 1.5, 1.6]), np.array([1.7, 1.8, 1.9]), np.array([2.0, 2.1])))
print(BDF2(0.2, np.array([1.1]), np.array([1.4]), np.array([1.7]), np.array([2.0, 2.1])))"""
"""M = np.array([[1.1, 1.2],[1.2, 1.1]])
B = np.array([[1.4, 1.5],[1.5 ,1.4]])
C = np.array([1.7,1.8])
sol = np.array([2.0,3.0])"""
# print(OST(0.5, 0.2, M, np.array([B, B]), np.array([C, C]), sol))
"""print(AB2(0.2, M, [B, B], [C, C], [sol,sol]))
print(AM3(0.2, M, [B, B, B], [C, C, C], [sol,sol]))
print(BDF2(0.2, M, B, C, [sol,sol]))"""
