import numpy as np
import matplotlib.pyplot as plt

def LagrangeBasis(x,n,i,x_node):
    lagrange = 1
    for k in range(n+1):
        if(k == i):
            continue
        else: 
            lagrange = lagrange*((x-x_node[k])/(x_node[i]-x_node[k]))
    return lagrange

def LagrangePolynom(x,n,x_node,f_node):
    poly = 0
    for i in range(n+1):
        poly = poly+np.multiply(f_node[i],LagrangeBasis(x,n,i,x_node))
    return poly
def forwarddifference(polyh,poly,h):
    fx0 = (polyh-poly)/h
    return fx0
h = 10**(-5)
x_node = [0,1,2,3,4]
f_node = (x_node/(np.add(1,x_node)))**5
x = np.linspace(0,4,500)
fx = (x/(np.add(1,x)))**5
def plotlagrange(polynome,diff):
    plt.plot(x,polynome,color="yellow",linewidth=1)
    plt.plot(x_node,f_node,'o',color='black')
    plt.plot(x,diff,color="green",linewidth="1")
    plt.plot(x,fx,color="black",linewidth=1)
    plt.title('Lagrangian polynome of order 4')
    plt.xlabel('x')
    plt.ylabel('f(x)')
    plt.grid(linewidth=.05)
    plt.show()
plotlagrange(LagrangePolynom(x,4,x_node,f_node),forwarddifference(LagrangePolynom(x+h,4,x_node,f_node),LagrangePolynom(x,4,x_node,f_node),h))

