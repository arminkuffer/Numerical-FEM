import numpy as np
import matplotlib.pyplot as plt

def LagrangeBasis(x,n,i,x_node):
    lagrange = 1
    for k in range(n+1):
        if(k != i):
            lagrange = lagrange*((x-x_node[k])/(x_node[i]-x_node[k]))
    return lagrange

def LagrangePolynom(x,n,x_node,f_node):
    poly = 0
    for i in range(n+1):
        poly = poly+np.multiply(f_node[i],LagrangeBasis(x,n,i,x_node))
    return poly

def forwarddifference(h):
    polyh = LagrangePolynom(x+h,o,x_node,f_node)
    poly = LagrangePolynom(x,o,x_node,f_node)
    fx0 = (polyh-poly)/h
    return fx0

def LagrangeBasisDeriv(x,n,i,x_node):
    lagrange = 0
    for m in range(n+1):
        product = 1
        for k in range(n+1):
            if(k!=i and k!=m):
                product = product*(x-x_node[k])/(x_node[i]-x_node[k])
        if(m!=i):
            lagrange += (1/(x_node[i]-x_node[m]))*product 
    return lagrange

def LagrangePolynomDeriv(x,n,x_node,f_node):
    poly = 0
    for i in range(n+1):
        poly += np.multiply(f_node[i],LagrangeBasisDeriv(x,n,i,x_node,f_node))
    return poly

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

##############################################
h = 10**(-6)                                 #stepsize
x_node = [0,1,2,3,4]                         #data points
f_node = (x_node/(np.add(1,x_node)))**5      #data points
x = np.linspace(0,4,500)                     #function input
fx = (x/(np.add(1,x)))**5                    #exact output
o = 4                                        #order  
##############################################
print(LagrangePolynomDeriv(0.6,4,x_node,f_node))
plotlagrange(LagrangePolynom(x,o,x_node,f_node),LagrangePolynomDeriv(x,o,x_node,f_node))
forwarddifference(h)

