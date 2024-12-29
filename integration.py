import numpy as np
import matplotlib.pyplot as plt
import lagrange2D as lg2D
from itertools import product
from lagrange2D import linquadref,linquadderivref
def midpointrule(a,b):
    return (b-a)*f((a+b)/2)

def trapezoidalrule(a,b):
    return 1/2 *(b-a)*(f(a)+f(b))



def gx(n):
    return points[n-1]

def gw(n):
    return weights[n-1]

def f(x):
    return (x**5)/(1+x)**5

def gaussianquadrature1D(x,w,a,b):
    sum = 0
    for i in range(len(x)):
        xi = ((b-a)/2)*x[i]+(a+b)/2
        sum += w[i]*f(xi)
    return (b-a)/2*sum

def gx2dref(n):
    if(n == 1): return np.array([0, 0])
    else:
        pos = list(product(points[n-1],repeat = 2))
        return np.array(pos)
    
def gw2dref(n):
    if n == 1 : return [4.0]
    else: 
        gw = list(product(weights[n-1],repeat=2))
    for i in range(len(gw)):
        gw[i] = gw[i][0]*gw[i][1]
    return gw
        
def getxPos(nodes, xi , eta):
    pos = 0
    for i in range(4):  #4 weil es 4 Ansatzfunktionen für 2D gibt!
        pos += linquadref(xi,eta)[i]*nodes[i]
    return pos 

def getJacobian(nodes,xi,eta):
    jacobival = [0]*3
    jacobian = [0]*2
    for i in range(4):
        jacobian[0] += linquadderivref(xi,eta)[i][0]*nodes[i]
        jacobian[1] += linquadderivref(xi,eta)[i][1]*nodes[i]
    jacobival[0] = jacobian
    jacobival[1] = np.linalg.det(jacobian)
    jacobival[2]= np.linalg.inv(jacobian)
    return jacobival

## FKT 8 noch nicht fertig!

points = [[0],[-1/np.sqrt(3),1/np.sqrt(3)],[float(-np.sqrt(3/5)),0,float(np.sqrt(3/5))]]
weights = [[2],[1,1],[5/9,8/9,5/9]]
print(gaussianquadrature1D(gx(1),gw(1),0,4))
print(gaussianquadrature1D(gx(2),gw(2),0,4))
print(gaussianquadrature1D(gx(3),gw(3),0,4))
print(gx2dref(3))
print(gw2dref(3))
print(getxPos(np.array([[2,1],[4,1],[4,3],[2,2]]),0.577,-0.577))
print(getJacobian(np.array([[2,1],[4,1],[4,3],[2,2]]),0.577,-0.577))