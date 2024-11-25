import numpy as np
import matplotlib.pyplot as plt
import lagrange2D as lg2D
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


## FKT 5 6 7 8 fehlen!

points = [[0],[-1/np.sqrt(3),1/np.sqrt(3)],[float(-np.sqrt(3/5)),0,float(np.sqrt(3/5))]]
weights = [[2],[1,1],[5/9,8/9,5/9]]

print(gaussianquadrature1D(gx(1),gw(1),0,4))
print(gaussianquadrature1D(gx(2),gw(2),0,4))
print(gaussianquadrature1D(gx(3),gw(3),0,4))