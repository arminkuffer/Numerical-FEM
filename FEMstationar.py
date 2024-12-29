import numpy as np
from lagrange2D import linquadderivref,linquadref
from triagplot import quadplot
from integration import gx2dref,gw2dref,getJacobian

#Variablen
Lambda = 48
T1 = 600
T2 = 300
r = 0.02
b = 0.3
h = 0.3


nodes = np.array([[0,0],
                  [b/3,0],
                  [(2/3)*b,0],
                  [b,0],
                   [0,h/3],
                   [b/3,h/3],
                   [(2/3)*b,h/3],
                   [b,h/3],
                   [0,(2*h)/3],
                   [b/3,(2*h)/3],
                   [(2/3)*b,(2*h)/3],
                   [b-r*np.sin(np.pi/6),h-r*np.cos(np.pi/6)],
                   [b,h-r],
                   [b-r*np.cos(np.pi/6),h-r*np.sin(np.pi/6)],
                   [0,h],
                   [b/3,h],
                   [b/2,h],
                   [b-r,h]])
elements = np.array([[1,2,6,5],
                     [2,3,7,6],
                     [3,4,8,7],
                     [5,6,10,9],
                      [6,7,11,10],
                      [7,12,14,11],
                      [7,8,13,12],
                      [9,10,16,15],
                      [10,11,17,16],
                      [11,14,18,17]])
dbc = np.array([[1,T1],
                [2,T1],
                [3,T1],
                [4,T1],
                [12,T2],
                [13,T2],     
                [14,T2],
                [18,T2]])

def evaluate_stat(elenodes,gpx,gpw):
    A = np.zeros((4,4))
    for k in range(len(gpx)):
        for i in range(A.shape[0]):
            for j in range(A.shape[1]):
                print(np.multiply(linquadderivref(gpx[k][0],gpx[k][1])[i],getJacobian(elenodes,gpx[k][0],gpx[k][1])[2]))
                A[i][j] += Lambda*linquadderivref(gpx[k][0],gpx[k][1])[i]*getJacobian(elenodes,gpx[k][0],gpx[k][1])[2]*linquadderivref(gpx[k][0],gpx[k][1])[j]*getJacobian(elenodes,gpx[k][0],gpx[k][1])[2]*getJacobian(elenodes,gpx[k][0],gpx[k][1])[1]*gpw[k]
    return A

print(evaluate_stat(np.array([[0,0],[1,0],[1,2],[0,2]]),gx2dref(3),gw2dref(3)))