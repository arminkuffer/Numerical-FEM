import numpy as np
from lagrange2D import linquadderivref,linquadref
from triagplot import quadplot
from integration import gx2dref,gw2dref,getJacobian
import linear_solver as linsolve
"""def f(x):
    return np.arctan(x)
def derivf(x):
    return 1/(1+x**2)
def newton1D(x0):
    xn = x0+1
    n=0
    while np.abs(f(x0)-f(xn)) >= 10**(-12):
        xn = x0
        x0 = x0-f(x0)/derivf(x0)
        n+=1
    return x0,n

print(newton1D(1.0))"""
Lambda = 48
T1 = 600
T2 = 300
r = 0.02+0.06 #r = 0.08 kleinster Radius für den gilt T15 < 450 (Maximalwert für y = h)
b = 0.3
h = 0.3
c_1 = 10**6
c_2 = 10**3
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
                      [11,7,12,14],
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
T0 = T = np.array([600]*4 + [300]*14)

#Functions:
def evaluate_stat_nlin(elenodes, gpx, gpw,T):
    elemat = np.zeros((4, 4))
    elevec = np.zeros((4,))
    JacobianTk = np.zeros((4,4))
    for k in range(len(gpx)):
        N_jT = 0
        J, detJ, invJ = getJacobian(elenodes, gpx[k][0], gpx[k][1])
        N = linquadref(gpx[k][0],gpx[k][1])
        N_jT = sum(N[j] * T[j] for j in range(elemat.shape[1]))
        for i in range(elemat.shape[0]):
            for j in range(elemat.shape[1]):
                gradN_i = np.matmul(linquadderivref(gpx[k][0], gpx[k][1])[i], invJ)
                gradN_j = np.matmul(linquadderivref(gpx[k][0], gpx[k][1])[j], invJ)
                elemat[i][j] += Lambda * (gradN_i @ gradN_j) * detJ * gpw[k] 
                JacobianTk[i][j] += elemat[i][j] - N[i]*N[j]*(c_1*c_2)/((N_jT)**2) * np.exp(-(c_2)/(N_jT)) * detJ*gpw[k]
            elevec[i] += N[i]*c_1*np.exp(-(c_2/N_jT))*detJ*gpw[k]
    F = -(np.dot(elemat,T)-elevec)
    return JacobianTk, F

def assemble(elemat,elevec,sysmat,rhs,ele):
    for i in range(len(ele)):
        rhs[int(ele[i]-1)] += elevec[i]
        for j in range(len(ele)):	
            sysmat[int(ele[i]-1)][int(ele[j]-1)] += elemat[i][j]
    return sysmat,rhs


def assignDBC_nlin(sysmat,rhs,dbc):
    for i in range(dbc.shape[0]):
        rhs[int(dbc[i][0]-1)] = dbc[i][1]
        for k in range(sysmat.shape[1]):
                if k == int(dbc[i][0]-1):
                    sysmat[int(dbc[i][0]-1)][k] = 1
                else:
                    rhs[k] -= sysmat[k][int(dbc[i][0]-1)]*dbc[i][1]
                    sysmat[int(dbc[i][0]-1)][k] = 0   
                    sysmat[k][int(dbc[i][0]-1)] = 0      
    return sysmat,rhs

def solve(nodes,elements,dbc,T,tol,iter):
    r_k = np.ones((18,))
    k = 0
    while np.linalg.norm(r_k) > tol and k < iter:
        eval = []
        assemb = []
        for e in elements:
            eval = evaluate_stat_nlin(nodes[e-1],gx2dref(3),gw2dref(3),T[e-1])
            if(assemb == []):
                assemb = assemble(eval[0],eval[1],np.zeros((nodes.shape[0],nodes.shape[0])),np.zeros((nodes.shape[0],)),e)
            else:       
                assemb = assemble(eval[0],eval[1],assemb[0],assemb[1],e)
        assign = assignDBC_nlin(assemb[0],assemb[1],dbc)
        r_k = assign[1]
        deltaT = np.linalg.solve(assign[0],assign[1])
        T = T+deltaT
        k+=1
    return T,k

print(evaluate_stat_nlin([[0.0, 0.0], [0.0, 1.0], [1.0, 1.0], [1.0, 0.0]], [[0.5, 0.5]], [1.0], [1.0, 2.0, 2.0, 1.0]))
sol = solve(nodes,elements,dbc,T0,10**(-8),10)
print(sol)