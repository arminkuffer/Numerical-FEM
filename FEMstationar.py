import numpy as np
from lagrange2D import linquadderivref
from triagplot import quadplot
from integration import gx2dref,gw2dref,getJacobian

#Variables
Lambda = 48
T1 = 600
T2 = 300
r = 0.02+0.06 #r = 0.08 kleinster Radius für den gilt T15 < 450 (Maximalwert für y = h)
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

#Functions:
def evaluate_stat(elenodes, gpx, gpw):
    elemat = np.zeros((4, 4))
    for k in range(len(gpx)):
        J, detJ, invJ = getJacobian(elenodes, gpx[k][0], gpx[k][1])
        for i in range(elemat.shape[0]):
            for j in range(elemat.shape[1]):
                gradN_i = np.matmul(linquadderivref(gpx[k][0], gpx[k][1])[i], invJ)
                gradN_j = np.matmul(linquadderivref(gpx[k][0], gpx[k][1])[j], invJ)
                elemat[i][j] += Lambda * (gradN_i @ gradN_j) * detJ * gpw[k]   
    return elemat, np.zeros(4)

def assemble(elemat,elevec,sysmat,rhs,ele):
    for i in range(4):
        rhs[int(ele[i]-1)] += elevec[i]
        for j in range(4):	
            sysmat[int(ele[i]-1)][int(ele[j]-1)] += elemat[i][j]
    return sysmat,rhs

def assignDBC(sysmat,rhs,dbc):
    for i in range(dbc.shape[0]):
        rhs[int(dbc[i][0]-1)] = dbc[i][1]
        for k in range(sysmat.shape[1]):
                if k == int(dbc[i][0]-1):
                    sysmat[int(dbc[i][0]-1)][k] = 1
                else:
                    sysmat[int(dbc[i][0]-1)][k] = 0               
    return sysmat,rhs

def solve(nodes,elements,dbc):
    eval = []
    assemb = []
    for e in elements:
        eval = evaluate_stat(nodes[e-1],gx2dref(3),gw2dref(3))
        if(assemb == []):
            assemb = assemble(eval[0],eval[1],np.zeros((nodes.shape[0],nodes.shape[0])),np.zeros((nodes.shape[0],1)),e)
        else:
            assemb = assemble(eval[0],eval[1],assemb[0],assemb[1],e)
    assign = assignDBC(assemb[0],assemb[1],dbc)
    sol = np.linalg.solve(assign[0],assign[1])
    sol = np.array(sol).flatten()
    return sol


sol = solve(nodes,elements,dbc)
print(sol)
quadplot(nodes,elements-1,sol)


