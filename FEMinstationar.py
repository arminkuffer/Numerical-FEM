import numpy as np
from lagrange2D import linquadref,linquadderivref
from triagplot import quadplot
from integration import gx2dref,gw2dref,getJacobian
import ode as ode
import time as Time
import linear_solver as linsolve


#Variables
rho = 7800
c = 452
Lambda = 48
T0 = np.array([300,300,300,300])
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
def evaluate_instat(elenodes,gpx,gpw,elesol,eleosol,timInt_m,timestep,theta,firststep):
    elemat = np.zeros((4, 4))
    elevec = np.zeros((4,),dtype=object)
    if firststep == 1:
        timInt_m = 1
    for k in range(len(gpx)):
        J, detJ, invJ = getJacobian(elenodes, gpx[k][0], gpx[k][1])
        for i in range(elemat.shape[0]):
            for j in range(elemat.shape[1]):
                N_i = linquadref(gpx[k][0], gpx[k][1])[i]
                N_j = linquadref(gpx[k][0], gpx[k][1])[j]
                gradN_i = np.matmul(linquadderivref(gpx[k][0], gpx[k][1])[i], invJ)
                gradN_j = np.matmul(linquadderivref(gpx[k][0], gpx[k][1])[j], invJ)
                M = rho*c*(N_i * N_j) * detJ * gpw[k]
                B = -(Lambda*(gradN_i @ gradN_j) * detJ * gpw[k])
                C = 0
                if(timInt_m == 1):
                    x = ode.OST(theta,timestep,np.array(M),np.array([B,B]),np.array([0,0]),elesol[j])
                    elemat[i][j] += x[0]
                    elevec[i] += x[1]
                elif(timInt_m == 2):
                    x = ode.AB2(timestep, np.array(M), np.array([B, B]), np.array([C, C]), [elesol[j], eleosol[j]])
                    elemat[i][j] += x[0]
                    elevec[i] += x[1]
                elif(timInt_m == 3):
                    x = ode.AM3(timestep, np.array(M), np.array([B, B, B]), np.array([C, C, C]), [elesol[j], eleosol[j]])
                    elemat[i][j] += x[0]
                    elevec[i] += x[1]
                elif(timInt_m == 4):
                    x = ode.BDF2(timestep,M,B,0,[elesol[j],eleosol[j]])
                    elemat[i][j] += x[0]
                    elevec[i] += x[1]
    return elemat, elevec

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
    elesol = np.full((18, 1), 300)
    eleosol = np.zeros((18,1))
    timInt_m = 1
    deltat = 500
    time = int(5000/deltat)
    theta = 0.5
    firststep = 0
    for t in range(time):
        if t == 0:
            firststep = 1
        else:
            firststep = 0
        eval = []
        assemb = []
        for e in elements:
            if t == 0:
                eval = evaluate_instat(nodes[e-1],gx2dref(2),gw2dref(2),elesol[e-1],eleosol[e-1],timInt_m,deltat,theta,firststep)
            else:
                eval = evaluate_instat(nodes[e-1],gx2dref(2),gw2dref(2),elesol[e-1],eleosol[e-1],timInt_m,deltat,theta,firststep)
            if(assemb == []):
                assemb = assemble(eval[0],eval[1],np.zeros((nodes.shape[0],nodes.shape[0])),np.zeros((nodes.shape[0],1)),e)
            else:
                assemb = assemble(eval[0],eval[1],assemb[0],assemb[1],e)
        assign = assignDBC(assemb[0],assemb[1],dbc)
        T = np.linalg.solve(assign[0],assign[1])
        T = np.array(T).flatten()
        eleosol = elesol
        elesol = T
        if(T[14] > 450):
            print((t+1)*(deltat))
            return T    
    return T

start_time = Time.time()
sol = solve(nodes,elements,dbc)
end_time = Time.time()
elapsed_time = end_time - start_time
print(elapsed_time)
print(sol)
"""print(evaluate_instat(np.array([[0,0],[1,0],[1,2],[0,2]]),gx2dref(3),gw2dref(3),np.array([1,2,3,4]),np.array([0,0,0,0]),1,1000,0.66,1))"""  
quadplot(nodes,elements-1,sol)
