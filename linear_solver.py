import numpy as np
import time as t 
def solveGauss(A, b):
    n = len(b)
    x = np.zeros((n,))
    for i in range(n-1):
        for j in range(i+1,n):
            m = (A[j][i] / A[i][i])
            A[j] = A[j] - m*A[i]
            b[j] = b[j] - m*b[i]
    for i in range(n-1,-1,-1):
        x[i] = b[i] / A[i][i] 
        for j in range(i+1,n):
            x[i] -= (A[i][j] * x[j]) / A[i][i]

    return x

def solveG(A,b,x0,rtol,itermax):
    r = b-np.matmul(A,x0)
    k = 0
    x = x0
    while(np.linalg.vector_norm(r,ord=2)>rtol and k < itermax):
        v = np.matmul(A,r)
        alphak = (np.dot(r.T,r))/(np.dot(r.T,v))
        x += alphak*r
        r -= alphak*v
        k+=1
    print(k)
    return x

def solveCG(A,b,x0,rtol,itermax):
    r = b-np.matmul(A,x0)
    p = r
    k = 0
    x = x0
    r_new = np.zeros((len(b),))
    while(np.linalg.vector_norm(r,ord=2)>rtol and k < itermax):
        v = np.matmul(A,p)
        alpha = (np.dot(r.T,r))/(np.dot(p.T,v))
        x += alpha*p
        r_new = r-alpha*v
        beta = (np.dot(r_new.T,r_new))/(np.dot(r.T,r))
        p = r_new + beta*p
        r = r_new
        k+=1
    print(k)
    return x

n = 300
phi = 5.01
b = np.ones((n,))
A = np.diag(np.full((n),phi)) -2*np.diag(np.ones((n-1)),k=1) - 2*np.diag(np.ones((n-1)),k=-1)
time_start = t.time()
sol = solveGauss(A,b)
time_stop = t.time()
print(sol)
print(time_stop-time_start)
time_start = t.time()
sol = solveG(A,b,np.zeros((n,)),0.0000001,1000)
time_stop = t.time()
print(sol)
print(time_stop-time_start)


    
"""print(solveCG(np.array([[10.0, 2.0, 10.0],[ 2.0, 40.0, 8.0],[ 10.0, 8.0, 60.0]]),np.array([1.0, 1.0, 2.0]),np.array([0.0, 0.0, 0.0]), 10**(-7), 1000))
print(solveG(np.array([[10.0, 2.0, 10.0],[ 2.0, 40.0, 8.0],[ 10.0, 8.0, 60.0]]),np.array([1.0, 1.0, 2.0]),np.array([0.0, 0.0, 0.0]), 10**(-7), 1000))
print(solveGauss(np.array([[10.0,2.0,2.0],[3.0,4.0,4.0],[1.0,8.0,4.0]]),np.array([1.0,1.0,2.0])))"""